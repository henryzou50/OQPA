"""Werner state via a Stinespring dilation of the 2-qubit depolarising channel.

We prepare the 4-qubit pure state on ``S`` (system, qubits 0, 1) and
``E`` (environment, qubits 2, 3)

    |Psi>_SE = sqrt((1 + 3p)/4) |00>_S|00>_E
             + sqrt((1 - p)/4) |01>_S|01>_E
             + sqrt((1 - p)/4) |10>_S|10>_E
             + sqrt((1 - p)/4) |11>_S|11>_E ,

so that ``rho_S = Tr_E |Psi><Psi|`` is the Werner state

    rho_S = p |00><00| + (1 - p) I / 4 .

The dilation is realised by preparing ``S`` in the Schmidt vector
``(sqrt((1+3p)/4), sqrt((1-p)/4), sqrt((1-p)/4), sqrt((1-p)/4))`` and
then copying each system qubit onto its environment partner with a
CNOT.

Two implementations of the Schmidt preparation on ``S`` are provided:

* ``"direct"``   -- explicit ``RY``/``CNOT`` gates derived analytically
  from the Schmidt decomposition (3 ``RY`` + 2 ``CNOT``); no use of
  ``StatePreparation`` so the circuit needs no high-level transpilation
  pass to run on AER.
* ``"stateprep"``-- qiskit's high-level ``StatePreparation`` instruction.

Run::

    python werner_stinespring.py                          # default p sweep + table
    python werner_stinespring.py -p 0.3 0.6 0.9           # custom p values
    python werner_stinespring.py --draw                   # also print the circuit
    python werner_stinespring.py --method stateprep       # use StatePreparation
    python werner_stinespring.py --plot                   # F(rho_S, |00>) vs p plot

Requires ``qiskit >= 1.0``, ``qiskit-aer`` and (for ``--plot``)
``matplotlib``.
"""

from __future__ import annotations

import argparse
import sys
from typing import Sequence

import numpy as np

from qiskit import QuantumCircuit, QuantumRegister, transpile
from qiskit.circuit.library import StatePreparation
from qiskit.quantum_info import DensityMatrix, Statevector, state_fidelity
from qiskit_aer import AerSimulator

K = 2
DIM = 2 ** K
SYS_QUBITS = list(range(K))


METHODS = ("direct", "stateprep")


def system_schmidt_vector(p: float) -> np.ndarray:
    """Schmidt coefficients of |Psi>_SE on the 2-qubit system register."""
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"p must lie in [0, 1] (got {p}).")
    a = np.sqrt((1.0 + 3.0 * p) / 4.0)
    b = np.sqrt((1.0 - p) / 4.0)
    return np.array([a, b, b, b], dtype=complex)


def system_schmidt_preparation_direct(p: float) -> QuantumCircuit:
    """Two-qubit subcircuit preparing the Schmidt vector on ``S`` from |00>.

    The target state on ``S`` is

        |phi> = sqrt((1+3p)/4) |00>
              + sqrt((1-p)/4) (|01> + |10> + |11>) .

    Schmidt by qubit (qubit 1 = s_1 = MSB in qiskit, qubit 0 = s_0 = LSB):

        |phi> = sqrt((1+p)/2) |0>_{s1} ( sqrt((1+3p)/(2(1+p))) |0>_{s0}
                                       + sqrt((1-p)/(2(1+p))) |1>_{s0} )
              + sqrt((1-p)/2) |1>_{s1} ( |+>_{s0} ) ,

    which we build in two stages:

      1. ``RY(arccos p)`` on s_1 sets |0>_{s1} amplitude = sqrt((1+p)/2)
         and |1>_{s1} amplitude = sqrt((1-p)/2).
      2. A uniformly-controlled ``RY`` on s_0 (control = s_1) with angles
            theta_0 = arccos(2p/(1+p))    if s_1 = 0
            theta_1 = pi/2                if s_1 = 1
         decomposes into two CNOTs and two RY's as

            RY(alpha)  CNOT(s_1 -> s_0)  RY(beta)  CNOT(s_1 -> s_0)

         with ``alpha = (theta_0 + theta_1)/2`` and
         ``beta = (theta_0 - theta_1)/2``.

    Total gate count: 3 ``RY`` + 2 ``CNOT``.
    """
    if not 0.0 <= p <= 1.0:
        raise ValueError(f"p must lie in [0, 1] (got {p}).")

    sys_reg = QuantumRegister(K, name="s")
    qc = QuantumCircuit(sys_reg, name="schmidt_direct")

    theta_top = float(np.arccos(p))                          # rotation on s_1
    theta_if0 = float(np.arccos((2.0 * p) / (1.0 + p)))      # s_1 = 0 branch
    theta_if1 = np.pi / 2.0                                  # s_1 = 1 branch
    alpha = 0.5 * (theta_if0 + theta_if1)
    beta = 0.5 * (theta_if0 - theta_if1)

    qc.ry(theta_top, sys_reg[1])
    qc.ry(alpha, sys_reg[0])
    qc.cx(sys_reg[1], sys_reg[0])
    qc.ry(beta, sys_reg[0])
    qc.cx(sys_reg[1], sys_reg[0])
    return qc


def werner_dilation_circuit(p: float, *, method: str = "direct") -> QuantumCircuit:
    """Build the 4-qubit Stinespring-dilation circuit.

    Convention: qubits ``[0, 1]`` are the system ``S``, qubits
    ``[2, 3]`` are the environment ``E`` with ``e_i`` paired to ``s_i``.

    Parameters
    ----------
    p : float
        Werner-state parameter in ``[0, 1]``.  ``p = 1`` is noiseless,
        ``p = 0`` is fully depolarised.
    method : {"direct", "stateprep"}, default ``"direct"``
        How the Schmidt vector on ``S`` is prepared.  ``"direct"`` uses
        only ``RY``/``CNOT`` gates derived analytically;
        ``"stateprep"`` uses qiskit's ``StatePreparation`` instruction.
        Both yield the same |Psi>_SE up to global phase.
    """
    if method not in METHODS:
        raise ValueError(
            f"unknown method {method!r}; expected one of {METHODS}."
        )

    sys_reg = QuantumRegister(K, name="s")
    env_reg = QuantumRegister(K, name="e")
    qc = QuantumCircuit(sys_reg, env_reg, name=f"werner_dilation_{method}")

    if method == "direct":
        qc.compose(
            system_schmidt_preparation_direct(p),
            qubits=sys_reg[:],
            inplace=True,
        )
    else:  # "stateprep"
        qc.append(StatePreparation(system_schmidt_vector(p)), sys_reg[:])

    for i in range(K):
        qc.cx(sys_reg[i], env_reg[i])
    return qc


def analytic_werner_state(p: float) -> DensityMatrix:
    """``rho = p |00><00| + (1 - p) I / 4``."""
    rho = ((1.0 - p) / DIM) * np.eye(DIM, dtype=complex)
    rho[0, 0] += p
    return DensityMatrix(rho)


def simulate_reduced_state(
    p: float,
    *,
    method: str = "direct",
    backend: AerSimulator | None = None,
) -> tuple[QuantumCircuit, DensityMatrix, DensityMatrix]:
    """Run the dilation on AER and return ``(circuit, rho_aer, rho_exact)``."""
    qc = werner_dilation_circuit(p, method=method)

    sim_qc = qc.copy()
    sim_qc.save_density_matrix(qubits=SYS_QUBITS, label="rho_S")

    if backend is None:
        backend = AerSimulator(method="statevector")
    # The "direct" method emits only RY/CNOT and runs on AER as-is; the
    # "stateprep" method emits a high-level StatePreparation instruction
    # that has to be transpiled to AER's basis first.  Transpiling
    # unconditionally is a small, safe overhead.
    sim_qc = transpile(sim_qc, backend=backend, optimization_level=1)
    result = backend.run(sim_qc, shots=1).result()
    rho_aer = DensityMatrix(result.data(0)["rho_S"])
    rho_exact = analytic_werner_state(p)
    return qc, rho_aer, rho_exact


def _format_complex_matrix(mat: np.ndarray, *, prec: int = 4) -> str:
    fmt = f"{{:+.{prec}f}}"
    rows = []
    for row in mat:
        cells = []
        for z in row:
            re = fmt.format(float(np.real(z)))
            im = fmt.format(float(np.imag(z)))
            cells.append(f"{re}{im}j")
        rows.append("  ".join(cells))
    return "\n".join(rows)


def _print_circuit(circ: QuantumCircuit) -> None:
    """Print the circuit, falling back to ASCII on Windows code pages."""
    try:
        print(circ.draw(output="text"))
    except UnicodeEncodeError:
        ascii_drawing = circ.draw(output="text").single_string().encode(
            sys.stdout.encoding or "ascii", errors="replace"
        )
        print(ascii_drawing.decode(sys.stdout.encoding or "ascii", errors="replace"))


def run_table(
    ps: Sequence[float], *, method: str, draw: bool, show_rho: bool
) -> int:
    """Per-``p`` table of ``||rho_AER - rho_exact||_F`` and purity."""
    backend = AerSimulator(method="statevector")

    if draw:
        circ = werner_dilation_circuit(ps[0], method=method)
        print(
            f"Dilation circuit (method={method!r}) for p={ps[0]} "
            "(qubits 0,1 = system S, 2,3 = environment E):"
        )
        _print_circuit(circ)
        print()

    print(
        f"k = {K} system qubits  ->  {2 * K}-qubit dilation "
        f"(Schmidt-preparation method = {method!r}).\n"
    )
    header = f"{'p':>6}  {'||rho_AER - rho_exact||_F':>26}  {'purity(rho_S)':>14}"
    print(header)
    print("-" * len(header))

    worst = 0.0
    for p in ps:
        _, rho_aer, rho_exact = simulate_reduced_state(
            p, method=method, backend=backend
        )
        diff = float(np.linalg.norm(rho_aer.data - rho_exact.data))
        purity = float(np.real(np.trace(rho_aer.data @ rho_aer.data)))
        worst = max(worst, diff)
        print(f"{p:6.3f}  {diff:26.3e}  {purity:14.6f}")

        if show_rho:
            print("  rho_S (from AER):")
            print("    " + _format_complex_matrix(rho_aer.data).replace("\n", "\n    "))
            print()

    print()
    if worst < 1e-9:
        print(f"OK: max ||rho_AER - rho_exact||_F over the sweep = {worst:.2e}.")
        return 0
    print(
        f"WARNING: max ||rho_AER - rho_exact||_F = {worst:.2e} "
        "is larger than expected (1e-9).",
        file=sys.stderr,
    )
    return 1


def run_fidelity_plot(
    n_points: int = 21,
    *,
    method: str = "direct",
    save_path: str | None = None,
    show: bool = True,
) -> int:
    """Sweep ``p`` over ``n_points`` from 0 to 1, plot ``F(rho_S, |00>)``."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(
            "matplotlib is required for --plot. Install with `pip install matplotlib`.",
            file=sys.stderr,
        )
        return 1

    backend = AerSimulator(method="statevector")
    ket_00 = Statevector.from_label("00")

    ps = np.linspace(0.0, 1.0, n_points)
    fid_aer = np.empty(n_points)
    for idx, p in enumerate(ps):
        _, rho_aer, _ = simulate_reduced_state(
            float(p), method=method, backend=backend
        )
        fid_aer[idx] = float(state_fidelity(rho_aer, ket_00))

    # Analytic curve: F(rho, |00>) = <00|rho|00> = (1 + 3p) / 4.
    p_fine = np.linspace(0.0, 1.0, 401)
    fid_analytic = (1.0 + 3.0 * p_fine) / 4.0

    print(f"{'p':>6}  {'F(rho_S, |00>) AER':>22}  {'F analytic (1+3p)/4':>22}")
    print("-" * 56)
    for p, f in zip(ps, fid_aer):
        f_th = (1.0 + 3.0 * float(p)) / 4.0
        print(f"{float(p):6.3f}  {f:22.6f}  {f_th:22.6f}")

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.plot(
        p_fine,
        fid_analytic,
        color="C0",
        linewidth=2,
        label=r"analytic  $F=(1+3p)/4$",
    )
    ax.plot(
        ps,
        fid_aer,
        marker="o",
        linestyle="",
        color="C3",
        markersize=6,
        label="AER simulation",
    )
    ax.set_xlabel(r"depolarising parameter $p$")
    ax.set_ylabel(r"fidelity  $F(\rho_S,\,|00\rangle)$")
    ax.set_title(
        r"Werner-state fidelity with $|00\rangle$ vs. $p$ "
        rf"(Stinespring dilation, $k=2$, method={method!r})"
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.2, 1.02)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="lower right")
    fig.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150)
        print(f"\nSaved plot to {save_path}")
    if show:
        plt.show()
    plt.close(fig)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    # Best-effort: switch stdout to UTF-8 so qiskit's box-drawing
    # characters render on Windows consoles (cp1252).
    try:
        sys.stdout.reconfigure(encoding="utf-8")  # type: ignore[attr-defined]
    except (AttributeError, OSError):
        pass

    parser = argparse.ArgumentParser(
        description=(
            "Simulate a 2-qubit Werner state by injecting depolarising noise "
            "through a Stinespring dilation, then verify on Qiskit Aer."
        ),
    )
    parser.add_argument(
        "-p",
        type=float,
        nargs="+",
        default=[0.0, 0.25, 0.5, 0.75, 1.0],
        help=(
            "depolarising parameter(s) for the table mode; p=1 => noiseless, "
            "p=0 => fully depolarised (default: 0 0.25 0.5 0.75 1)"
        ),
    )
    parser.add_argument(
        "--method",
        choices=METHODS,
        default="direct",
        help=(
            "how the Schmidt vector on S is prepared: 'direct' uses an "
            "explicit RY/CNOT decomposition (default); 'stateprep' uses "
            "qiskit's StatePreparation high-level instruction"
        ),
    )
    parser.add_argument(
        "--draw",
        action="store_true",
        help="print the dilation circuit for the first p value",
    )
    parser.add_argument(
        "--show-rho",
        action="store_true",
        help="print the AER reduced density matrix for each p",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help=(
            "sweep 21 points of p in [0, 1] and plot F(rho_S, |00>) "
            "vs p instead of running the table mode"
        ),
    )
    parser.add_argument(
        "--save-plot",
        type=str,
        default=None,
        metavar="PATH",
        help="when used with --plot, also save the figure to PATH (e.g. fidelity.png)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="when used with --plot, do not open an interactive window",
    )
    args = parser.parse_args(argv)

    if args.plot:
        return run_fidelity_plot(
            n_points=21,
            method=args.method,
            save_path=args.save_plot,
            show=not args.no_show,
        )
    return run_table(
        args.p, method=args.method, draw=args.draw, show_rho=args.show_rho
    )


if __name__ == "__main__":
    raise SystemExit(main())
