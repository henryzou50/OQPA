"""
QPA Circuit Engine
==================
Builds Quantum Phase Amplification circuits using two strategies:

  HybridStrategy    - single dynamic circuit with classical feedback (if_test)
  UnrolledStrategy  - set of static circuits, one per execution path

Both strategies implement the same QPA protocol:
  1. Prepare n_registers (odd, 2i+1) of k qubits each
  2. For each trial: perform parallel Schur tests on register pairs
  3. Rotate surviving registers via cyclic permutation
  4. Measure the reserve register at the end
"""

from itertools import product
from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

import numpy as np

from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.circuit import Parameter
from qiskit.circuit.library import RZGate, RXGate
from qiskit_aer.noise import depolarizing_error


# ==========================================
# CORE OPERATIONS
# ==========================================

def apply_schur_test(qc, ancilla, reg_a, reg_b, k, include_measure=False, res_bit=None):
    """
    Schur (swap) test: H - CSWAP^k - H between two k-qubit registers.

    When include_measure=True, also measures into res_bit and resets
    the ancilla afterward (used by HybridStrategy).
    """
    qc.h(ancilla)
    for i in range(k):
        qc.cswap(ancilla, reg_a[i], reg_b[i])
    qc.h(ancilla)
    if include_measure:
        qc.measure(ancilla, res_bit)
        qc.reset(ancilla)
        return res_bit


def apply_cyclic_rotation(qc, regs, k):
    """Cyclic rotation of data across QuantumRegister objects using SWAPs."""
    N = len(regs)
    for i in range(k):
        for j in range(N - 1, 0, -1):
            qc.swap(regs[j][i], regs[j - 1][i])


def apply_cyclic_rotation_indices(qc, indices, k):
    """
    Cyclic rotation using register indices.
    Assumes data registers are the first N*k qubits in the circuit.
    """
    for i in range(len(indices) - 1, 0, -1):
        idx_j = indices[i]
        idx_prev = indices[i - 1]
        for b in range(k):
            q_j = qc.qubits[idx_j * k + b]
            q_prev = qc.qubits[idx_prev * k + b]
            qc.swap(q_j, q_prev)


# ==========================================
# REGISTER MANAGEMENT
# ==========================================

class QPARegisters:
    """Manages quantum and classical registers for the QPA protocol."""

    def __init__(self, n_registers, k, n_trials, use_ancilla_pool=False, no_reset=False,
                 match_dynamic=False):
        self.n = n_registers
        self.k = k
        self.n_trials = n_trials
        self.no_reset = no_reset
        self.match_dynamic = match_dynamic

        self.qr_data = [QuantumRegister(k, f"R{i+1}") for i in range(n_registers)]
        self.qr_ancilla = []
        self.cr_ancilla = None
        self.cr_readout = ClassicalRegister(k, "readout")
        self.use_ancilla_pool = use_ancilla_pool

        # Extra registers only used when match_dynamic=True
        self.qr_ctrl_par = None
        self.qr_ctrl_rec = None
        self.cr_pool = None
        self.cr_rec = None

        if use_ancilla_pool:
            num_concurrent = n_registers // 2

            if match_dynamic:
                # Mirror HybridStrategy register structure
                self.qr_ctrl_par = QuantumRegister(num_concurrent, "ctrl_par")
                self.qr_ctrl_rec = QuantumRegister(1, "ctrl_rec")
                self.qr_ancilla = [self.qr_ctrl_par[i] for i in range(num_concurrent)]
                self.cr_pool = [ClassicalRegister(num_concurrent, f"res_t{t}")
                                for t in range(n_trials)]
                self.cr_rec = ClassicalRegister(n_trials, "res_rec")
                # Build flat clbit list so cr_ancilla[cl_idx] still works
                self._cr_ancilla_flat = []
                for t in range(n_trials):
                    for i in range(num_concurrent):
                        self._cr_ancilla_flat.append(self.cr_pool[t][i])
                self.cr_ancilla = self._cr_ancilla_flat
            else:
                if no_reset:
                    total = num_concurrent * n_trials
                    self.qr_ancilla = [QuantumRegister(1, f"anc_{i}") for i in range(total)]
                else:
                    self.qr_ancilla = [QuantumRegister(1, f"anc_{i}") for i in range(num_concurrent)]
                max_measurements = n_trials * num_concurrent
                self.cr_ancilla = ClassicalRegister(max_measurements, "anc_meas")

    def get_circuit_registers(self) -> list:
        regs = [*self.qr_data]
        if self.use_ancilla_pool:
            if self.match_dynamic:
                regs.extend([self.qr_ctrl_par, self.qr_ctrl_rec])
                regs.extend(self.cr_pool)
                regs.extend([self.cr_rec, self.cr_readout])
            else:
                regs.extend(self.qr_ancilla)
                regs.append(self.cr_readout)
                regs.append(self.cr_ancilla)
        return regs


# ==========================================
# NOISE STRATEGIES
# ==========================================

class NoiseStrategy(ABC):
    @abstractmethod
    def apply_noise(self, qc, registers, epsilon):
        pass

    def generate_bindings(self, circuit, num_randomizations, epsilon):
        return np.empty((num_randomizations, 0))


class ParameterizedPauliTwirlingStrategy(NoiseStrategy):
    """Applies parameterized RZ/RX rotations to simulate Pauli twirling noise."""

    def __init__(self, k):
        self.k = k

    def apply_noise(self, qc, registers, epsilon):
        current_param_count = len(qc.parameters)
        for reg in registers:
            register_uid = current_param_count
            current_param_count += 2 * len(reg)
            for i, q in enumerate(reg):
                theta_x = Parameter(f"twirl_{register_uid}_q{i}_x")
                theta_z = Parameter(f"twirl_{register_uid}_q{i}_z")
                qc.append(RZGate(theta_z), [q])
                qc.append(RXGate(theta_x), [q])

    def generate_bindings(self, circuit, num_randomizations, epsilon):
        rng = np.random.default_rng()
        circuit_params = circuit.parameters
        num_params = len(circuit_params)
        if num_params == 0:
            return np.empty((num_randomizations, 0))

        bindings_matrix = np.zeros((num_randomizations, num_params))
        reg_map = {}

        for i, param in enumerate(circuit_params):
            name = param.name
            if name.startswith("twirl_"):
                parts = name.split('_')
                if len(parts) >= 4:
                    r_uid = int(parts[1])
                    q_idx = int(parts[2][1:])
                    p_type = parts[3]
                    reg_map.setdefault(r_uid, {}).setdefault(q_idx, {})[p_type] = i

        vals_x = [0.0, np.pi, 0.0, np.pi]
        vals_z = [0.0, 0.0, np.pi, np.pi]

        for r_uid, qubits_map in reg_map.items():
            k = len(qubits_map)
            num_paulis = 4 ** k
            is_error = rng.random(size=num_randomizations) < epsilon
            error_indices = rng.integers(0, num_paulis, size=num_randomizations)
            final_indices = np.where(is_error, error_indices, 0)

            current_val = final_indices
            for q_idx in sorted(qubits_map.keys()):
                p_choice = current_val % 4
                current_val = current_val // 4
                p_indices = qubits_map[q_idx]
                if 'x' in p_indices:
                    bindings_matrix[:, p_indices['x']] = [vals_x[c] for c in p_choice]
                if 'z' in p_indices:
                    bindings_matrix[:, p_indices['z']] = [vals_z[c] for c in p_choice]

        return bindings_matrix


# ==========================================
# BASE STRATEGY
# ==========================================

class QPAStrategy(ABC):
    """
    Base class for QPA circuit building strategies.

    Parameters
    ----------
    k : int
        Number of qubits per data register.
    n_trials : int
        Number of Schur-test rounds.
    n_registers : int
        Total data registers (must be odd, 2i+1). Defaults to 3.
    """

    def __init__(self, k, n_trials, n_registers=None):
        self.k = k
        self.n_trials = n_trials
        self.n_registers = n_registers if n_registers else 3
        assert self.n_registers % 2 == 1, "n_registers must be odd (2i+1)"
        self.noise_strategy: Optional[NoiseStrategy] = None

    def set_noise_strategy(self, strategy: NoiseStrategy):
        self.noise_strategy = strategy

    @abstractmethod
    def build_circuit(self, epsilon=0.0):
        """Build and return the QPA circuit(s)."""
        pass

    def apply_global_noise(self, qc, registers, epsilon):
        """Apply depolarizing noise to all data registers (simulator only)."""
        if epsilon <= 0:
            return
        noise = depolarizing_error(epsilon, self.k)
        for reg in registers:
            qc.append(noise, reg)


# ==========================================
# DYNAMIC CIRCUIT STRATEGY
# ==========================================

class HybridStrategy(QPAStrategy):
    """
    Dynamic QPA circuit with classical feedback.

    Produces a single QuantumCircuit that uses mid-circuit measurement and
    ``if_test`` to branch at runtime based on Schur test outcomes.
    This is the most compact representation but requires dynamic-circuit
    support from the backend.

    Usage
    -----
    >>> strategy = HybridStrategy(k=2, n_trials=3, n_registers=5)
    >>> qc = strategy.build_circuit(epsilon=0.1)
    """

    def __init__(self, k, n_trials, n_registers=5):
        super().__init__(k, n_trials, n_registers)

    def build_circuit(self, epsilon=0.0):
        qr_data = [QuantumRegister(self.k, f"R{i+1}") for i in range(self.n_registers)]
        max_parallel_tests = (self.n_registers - 1) // 2
        qr_ctrl_par = QuantumRegister(max_parallel_tests, "ctrl_par")
        qr_ctrl_rec = QuantumRegister(1, "ctrl_rec")
        cr_pool = [ClassicalRegister(max_parallel_tests, f"res_t{t}") for t in range(self.n_trials)]
        cr_rec = ClassicalRegister(self.n_trials, "res_rec")
        cr_final = ClassicalRegister(self.k, "readout")

        regs = [*qr_data, qr_ctrl_par, qr_ctrl_rec, *cr_pool, cr_rec, cr_final]
        qc = QuantumCircuit(*regs)

        self.apply_global_noise(qc, qr_data, epsilon)

        initial_reserve = qr_data[-1]
        initial_pairs = [(qr_data[i], qr_data[i + 1]) for i in range(0, self.n_registers - 1, 2)]

        self._build_recursive_layer(
            qc, initial_pairs, initial_reserve,
            qr_ctrl_par, qr_ctrl_rec,
            cr_pool, cr_rec, 0,
        )

        qc.measure(initial_reserve, cr_final)
        return qc

    def _build_recursive_layer(self, qc, current_pairs, reserve_reg,
                               qr_par, qr_rec, cr_pool, cr_rec, current_trial):
        if current_trial >= self.n_trials:
            return

        num_pairs = len(current_pairs)
        current_cr = cr_pool[current_trial]

        # 1. Parallel Schur tests
        for i in range(num_pairs):
            apply_schur_test(
                qc, qr_par[i], current_pairs[i][0], current_pairs[i][1],
                self.k, include_measure=True, res_bit=current_cr[i],
            )

        # 2. Branch on every outcome combination
        outcomes = list(product([0, 1], repeat=num_pairs))

        for outcome in outcomes:
            conditions = [(current_cr[i], val) for i, val in enumerate(outcome)]

            def apply_conditions(cond_list, block_func):
                if not cond_list:
                    block_func()
                    return
                head, *tail = cond_list
                with qc.if_test(head):
                    apply_conditions(tail, block_func)

            def logic_block(outcome=outcome):
                if all(v == 0 for v in outcome):
                    # All tests passed: rotate all registers and continue
                    all_regs = []
                    for p in current_pairs:
                        all_regs.extend(p)
                    all_regs.append(reserve_reg)
                    apply_cyclic_rotation(qc, all_regs, self.k)
                    self._build_recursive_layer(
                        qc, current_pairs, reserve_reg,
                        qr_par, qr_rec, cr_pool, cr_rec, current_trial + 1,
                    )
                else:
                    # Some tests failed: keep only surviving pairs
                    surviving_pairs = [
                        current_pairs[i] for i, val in enumerate(outcome) if val == 0
                    ]
                    if surviving_pairs:
                        active_regs = []
                        for p in surviving_pairs:
                            active_regs.extend(p)
                        active_regs.append(reserve_reg)
                        apply_cyclic_rotation(qc, active_regs, self.k)
                        self._build_recursive_layer(
                            qc, surviving_pairs, reserve_reg,
                            qr_par, qr_rec, cr_pool, cr_rec, current_trial + 1,
                        )

            apply_conditions(conditions, logic_block)


# Backward-compatible alias
HybridNRegStrategy = HybridStrategy


class HybridStrategy_n5_ntrials_1(QPAStrategy):
    """
    Dynamic QPA circuit for N=5 with a single trial round, using flat
    (non-nested) ``if_test`` on the 2-bit classical register value.

    With N=5 and T=1, there are 2 parallel Schur test pairs producing a
    2-bit result register.  Instead of nesting ``if_test`` per bit, this
    strategy branches on the combined register value (0–3):

        value 0 (00): both pass  → full cyclic rotation of all 5 registers
        value 1 (01): pair 0 fail, pair 1 pass → rotate R3, R4, R5
        value 2 (10): pair 0 pass, pair 1 fail → rotate R1, R2, R5
        value 3 (11): both fail  → no operation (omitted)

    This avoids the unnecessary nested-if structure that HybridStrategy
    produces for this case, resulting in a simpler circuit for transpilation.

    Parameters
    ----------
    k : int
        Number of qubits per data register (qudit dimension d = 2^k).
    """

    def __init__(self, k):
        super().__init__(k, n_trials=1, n_registers=5)

    def build_circuit(self, epsilon=0.0):
        qr_data = [QuantumRegister(self.k, f"R{i+1}") for i in range(5)]
        qr_ctrl_par = QuantumRegister(2, "ctrl_par")
        qr_ctrl_rec = QuantumRegister(1, "ctrl_rec")
        cr_par = ClassicalRegister(2, "res_t0")
        cr_rec = ClassicalRegister(1, "res_rec")
        cr_final = ClassicalRegister(self.k, "readout")

        qc = QuantumCircuit(
            *qr_data, qr_ctrl_par, qr_ctrl_rec,
            cr_par, cr_rec, cr_final,
        )

        self.apply_global_noise(qc, qr_data, epsilon)

        reserve = qr_data[4]  # R5
        pair0 = (qr_data[0], qr_data[1])  # (R1, R2)
        pair1 = (qr_data[2], qr_data[3])  # (R3, R4)

        # Parallel Schur tests
        apply_schur_test(
            qc, qr_ctrl_par[0], pair0[0], pair0[1],
            self.k, include_measure=True, res_bit=cr_par[0],
        )
        apply_schur_test(
            qc, qr_ctrl_par[1], pair1[0], pair1[1],
            self.k, include_measure=True, res_bit=cr_par[1],
        )

        # Value 0 (00): both pass → rotate all [R1, R2, R3, R4, R5]
        with qc.if_test((cr_par, 0)):
            apply_cyclic_rotation(qc, [*pair0, *pair1, reserve], self.k)

        # Value 2 (10): pair 0 pass, pair 1 fail → rotate [R1, R2, R5]
        with qc.if_test((cr_par, 2)):
            apply_cyclic_rotation(qc, [pair0[0], pair0[1], reserve], self.k)

        # Value 1 (01): pair 0 fail, pair 1 pass → rotate [R3, R4, R5]
        with qc.if_test((cr_par, 1)):
            apply_cyclic_rotation(qc, [pair1[0], pair1[1], reserve], self.k)

        # Value 3 (11): both fail → no-op (omitted entirely)

        qc.measure(reserve, cr_final)
        return qc


# ==========================================
# STATIC CIRCUIT STRATEGY
# ==========================================

class UnrolledStrategy(QPAStrategy):
    """
    Static QPA circuits -- one circuit per execution path.

    Enumerates all possible Schur-test outcome combinations and produces
    a separate static circuit for each path.  No dynamic control flow is
    needed, so these circuits run on any backend.

    The ``build_longest_path`` method produces only the single path where
    every Schur test succeeds (outcome 0), which is the deepest static
    circuit and matches the "else" branch of the dynamic HybridStrategy.

    Usage
    -----
    >>> strategy = UnrolledStrategy(k=2, n_trials=3, n_registers=5)
    >>> # All paths
    >>> all_paths = strategy.build_circuit(epsilon=0.0)
    >>> # Longest path only (all Schur tests pass)
    >>> longest = strategy.build_longest_path(epsilon=0.0)
    """

    def __init__(self, k, n_trials, n_registers=None, no_reset=False,
                 match_dynamic_registers=False):
        super().__init__(k, n_trials, n_registers)
        self.no_reset = no_reset
        self.match_dynamic_registers = match_dynamic_registers
        self.circuits_data: List[Dict[str, Any]] = []

    # ------ public API (also kept as `build` for backward compat) ------

    def build_circuit(self, epsilon=0.0):
        """Build static circuits for all execution paths."""
        return self.build(epsilon)

    def build(self, epsilon=0.0):
        """Build static circuits for all execution paths.

        Returns list of dicts with keys: 'circuit', 'conditions', 'metadata'.
        """
        self.circuits_data = []
        regs = QPARegisters(
            self.n_registers, self.k, self.n_trials,
            use_ancilla_pool=True, no_reset=self.no_reset,
            match_dynamic=self.match_dynamic_registers,
        )
        qc_template = QuantumCircuit(*regs.get_circuit_registers())

        if self.noise_strategy:
            self.noise_strategy.apply_noise(qc_template, regs.qr_data, epsilon)

        initial_pairs = [[i, i + 1] for i in range(0, self.n_registers - 1, 2)]
        reserve = self.n_registers - 1

        self._recurse(
            qc_template, initial_pairs, reserve,
            trial=0, conditions={}, total_meas_index=0, regs=regs,
        )
        return self.circuits_data

    def build_longest_path(self, epsilon=0.0):
        """
        Build only the longest execution path (all Schur tests pass every trial).

        This path corresponds to the heaviest branch of the dynamic
        HybridStrategy circuit and is useful for depth/gate-count analysis
        or for running the worst-case static circuit on hardware.

        Returns a single dict with keys: 'circuit', 'conditions', 'metadata'.
        """
        self.circuits_data = []
        regs = QPARegisters(
            self.n_registers, self.k, self.n_trials,
            use_ancilla_pool=True, no_reset=self.no_reset,
            match_dynamic=self.match_dynamic_registers,
        )
        qc_template = QuantumCircuit(*regs.get_circuit_registers())

        if self.noise_strategy:
            self.noise_strategy.apply_noise(qc_template, regs.qr_data, epsilon)

        initial_pairs = [[i, i + 1] for i in range(0, self.n_registers - 1, 2)]
        reserve = self.n_registers - 1

        self._recurse_longest(
            qc_template, initial_pairs, reserve,
            trial=0, conditions={}, total_meas_index=0, regs=regs,
        )
        return self.circuits_data[0]

    # ------ internal recursion ------

    def _recurse(self, current_qc, current_pairs, reserve_idx,
                 trial, conditions, total_meas_index, regs):
        if trial >= self.n_trials:
            self._finalize_path(current_qc, conditions, reserve_idx, regs)
            return

        num_pairs = len(current_pairs)
        outcomes = list(product([0, 1], repeat=num_pairs))

        for outcome in outcomes:
            branch_qc = current_qc.copy()
            branch_conditions = conditions.copy()

            self._apply_schur_tests(
                branch_qc, current_pairs, outcome,
                trial, total_meas_index, branch_conditions, regs,
            )

            surviving_pairs = [current_pairs[i] for i, res in enumerate(outcome) if res == 0]

            if surviving_pairs:
                survivor_flat = [idx for p in surviving_pairs for idx in p]
                active_indices = survivor_flat + [reserve_idx]
                apply_cyclic_rotation_indices(branch_qc, active_indices, self.k)

                new_pairs = [
                    [active_indices[i], active_indices[i + 1]]
                    for i in range(0, len(active_indices) - 1, 2)
                ]
                new_reserve = active_indices[-1]
                self._recurse(
                    branch_qc, new_pairs, new_reserve, trial + 1,
                    branch_conditions, total_meas_index + num_pairs, regs,
                )
            else:
                self._finalize_path(branch_qc, branch_conditions, reserve_idx, regs)

    def _recurse_longest(self, current_qc, current_pairs, reserve_idx,
                         trial, conditions, total_meas_index, regs):
        """Like _recurse but only follows the all-success (all 0) branch."""
        if trial >= self.n_trials:
            self._finalize_path(current_qc, conditions, reserve_idx, regs, tag="longest_path")
            return

        num_pairs = len(current_pairs)
        outcome = tuple(0 for _ in range(num_pairs))  # all pass

        branch_qc = current_qc.copy()
        branch_conditions = conditions.copy()

        self._apply_schur_tests(
            branch_qc, current_pairs, outcome,
            trial, total_meas_index, branch_conditions, regs,
        )

        # All pairs survive
        all_flat = [idx for p in current_pairs for idx in p]
        active_indices = all_flat + [reserve_idx]
        apply_cyclic_rotation_indices(branch_qc, active_indices, self.k)

        new_pairs = [
            [active_indices[i], active_indices[i + 1]]
            for i in range(0, len(active_indices) - 1, 2)
        ]
        new_reserve = active_indices[-1]

        self._recurse_longest(
            branch_qc, new_pairs, new_reserve, trial + 1,
            branch_conditions, total_meas_index + num_pairs, regs,
        )

    def _apply_schur_tests(self, qc, pairs, outcome, trial, total_meas_index, conditions, regs):
        """Apply Schur tests + measurements for one trial on the given pairs."""
        for i, pair in enumerate(pairs):
            if self.no_reset:
                num_concurrent = self.n_registers // 2
                anc_idx = trial * num_concurrent + i
                anc_qubit = regs.qr_ancilla[anc_idx]
            else:
                anc_qubit = regs.qr_ancilla[i]

            rA_idx, rB_idx = pair
            rA = qc.qubits[rA_idx * self.k: (rA_idx + 1) * self.k]
            rB = qc.qubits[rB_idx * self.k: (rB_idx + 1) * self.k]

            apply_schur_test(qc, anc_qubit, rA, rB, self.k)

            cl_idx = total_meas_index + i
            qc.measure(anc_qubit, regs.cr_ancilla[cl_idx])

            if not self.no_reset:
                qc.reset(anc_qubit)

            conditions[self.k + cl_idx] = outcome[i]

    def _finalize_path(self, qc, conditions, reserve_idx, regs, tag=None):
        final_qc = qc.copy()
        r_qubits = final_qc.qubits[reserve_idx * self.k: (reserve_idx + 1) * self.k]
        for i in range(self.k):
            final_qc.measure(r_qubits[i], regs.cr_readout[i])

        path_name = tag or f"path_{len(self.circuits_data)}"
        self.circuits_data.append({
            'circuit': final_qc,
            'conditions': conditions,
            'metadata': {'type': 'unrolled', 'path_name': path_name},
        })


# Keep old name available for backward compatibility
CircuitGenerationStrategy = QPAStrategy
