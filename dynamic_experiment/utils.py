from collections import defaultdict
from qiskit import QuantumCircuit
from qiskit.qasm3 import dumps, loads as parse
from qiskit.transpiler import Layout
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

def strip_if_keep_else(qasm: str) -> str:
    """
    Remove all `if (...) { ... } else { ... }` blocks from a QASM 3 string,
    keeping only the contents of the `else` block. If an `if` has no `else`,
    the entire `if` statement is removed.

    This version avoids adding blank lines after the `else` body.
    """
    i, n = 0, len(qasm)
    out = []

    def peek(offset=0):
        j = i + offset
        return qasm[j] if 0 <= j < n else ''

    def is_word_char(ch):
        return ch.isalnum() or ch == '_'

    def skip_ws(j):
        while j < n and qasm[j].isspace():
            j += 1
        return j

    def parse_paren_block(j, open_ch='(', close_ch=')'):
        depth, k = 1, j + 1
        inner_start = k
        while k < n and depth:
            if qasm[k] == open_ch:
                depth += 1
            elif qasm[k] == close_ch:
                depth -= 1
            k += 1
        if depth != 0:
            raise ValueError("Unbalanced parentheses/braces in input.")
        return k, inner_start, k - 1

    def parse_brace_block(j):
        return parse_paren_block(j, '{', '}')

    while i < n:
        if (peek() == 'i' and peek(1) == 'f'
                and (i == 0 or not is_word_char(qasm[i - 1]))
                and not is_word_char(peek(2))):
            j = i + 2
            j = skip_ws(j)
            if j >= n or qasm[j] != '(':
                out.append(qasm[i])
                i += 1
                continue
            j, _, _ = parse_paren_block(j, '(', ')')
            j = skip_ws(j)
            if j >= n or qasm[j] != '{':
                out.append(qasm[i])
                i += 1
                continue
            j_after_then, _, _ = parse_brace_block(j)
            j = skip_ws(j_after_then)

            if qasm[j:j+4] == 'else':
                k = skip_ws(j + 4)
                if k < n and qasm[k] == '{':
                    k_after_else, else_start, else_end = parse_brace_block(k)
                    else_body = qasm[else_start:else_end].strip("\n\r ")
                    # Dedent and reformat cleanly
                    lines = [line.strip() for line in else_body.splitlines() if line.strip()]
                    out.append("\n".join(lines))
                    i = k_after_else
                    continue
            i = j  # drop the if if no else
            continue
        out.append(qasm[i])
        i += 1

    # Post-process: collapse multiple blank lines into one
    result = ''.join(out)
    clean_lines = [line.rstrip() for line in result.splitlines()]
    final_lines = []
    prev_blank = False
    for line in clean_lines:
        if line.strip() == "":
            if not prev_blank:
                final_lines.append("")
            prev_blank = True
        else:
            final_lines.append(line)
            prev_blank = False
    return "\n".join(final_lines)

def strip_if_keep_then(qasm: str) -> str:
    """
    Remove all `if (...) { ... } else { ... }` blocks from a QASM 3 string,
    keeping only the contents of the `if` (then) block.
    If an `if` has no `else`, the `if` body is kept.
    """
    i, n = 0, len(qasm)
    out = []

    def peek(offset=0):
        j = i + offset
        return qasm[j] if 0 <= j < n else ''

    def is_word_char(ch):
        return ch.isalnum() or ch == '_'

    def skip_ws(j):
        while j < n and qasm[j].isspace():
            j += 1
        return j

    def parse_paren_block(j, open_ch='(', close_ch=')'):
        depth, k = 1, j + 1
        inner_start = k
        while k < n and depth:
            if qasm[k] == open_ch:
                depth += 1
            elif qasm[k] == close_ch:
                depth -= 1
            k += 1
        if depth != 0:
            raise ValueError("Unbalanced parentheses/braces in input.")
        return k, inner_start, k - 1

    def parse_brace_block(j):
        return parse_paren_block(j, '{', '}')

    while i < n:
        if (peek() == 'i' and peek(1) == 'f'
                and (i == 0 or not is_word_char(qasm[i - 1]))
                and not is_word_char(peek(2))):

            j = i + 2
            j = skip_ws(j)

            if j >= n or qasm[j] != '(':
                out.append(qasm[i])
                i += 1
                continue

            # Skip condition (...)
            j, _, _ = parse_paren_block(j, '(', ')')
            j = skip_ws(j)

            if j >= n or qasm[j] != '{':
                out.append(qasm[i])
                i += 1
                continue

            # Extract THEN block
            j_after_then, then_start, then_end = parse_brace_block(j)
            then_body = qasm[then_start:then_end].strip()

            # Append cleaned THEN body
            lines = [line.strip() for line in then_body.splitlines() if line.strip()]
            out.append("\n".join(lines))

            # Skip optional ELSE block entirely
            j = skip_ws(j_after_then)
            if qasm[j:j+4] == 'else':
                k = skip_ws(j + 4)
                if k < n and qasm[k] == '{':
                    j_after_else, _, _ = parse_brace_block(k)
                    i = j_after_else
                    continue

            i = j_after_then
            continue

        out.append(qasm[i])
        i += 1

    return ''.join(out)

def get_dynamic_part(qc):
    """ Copies only the dynamics part of the circuit. """
    out = QuantumCircuit(*qc.qregs, *qc.cregs, name=f"{qc.name}_only_ifelse")
    out.global_phase = qc.global_phase
    if qc.metadata is not None:
        out.metadata = dict(qc.metadata)

    for instr in qc.data:               # instr is a CircuitInstruction
        op = instr.operation
        if getattr(op, "name", "") != "if_else":
            continue                    # drop the control-flow instruction
        out.append(op, instr.qubits, instr.clbits)
    return out

def convert_dynamic_to_else(qc):
    """
    Convert a quantum circuit with dynamics to a circuit with `if` statements
    that can be executed on a classical simulator.
    """
    qasm_string = dumps(qc)
    qasm_string = strip_if_keep_else(qasm_string)
    qc_result = parse(qasm_string)

    return qc_result

def convert_dynamic_to_if(qc):
    qasm_string = dumps(qc)
    qasm_string = strip_if_keep_then(qasm_string)
    qc_result = parse(qasm_string)
    
    return qc_result

def get_layout_from_static(qc_dynamic, backend, branch='else', optimization_level=1, seed_transpiler=42):
    """
    Transpile the static (if or else) branch of a dynamic circuit to extract
    a good initial virtual layout for routing the full dynamic circuit.

    Args:
        qc_dynamic: Dynamic QuantumCircuit with if/else control flow.
        backend: Target backend.
        branch: 'else' or 'if' — which branch to use as the static proxy.
        optimization_level: Optimization level for the preset pass manager.
        seed_transpiler: Seed for reproducibility.

    Returns:
        Layout object suitable for use as initial_layout.
    """
    if branch == 'else':
        qc_static = convert_dynamic_to_else(qc_dynamic)
    elif branch == 'if':
        qc_static = convert_dynamic_to_if(qc_dynamic)
    else:
        raise ValueError(f"branch must be 'if' or 'else', got {branch!r}")

    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
    )
    qc_static_tr = pm.run(qc_static)
    layout = qc_static_tr.layout.initial_virtual_layout(filter_ancillas=True)

    # QASM round-trip creates new qubit objects; remap virtual qubits from the
    # static circuit back to qc_dynamic's qubits by matching global qubit index.
    static_to_dynamic = {s: qc_dynamic.qubits[i] for i, s in enumerate(qc_static.qubits)}
    remapped = Layout({
        phys: static_to_dynamic[virt]
        for phys, virt in layout.get_physical_bits().items()
        if virt in static_to_dynamic
    })
    return remapped


def transpile_dynamic_with_static_layout(qc_dynamic, backend, branch='else', optimization_level=3, seed_transpiler=42):
    """
    Transpile a dynamic circuit using a layout derived from its static branch.

    Qiskit's routing algorithms perform better on static circuits. This function
    extracts a good initial layout from the if or else branch, then uses it as
    the starting point when transpiling the full dynamic circuit.

    Args:
        qc_dynamic: Dynamic QuantumCircuit with if/else control flow.
        backend: Target backend.
        branch: 'else' or 'if' — which branch to use as the static proxy.
        optimization_level: Optimization level for the preset pass manager.
        seed_transpiler: Seed for reproducibility.

    Returns:
        (qc_transpiled, initial_layout): The transpiled dynamic circuit and the
        layout that was used.
    """
    initial_layout = get_layout_from_static(
        qc_dynamic, backend, branch=branch,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
    )
    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
        initial_layout=initial_layout,
    )
    qc_transpiled = pm.run(qc_dynamic)
    return qc_transpiled, initial_layout


def transpile_dynamic_best_of_n(
    qc_dynamic,
    backend,
    branch: str = 'else',
    optimization_level: int = 1,
    n_trials: int = 20,
    seeds=None,
    depth_branch: str = 'else',
    filter_fn=None,
):
    """
    Transpile a dynamic circuit multiple times with different seeds and return
    the result with the lowest static depth.

    For each trial the full static-layout strategy is used:
    a layout is extracted from the static branch, then the dynamic circuit is
    re-transpiled with that layout.

    Args:
        qc_dynamic: Dynamic QuantumCircuit with if/else control flow.
        backend: Target backend.
        branch: Branch used for layout extraction ('else' or 'if').
        optimization_level: Optimization level for the preset pass managers.
        n_trials: Number of random seeds to try (ignored when seeds is given).
        seeds: Explicit list of integer seeds. If None, uses range(n_trials).
        depth_branch: Branch used to evaluate depth ('else' or 'if').
        filter_fn: Optional callable for QuantumCircuit.depth() filtering
                   (e.g. ``lambda x: x.operation.nulamm_qubits == 2``).

    Returns:
        dict with keys:
            'circuit'    – best transpiled QuantumCircuit
            'seed'       – seed that produced the best result
            'best_depth' – lowest static depth found
            'avg_depth'  – mean static depth across all trials
            'all_depths' – list of (seed, depth) for every trial
    """
    if seeds is None:
        seeds = list(range(n_trials))

    all_depths = []
    best_depth = None
    best_circuit = None
    best_seed = None

    for seed in seeds:
        qc_tr, _ = transpile_dynamic_with_static_layout(
            qc_dynamic,
            backend,
            branch=branch,
            optimization_level=optimization_level,
            seed_transpiler=seed,
        )
        depth = get_static_depth(qc_tr, branch=depth_branch, filter_fn=filter_fn)
        all_depths.append((seed, depth))

        if best_depth is None or depth < best_depth:
            best_depth = depth
            best_circuit = qc_tr
            best_seed = seed

    avg_depth = sum(d for _, d in all_depths) / len(all_depths)

    return {
        'circuit': best_circuit,
        'seed': best_seed,
        'best_depth': best_depth,
        'avg_depth': avg_depth,
        'all_depths': all_depths,
    }


def _flatten_branch(qc: QuantumCircuit, branch: str, reps: int) -> QuantumCircuit:
    """Apply branch flattening `reps` times to handle nested if_else blocks."""
    if branch == 'else':
        convert = convert_dynamic_to_else
    elif branch == 'if':
        convert = convert_dynamic_to_if
    else:
        raise ValueError(f"branch must be 'if' or 'else', got {branch!r}")
    qc_static = qc
    for _ in range(reps):
        qc_static = convert(qc_static)
    return qc_static


def get_static_gate_counts(qc: QuantumCircuit, branch: str = 'else', reps: int = 1) -> dict:
    """
    Return gate counts for a dynamic circuit by first converting it to a
    static circuit via the specified branch.

    Unlike QuantumCircuit.count_ops(), this descends into if_else blocks so
    the result reflects the actual gates that would be executed on the chosen
    execution path, rather than listing 'if_else' as a single opaque operation.

    For circuits with nested if_else blocks, increase `reps` to flatten
    multiple levels of branching in one call.

    Args:
        qc: QuantumCircuit (may contain if_else control flow).
        branch: 'else' or 'if' — which branch to flatten into.
        reps: Number of times to apply the branch conversion (default 1).
              Use reps > 1 for circuits with nested if_else blocks.

    Returns:
        dict mapping gate name → count.
    """
    return dict(_flatten_branch(qc, branch, reps).count_ops())


def get_static_depth(qc: QuantumCircuit, branch: str = 'else', filter_fn=None, reps: int = 1) -> int:
    """
    Return a more accurate depth for a dynamic circuit by first converting it
    to a static circuit via the specified branch, then measuring depth.

    For QPA-style circuits the 'else' branch is the heavier one (it contains
    the ancilla SWAPs), so it gives a truer upper-bound on execution depth.

    For circuits with nested if_else blocks, increase `reps` to flatten
    multiple levels of branching in one call.

    Args:
        qc: QuantumCircuit (may contain if_else control flow).
        branch: 'else' or 'if' — which branch to flatten into.
        filter_fn: Optional callable passed to QuantumCircuit.depth() to
                   filter which instructions count (e.g. 2Q-gate filter).
        reps: Number of times to apply the branch conversion (default 1).
              Use reps > 1 for circuits with nested if_else blocks.

    Returns:
        Depth of the static circuit (int).
    """
    qc_static = _flatten_branch(qc, branch, reps)
    return qc_static.depth(filter_fn) if filter_fn is not None else qc_static.depth()


def get_layout_from_unrolled(qc_dynamic, backend, k, n_trials, n_registers,
                              optimization_level=1, seed_transpiler=42):
    """
    Transpile an UnrolledStrategy longest-path circuit to extract a good
    initial virtual layout for routing the full dynamic circuit.

    Instead of flattening the dynamic circuit's if/else branches, this builds
    a separate static circuit via UnrolledStrategy.build_longest_path with
    match_dynamic_registers=True so the qubit registers match exactly.

    Args:
        qc_dynamic: Dynamic QuantumCircuit (HybridStrategy) to eventually transpile.
        backend: Target backend.
        k: Qubits per data register.
        n_trials: Number of Schur-test rounds.
        n_registers: Number of data registers (must be odd).
        optimization_level: Optimization level for the preset pass manager.
        seed_transpiler: Seed for reproducibility.

    Returns:
        Layout object suitable for use as initial_layout on qc_dynamic.
    """
    from qpa_engine import UnrolledStrategy

    strat = UnrolledStrategy(
        k=k, n_trials=n_trials, n_registers=n_registers,
        match_dynamic_registers=True,
    )
    qc_static = strat.build_longest_path()['circuit']

    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
    )
    qc_static_tr = pm.run(qc_static)
    layout = qc_static_tr.layout.initial_virtual_layout(filter_ancillas=True)

    # Remap virtual qubits from the unrolled circuit back to the dynamic circuit's
    # qubits by matching global qubit index (registers match by construction).
    static_to_dynamic = {s: qc_dynamic.qubits[i] for i, s in enumerate(qc_static.qubits)}
    remapped = Layout({
        phys: static_to_dynamic[virt]
        for phys, virt in layout.get_physical_bits().items()
        if virt in static_to_dynamic
    })
    return remapped


def transpile_dynamic_with_unrolled_layout(qc_dynamic, backend, k, n_trials, n_registers,
                                            optimization_level=3, seed_transpiler=42):
    """
    Transpile a dynamic circuit using a layout derived from the UnrolledStrategy
    longest-path static circuit.

    Args:
        qc_dynamic: Dynamic QuantumCircuit (HybridStrategy).
        backend: Target backend.
        k: Qubits per data register.
        n_trials: Number of Schur-test rounds.
        n_registers: Number of data registers (must be odd).
        optimization_level: Optimization level for the preset pass manager.
        seed_transpiler: Seed for reproducibility.

    Returns:
        (qc_transpiled, initial_layout): The transpiled dynamic circuit and the
        layout that was used.
    """
    initial_layout = get_layout_from_unrolled(
        qc_dynamic, backend,
        k=k, n_trials=n_trials, n_registers=n_registers,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
    )
    pm = generate_preset_pass_manager(
        backend=backend,
        optimization_level=optimization_level,
        seed_transpiler=seed_transpiler,
        initial_layout=initial_layout,
    )
    qc_transpiled = pm.run(qc_dynamic)
    return qc_transpiled, initial_layout


def transpile_dynamic_best_of_n_unrolled(
    qc_dynamic,
    backend,
    k: int,
    n_trials_qpa: int,
    n_registers: int,
    optimization_level: int = 1,
    n_trials: int = 20,
    seeds=None,
    depth_branch: str = 'if',
    filter_fn=None,
    reps: int = 1,
):
    """
    Transpile a dynamic circuit multiple times with different seeds using the
    unrolled-layout strategy and return the result with the lowest static depth.

    Args:
        qc_dynamic: Dynamic QuantumCircuit (HybridStrategy).
        backend: Target backend.
        k: Qubits per data register.
        n_trials_qpa: Number of QPA Schur-test rounds (passed to UnrolledStrategy).
        n_registers: Number of data registers.
        optimization_level: Optimization level for the preset pass managers.
        n_trials: Number of random seeds to try (ignored when seeds is given).
        seeds: Explicit list of integer seeds. If None, uses range(n_trials).
        depth_branch: Branch used to evaluate depth ('else' or 'if').
        filter_fn: Optional callable for depth filtering.
        reps: Number of times to apply branch conversion for depth evaluation.

    Returns:
        dict with keys:
            'circuit'    - best transpiled QuantumCircuit
            'seed'       - seed that produced the best result
            'best_depth' - lowest static depth found
            'avg_depth'  - mean static depth across all trials
            'all_depths' - list of (seed, depth) for every trial
    """
    if seeds is None:
        seeds = list(range(n_trials))

    all_depths = []
    best_depth = None
    best_circuit = None
    best_seed = None

    for seed in seeds:
        qc_tr, _ = transpile_dynamic_with_unrolled_layout(
            qc_dynamic, backend,
            k=k, n_trials=n_trials_qpa, n_registers=n_registers,
            optimization_level=optimization_level,
            seed_transpiler=seed,
        )
        depth = get_static_depth(qc_tr, branch=depth_branch, filter_fn=filter_fn, reps=reps)
        all_depths.append((seed, depth))

        if best_depth is None or depth < best_depth:
            best_depth = depth
            best_circuit = qc_tr
            best_seed = seed

    avg_depth = sum(d for _, d in all_depths) / len(all_depths)

    return {
        'circuit': best_circuit,
        'seed': best_seed,
        'best_depth': best_depth,
        'avg_depth': avg_depth,
        'all_depths': all_depths,
    }


def analyze_circuit(qc: QuantumCircuit):

    # Get the depth
    depth = qc.depth()
    depth_2q = qc.depth(lambda x: x.operation.num_qubits == 2)

    # Get the number of gates
    gate_counts = defaultdict(int)
    for inst in qc.data:
        gate_counts[len(inst.qubits)] += 1

    # Get types of gates
    gates = qc.count_ops()
    
    # print the results
    print(f"Depth of the circuit: {depth}")
    print(f"Depth of 2-qubit gates: {depth_2q}")
    print("Gate counts:")
    for num_qubits, count in gate_counts.items():
        print(f"  {num_qubits}-qubit gates: {count}")
    print("Types of gates:")
    for gate, count in gates.items():
        print(f"  {gate}: {count}")
    #print(f"Estimated duration of the circuit: {duration} microseconds")

    results = {
        "depth": depth,
        "depth_2q": depth_2q,
        "gate_counts": dict(gate_counts),
        "gates": dict(gates),
    }
    return (qc, results)

