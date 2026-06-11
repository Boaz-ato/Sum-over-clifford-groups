from __future__ import annotations

import math

import qecc as q


def generate_clifford_groups(n: int = 1, with_phases: bool = False) -> None:
    """Generate the full Clifford group for a given number of qubits and write it to disk.

    The Clifford group consists of all unitary operations that map Pauli operators
    to Pauli operators under conjugation. This function iterates over all group
    elements using the ``qecc`` library, optionally multiplies each element by
    8 equally-spaced phase factors (multiples of π/4), and writes both the real
    and imaginary parts of each matrix to a tab-separated text file.

    Supported sizes:
        - n=1: 24 Clifford gates, matrix size 2×2
        - n=2: 11520 Clifford gates, matrix size 4×4

    Args:
        n: Number of qubits (must be 1 or 2).
        with_phases: If True, each Clifford element is multiplied by all 8
            eighth-roots of unity (e^{iπk/4} for k=1,...,8), expanding the
            output by a factor of 8. Defaults to False.

    Raises:
        ValueError: If n is not 1 or 2.

    Side effects:
        Writes to ``datasets/full_cliffords/Clifford_group_<n>.txt`` (or
        ``Clifford_group_<n>_with_phases.txt`` if ``with_phases=True``).
        Each matrix is written as two consecutive blocks of rows: first the
        real part, then the imaginary part, each in row-major order.
    """
    n = int(n)
    if n not in (1, 2):
        raise ValueError("n must be 1 or 2")

    iterator = q.clifford_group(n, True)

    if with_phases:
        phases = [
            complex(math.cos(k * math.pi / 4), math.sin(k * math.pi / 4))
            for k in range(1, 9)
        ]
    else:
        phases = [1]

    num_clifford_gates = 24 if n == 1 else 11520

    output_path = (
        f'datasets/full_cliffords/Clifford_group_{n}_with_phases.txt'
        if with_phases
        else f'datasets/full_cliffords/Clifford_group_{n}.txt'
    )

    nmatrix = 2 ** n
    with open(output_path, 'w') as f:
        for _ in range(num_clifford_gates):
            base_matrix = next(iterator).as_unitary()
            for phase in phases:
                matrix = base_matrix * phase
                for i in range(nmatrix):
                    for j in range(nmatrix):
                        sep = '\n' if j == nmatrix - 1 else '\t'
                        f.write(f'{matrix[i][j].real:f}{sep}')
                for i in range(nmatrix):
                    for j in range(nmatrix):
                        sep = '\n' if j == nmatrix - 1 else '\t'
                        f.write(f'{matrix[i][j].imag:f}{sep}')


if __name__ == "__main__":
    generate_clifford_groups(2, with_phases=False)
