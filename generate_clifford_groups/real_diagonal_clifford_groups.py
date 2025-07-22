import numpy as np


def one_qubit_diagonal_group():
    """Generate 1-qubit real diagonal Clifford group """
    P_diag = np.array([1, complex(0, 1)], dtype=complex)

    diagonal_group = []
    for i in range(1, 5):
        P_power_diag = P_diag ** i
        if np.allclose(P_power_diag.imag, 0):
            diagonal_group.append(P_power_diag.real)

    return np.array(diagonal_group)


def two_qubit_diagonal_group():
    """Generate 2-qubit real diagonal Clifford group"""
    P_diag = np.array([1, complex(0, 1)], dtype=complex)
    I_diag = np.array([1, 1], dtype=complex)
    CZ_diag = np.array([1, 1, 1, -1], dtype=complex)

    P1_diag = np.kron(P_diag, I_diag)
    P2_diag = np.kron(I_diag, P_diag)
    P1P2_diag = np.kron(P_diag, P_diag)

    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 3):
                    result_diag = (P1_diag ** i) * (P2_diag ** j) * (P1P2_diag ** k) * (CZ_diag ** l)
                    if np.allclose(result_diag.imag, 0):
                        diagonal_group.append(result_diag.real)

    return np.array(diagonal_group)


def create_cz_diagonal(qubit1, qubit2, n_qubits):
    """Create diagonal elements of CZ gate between qubit1 and qubit2 for n_qubits system."""
    dim = 2 ** n_qubits
    diagonal = np.ones(dim, dtype=complex)

    for i in range(dim):
        # Check if both qubits are in state |1> (using 1-based indexing)
        if ((i >> (n_qubits - qubit1)) & 1) and ((i >> (n_qubits - qubit2)) & 1):
            diagonal[i] = -1

    return diagonal


def create_pauli_z_diagonal(qubit, n_qubits):
    """Create diagonal elements of Pauli-Z gate on specific qubit for n_qubits system."""
    dim = 2 ** n_qubits
    P_single = np.array([1, complex(0, 1)], dtype=complex)
    I_single = np.array([1, 1], dtype=complex)

    # Build tensor product: I ⊗ ... ⊗ P ⊗ ... ⊗ I
    result = np.array([1], dtype=complex)
    for q in range(1, n_qubits + 1):
        if q == qubit:
            result = np.kron(result, P_single)
        else:
            result = np.kron(result, I_single)

    return result


def three_qubit_diagonal_group():
    """Generate 3-qubit real diagonal Clifford group"""
    n_qubits = 3

    # Create P gates for each qubit
    P_gates = []
    for qubit in range(1, n_qubits + 1):
        P_gates.append(create_pauli_z_diagonal(qubit, n_qubits))

    # Create CZ gates
    CZ12_diag = create_cz_diagonal(1, 2, n_qubits)
    CZ23_diag = create_cz_diagonal(2, 3, n_qubits)
    CZ13_diag = create_cz_diagonal(1, 3, n_qubits)

    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 3):
                    for m in range(1, 3):
                        for n in range(1, 3):
                            result_diag = (P_gates[0] ** i) * (P_gates[1] ** j) * (P_gates[2] ** k) * \
                                          (CZ12_diag ** l) * (CZ23_diag ** n) * (CZ13_diag ** m)
                            if np.allclose(result_diag.imag, 0):
                                diagonal_group.append(result_diag.real)

    return np.array(diagonal_group)


def four_qubit_diagonal_group():
    """Generate 4-qubit real diagonal Clifford group"""
    n_qubits = 4

    # Create P gates for each qubit
    P_gates = []
    for qubit in range(1, n_qubits + 1):
        P_gates.append(create_pauli_z_diagonal(qubit, n_qubits))

    # Create all CZ gates (6 pairs for 4 qubits)
    cz_gates = []
    pairs = [(1, 2), (2, 3), (3, 4), (1, 3), (1, 4), (2, 4)]
    for q1, q2 in pairs:
        cz_gates.append(create_cz_diagonal(q1, q2, n_qubits))

    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 3):
                        for n in range(1, 3):
                            for o in range(1, 3):
                                for p in range(1, 3):
                                    for q in range(1, 3):
                                        for r in range(1, 3):
                                            result_diag = (P_gates[0] ** i) * (P_gates[1] ** j) * \
                                                          (P_gates[2] ** k) * (P_gates[3] ** l) * \
                                                          (cz_gates[0] ** m) * (cz_gates[1] ** n) * \
                                                          (cz_gates[2] ** o) * (cz_gates[3] ** p) * \
                                                          (cz_gates[4] ** q) * (cz_gates[5] ** r)
                                            if np.allclose(result_diag.imag, 0):
                                                diagonal_group.append(result_diag.real)

    return np.array(diagonal_group)


def five_qubit_diagonal_group():
    """Generate 5-qubit real diagonal Clifford group."""
    n_qubits = 5

    # Create P gates for each qubit
    P_gates = []
    for qubit in range(1, n_qubits + 1):
        P_gates.append(create_pauli_z_diagonal(qubit, n_qubits))

    # Create all CZ gates (10 pairs for 5 qubits)
    cz_gates = []
    for i in range(1, n_qubits + 1):
        for j in range(i + 1, n_qubits + 1):
            cz_gates.append(create_cz_diagonal(i, j, n_qubits))

    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 5):
                        # For each CZ gate, use powers 1 and 2
                        for cz_powers in np.ndindex(*([2] * 10)):  # 10 CZ gates
                            result_diag = (P_gates[0] ** i) * (P_gates[1] ** j) * \
                                          (P_gates[2] ** k) * (P_gates[3] ** l) * (P_gates[4] ** m)

                            for idx, power in enumerate(cz_powers):
                                result_diag *= (cz_gates[idx] ** (power + 1))

                            if np.allclose(result_diag.imag, 0):
                                diagonal_group.append(result_diag.real)

    return np.array(diagonal_group)


def six_qubit_diagonal_group():
    """Generate 6-qubit real diagonal Clifford group"""
    n_qubits = 6

    # Create P gates for each qubit
    P_gates = []
    for qubit in range(1, n_qubits + 1):
        P_gates.append(create_pauli_z_diagonal(qubit, n_qubits))

    # Create all CZ gates (15 pairs for 6 qubits)
    cz_gates = []
    for i in range(1, n_qubits + 1):
        for j in range(i + 1, n_qubits + 1):
            cz_gates.append(create_cz_diagonal(i, j, n_qubits))

    diagonal_group = []
    for i in range(1, 5):
        for j in range(1, 5):
            for k in range(1, 5):
                for l in range(1, 5):
                    for m in range(1, 5):
                        for n in range(1, 5):
                            # For each CZ gate, use powers 1 and 2
                            for cz_powers in np.ndindex(*([2] * 15)):  # 15 CZ gates
                                result_diag = (P_gates[0] ** i) * (P_gates[1] ** j) * \
                                              (P_gates[2] ** k) * (P_gates[3] ** l) * \
                                              (P_gates[4] ** m) * (P_gates[5] ** n)

                                for idx, power in enumerate(cz_powers):
                                    result_diag *= (cz_gates[idx] ** (power + 1))

                                if np.allclose(result_diag.imag, 0):
                                    diagonal_group.append(result_diag.real)

    return np.array(diagonal_group)


if __name__ == '__main__':

    # 1-qubit
    # group1 = one_qubit_diagonal_group()
    # np.savetxt('real_diagonal_cliffords/one_qubit_real_diagonal_clifford_groups.txt', group1)
    # print(f"Number of Clifford elements: {group1.shape[0]}")

    # # 2-qubit
    # group2 = two_qubit_diagonal_group()
    # np.savetxt('real_diagonal_cliffords/two_qubit_real_diagonal_clifford_groups.txt', group2)
    # print("Number of 2-qubit Clifford elements:", group2.shape[0])


    # # 3-qubit
    # group3 = three_qubit_diagonal_group()
    # np.savetxt('real_diagonal_cliffords/three_qubit_real_diagonal_clifford_groups.txt', group3)
    # print("Number of 3-qubit Clifford elements:", group3.shape[0])

    # # 4-qubit
    # group4 = four_qubit_diagonal_group()
    # np.savetxt('real_diagonal_cliffords/four_qubit_real_diagonal_clifford_groups.txt', group4)
    # print("Number of 4-qubit Clifford elements:", group4.shape[0])

    # # 5-qubit
    # group5 = five_qubit_diagonal_group()
    # np.savetxt('five_qubit_real_diagonal_clifford_groups.txt', group5)
    # print("Number of 5-qubit Clifford elements:", group5.shape[0])
    
    # 6-qubit
    group6 = six_qubit_diagonal_group()
    np.savetxt('six_qubit_real_diagonal_clifford_groups.txt', group6)
    print("Number of 6-qubit Clifford elements:", group6.shape[0])
