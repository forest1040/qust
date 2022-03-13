use crate::state::QuantumState;
use num_complex::Complex64;

pub fn single_qubit_dense_matrix_gate(
    target_qubit_index: usize,
    matrix: Vec<Complex64>,
    state: &QuantumState,
    dim: u64,
) {
    let loop_dim = dim / 2;
    let mask = 1 << target_qubit_index;
    let mask_low = mask - 1;
    let mask_high = !mask_low;
    for state_index in 0..loop_dim {
        let basis_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);
        let basis_1 = basis_0 + mask;
        let basis_0 = basis_0 as usize;
        let basis_1 = basis_1 as usize;
        let cval_0 = state.get_vector()[basis_0];
        let cval_1 = state.get_vector()[basis_1];
        state.get_vector()[basis_0] = matrix[0] * cval_0 + matrix[1] * cval_1;
        state.get_vector()[basis_1] = matrix[2] * cval_0 + matrix[3] * cval_1;
    }
}
