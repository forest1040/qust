use crate::state::QuantumState;
use num_complex::Complex64;

pub fn single_qubit_dense_matrix_gate(
    target_qubit_index: usize,
    matrix: Vec<Complex64>,
    state: &mut QuantumState,
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

pub fn single_qubit_control_single_qubit_dense_matrix_gate(
    control_qubit_index: u64,
    control_value: u64,
    target_qubit_index: u64,
    matrix: Vec<Complex64>,
    state: &mut QuantumState,
    dim: u64,
) {
    let loop_dim = dim / 4;
    let target_mask = 1 << target_qubit_index;
    let control_mask = 1 << control_qubit_index;
    let min_qubit_index = get_min_ui(control_qubit_index, target_qubit_index);
    let max_qubit_index = get_max_ui(control_qubit_index, target_qubit_index);
    let min_qubit_mask = 1 << min_qubit_index;
    let max_qubit_mask = 1 << (max_qubit_index - 1);
    let low_mask = min_qubit_mask - 1;
    let mid_mask = (max_qubit_mask - 1) ^ low_mask;
    let high_mask = !(max_qubit_mask - 1);

    if target_qubit_index == 0 {
        for state_index in 0..loop_dim {
            let basis_index_0 = (state_index & low_mask)
                + ((state_index & mid_mask) << 1)
                + ((state_index & high_mask) << 2)
                + control_mask * control_value;
            let basis_index_1 = basis_index_0 + 1;

            // fetch values
            let cval0 = state.get_vector()[basis_index_0 as usize];
            let cval1 = state.get_vector()[basis_index_1 as usize];
            // set values
            state.get_vector()[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
            state.get_vector()[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
        }
    } else if control_qubit_index == 0 {
        for state_index in 0..loop_dim {
            let basis_index_0 = (state_index & low_mask)
                + ((state_index & mid_mask) << 1)
                + ((state_index & high_mask) << 2)
                + control_mask * control_value;
            let basis_index_1 = basis_index_0 + target_mask;

            // fetch values
            let cval0 = state.get_vector()[basis_index_0 as usize];
            let cval1 = state.get_vector()[basis_index_1 as usize];
            // set values
            state.get_vector()[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
            state.get_vector()[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
        }
    } else {
        for state_index in (0..loop_dim).step_by(2) {
            let basis_index_0 = (state_index & low_mask)
                + ((state_index & mid_mask) << 1)
                + ((state_index & high_mask) << 2)
                + control_mask * control_value;
            let basis_index_1 = basis_index_0 + target_mask;

            // fetch values
            let cval0 = state.get_vector()[basis_index_0 as usize];
            let cval1 = state.get_vector()[basis_index_1 as usize];
            let cval2 = state.get_vector()[(basis_index_0 + 1) as usize];
            let cval3 = state.get_vector()[(basis_index_1 + 1) as usize];
            // set values
            state.get_vector()[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
            state.get_vector()[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
            state.get_vector()[(basis_index_0 + 1) as usize] =
                matrix[0] * cval2 + matrix[1] * cval3;
            state.get_vector()[(basis_index_1 + 1) as usize] =
                matrix[2] * cval2 + matrix[3] * cval3;
        }
    }
}

#[inline]
fn get_min_ui(index_0: u64, index_1: u64) -> u64 {
    if index_0 < index_1 {
        index_0
    } else {
        index_1
    }
}

#[inline]
fn get_max_ui(index_0: u64, index_1: u64) -> u64 {
    if index_0 > index_1 {
        index_0
    } else {
        index_1
    }
}
