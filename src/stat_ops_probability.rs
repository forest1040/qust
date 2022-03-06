use num_complex::Complex64;

pub fn M0_prob(target_qubit_index: usize, state: &Vec<Complex64>, dim: u64) -> f64 {
    let loop_dim = dim / 2;
    let mask: usize = 1 << target_qubit_index;
    let mut sum = 0.;
    for i in 0..loop_dim as usize {
        let basis_1 = insert_zero_to_basis_index(i, mask, target_qubit_index) ^ mask;
        //sum += pow(_cabs(state[basis_1]), 2);
        sum += state[basis_1].powi(2).re;
    }
    sum
}

#[inline]
pub fn insert_zero_to_basis_index(
    basis_index: usize,
    basis_mask: usize,
    qubit_index: usize,
) -> usize {
    // ITYPE temp_basis = (basis_index >> qubit_index) << (qubit_index + 1);
    // return temp_basis + basis_index % basis_mask;
    let temp_basis = (basis_index >> qubit_index) << (qubit_index + 1);
    temp_basis + basis_index % basis_mask
}
