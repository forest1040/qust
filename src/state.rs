//use crate::binary_search::BinarySearch;
use crate::qsort;
use crate::random::Xor128;
use crate::stat_ops_probability;
//use num_complex::{Complex, Complex64};
use num_complex::Complex64;
//use rand::{seq::IteratorRandom, thread_rng};
//use rand::Rng;
use std::fmt;

use rand::distributions::WeightedIndex;
use rand::prelude::{thread_rng, Distribution};

//use std::collections::HashMap;
//use ndarray::{Array, Vector};
//use ndarray::{Array, Ix1};

pub struct QuantumState {
    dim: u64,
    qubit_count: u32,
    //classical_register: HashMap<String, i32>,
    state_vector: Vec<Complex64>,
}

impl QuantumState {
    pub fn new(qubit_count: u32) -> QuantumState {
        let dim = 1 << qubit_count;
        let state_vector = QuantumState::create_init_state(dim);
        QuantumState {
            qubit_count,
            dim,
            state_vector,
        }
    }

    fn create_init_state(dim: u64) -> Vec<Complex64> {
        //let mut state_vector = Vec::new();
        let mut state_vector = Vec::with_capacity(dim as usize);
        for i in 0..dim {
            let v = if i == 0 {
                Complex64::new(1., 0.)
            } else {
                Complex64::new(0., 0.)
            };
            state_vector.push(v);
        }
        state_vector
    }

    fn initialize_quantum_state(&mut self, dim: u64) {
        self.state_vector = QuantumState::create_init_state(dim);
    }

    pub fn set_zero_state(&mut self) {
        self.initialize_quantum_state(self.dim);
    }

    pub fn set_computational_basis(&mut self, comp_basis: usize) {
        // TODO: validate
        // if (comp_basis >= (ITYPE)(1ULL << this->qubit_count)) {
        //     throw std::invalid_argument("basis index >= 2^n");
        // }
        self.set_zero_state();
        self.state_vector[0] = Complex64::new(0., 0.);
        self.state_vector[comp_basis] = Complex64::new(1., 0.);
    }

    pub fn set_haar_random_state(&mut self, seed: i64) {
        self.initialize_haar_random_state_with_seed(seed)
    }

    fn initialize_haar_random_state_with_seed(&mut self, seed: i64) {
        let mut xor = Xor128::from_seed(seed);
        let mut norm = 0.;
        for i in 0..self.dim as usize {
            let r1 = xor.random_normal();
            let r2 = xor.random_normal();
            //state[index] = r1 + 1.i * r2;
            self.state_vector[i] = Complex64::new(r1, r2);
            // TODO: self.state_vector[i].norm()
            norm += r1 * r1 + r2 * r2;
        }
        let norm = norm.sqrt();
        for i in 0..self.dim as usize {
            // state[index] /= norm;
            self.state_vector[i] /= norm;
        }
    }

    pub fn get_zero_probability(&self, target_qubit_index: usize) -> f64 {
        // TODO: validate
        // if (target_qubit_index >= this->qubit_count) {
        //     throw std::invalid_argument("qubit index >= num_qubit");
        // }
        stat_ops_probability::m0_prob(target_qubit_index, &self.get_vector(), self.dim)
    }

    pub fn sampling(&self, sampling_count: usize) -> Vec<usize> {
        let mut sum = 0.;
        let mut weights = Vec::with_capacity(self.dim as usize + 1);
        weights.push(0.);
        for i in 0..self.dim as usize {
            //sum += norm(ptr[i]);
            sum += self.state_vector[i].norm();
            weights.push(sum);
        }
        let dist = WeightedIndex::new(&weights).unwrap();
        let mut rng = thread_rng();
        let mut result: Vec<usize> = Vec::with_capacity(sampling_count);
        for _ in 0..sampling_count {
            result.push(dist.sample(&mut rng) - 1);
        }
        result
    }

    // pub fn inner_product(
    //     state_bra: &Array<Complex64, Ix1>,
    //     state_ket: &Array<Complex64, Ix1>,
    // ) -> Complexf64 {
    //     let vs = state_ket.iter().map(|x| x.conj()).collect();
    //     let result = state_bra.dot(&vs);
    //     //result
    // }

    pub fn inner_product(state_bra: QuantumState, state_ket: QuantumState) -> Complex64 {
        // double real_sum = 0.;
        // double imag_sum = 0.;
        // ITYPE index;
        // for (index = 0; index < dim; ++index) {
        //     CTYPE value;
        //     value += conj(state_bra[index]) * state_ket[index];
        //     real_sum += _creal(value);
        //     imag_sum += _cimag(value);
        // }
        // return real_sum + 1.i * imag_sum;
        let mut real_sum = 0.;
        let mut imag_sum = 0.;
        for i in 0..state_bra.dim as usize {
            let value = state_bra.state_vector[i].conj() * state_ket.state_vector[i];
            real_sum += value.re;
            imag_sum += value.im;
        }
        Complex64::new(real_sum, imag_sum)
    }

    pub fn get_vector(&self) -> &Vec<Complex64> {
        &self.state_vector
    }

    pub fn get_dim(&self) -> u64 {
        self.dim
    }

    pub fn load(&mut self, state: Vec<Complex64>) {
        // TODO: validate
        // if (_state.size() != _dim) {
        //     throw std::invalid_argument("inconsistent dim");
        // }
        self.state_vector = state
    }

    pub fn single_qubit_dense_matrix_gate(
        &mut self,
        target_qubit_index: usize,
        matrix: Vec<Complex64>,
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
            let cval_0 = self.state_vector[basis_0];
            let cval_1 = self.state_vector[basis_1];
            self.state_vector[basis_0] = matrix[0] * cval_0 + matrix[1] * cval_1;
            self.state_vector[basis_1] = matrix[2] * cval_0 + matrix[3] * cval_1;
        }
    }

    pub fn single_qubit_control_single_qubit_dense_matrix_gate(
        &mut self,
        control_qubit_index: u64,
        control_value: u64,
        target_qubit_index: u64,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        let loop_dim = dim / 4;
        let target_mask = 1 << target_qubit_index;
        let control_mask = 1 << control_qubit_index;
        let min_qubit_index = QuantumState::get_min_ui(control_qubit_index, target_qubit_index);
        let max_qubit_index = QuantumState::get_max_ui(control_qubit_index, target_qubit_index);
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
                let cval0 = self.state_vector[basis_index_0 as usize];
                let cval1 = self.state_vector[basis_index_1 as usize];
                // set values
                self.state_vector[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
            }
        } else if control_qubit_index == 0 {
            for state_index in 0..loop_dim {
                let basis_index_0 = (state_index & low_mask)
                    + ((state_index & mid_mask) << 1)
                    + ((state_index & high_mask) << 2)
                    + control_mask * control_value;
                let basis_index_1 = basis_index_0 + target_mask;

                // fetch values
                let cval0 = self.state_vector[basis_index_0 as usize];
                let cval1 = self.state_vector[basis_index_1 as usize];
                // set values
                self.state_vector[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
            }
        } else {
            for state_index in (0..loop_dim).step_by(2) {
                let basis_index_0 = (state_index & low_mask)
                    + ((state_index & mid_mask) << 1)
                    + ((state_index & high_mask) << 2)
                    + control_mask * control_value;
                let basis_index_1 = basis_index_0 + target_mask;

                // fetch values
                let cval0 = self.state_vector[basis_index_0 as usize];
                let cval1 = self.state_vector[basis_index_1 as usize];
                let cval2 = self.state_vector[(basis_index_0 + 1) as usize];
                let cval3 = self.state_vector[(basis_index_1 + 1) as usize];
                // set values
                self.state_vector[basis_index_0 as usize] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_index_1 as usize] = matrix[2] * cval0 + matrix[3] * cval1;
                self.state_vector[(basis_index_0 + 1) as usize] =
                    matrix[0] * cval2 + matrix[1] * cval3;
                self.state_vector[(basis_index_1 + 1) as usize] =
                    matrix[2] * cval2 + matrix[3] * cval3;
            }
        }
    }

    pub fn double_qubit_dense_matrix_gate_c(
        &mut self,
        target_qubit_index1: usize,
        target_qubit_index2: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        let min_qubit_index =
            QuantumState::get_min_ui(target_qubit_index1 as u64, target_qubit_index2 as u64);
        let max_qubit_index =
            QuantumState::get_max_ui(target_qubit_index1 as u64, target_qubit_index2 as u64);
        let min_qubit_mask = 1 << min_qubit_index;
        let max_qubit_mask = 1 << (max_qubit_index - 1);
        let low_mask = min_qubit_mask - 1;
        let mid_mask = (max_qubit_mask - 1) ^ low_mask;
        let high_mask = !(max_qubit_mask - 1);
        let target_mask1 = 1 << target_qubit_index1;
        let target_mask2 = 1 << target_qubit_index2;
        // loop variables
        let loop_dim = dim / 4;
        for state_index in 0..loop_dim {
            // create index
            let basis_0 = (state_index & low_mask)
                + ((state_index & mid_mask) << 1)
                + ((state_index & high_mask) << 2);
            let basis_0 = basis_0 as usize;

            // gather index
            let basis_1 = basis_0 + target_mask1;
            let basis_2 = basis_0 + target_mask2;
            let basis_3 = basis_1 + target_mask2;

            // fetch values
            let cval_0 = self.state_vector[basis_0];
            let cval_1 = self.state_vector[basis_1];
            let cval_2 = self.state_vector[basis_2];
            let cval_3 = self.state_vector[basis_3];

            // set values
            self.state_vector[basis_0] =
                matrix[0] * cval_0 + matrix[1] * cval_1 + matrix[2] * cval_2 + matrix[3] * cval_3;
            self.state_vector[basis_1] =
                matrix[4] * cval_0 + matrix[5] * cval_1 + matrix[6] * cval_2 + matrix[7] * cval_3;
            self.state_vector[basis_2] =
                matrix[8] * cval_0 + matrix[9] * cval_1 + matrix[10] * cval_2 + matrix[11] * cval_3;
            self.state_vector[basis_3] = matrix[12] * cval_0
                + matrix[13] * cval_1
                + matrix[14] * cval_2
                + matrix[15] * cval_3;
        }
    }

    pub fn multi_qubit_dense_matrix_gate_single(
        &mut self,
        target_qubit_index_list: &Vec<usize>,
        target_qubit_index_count: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        let mut sort_array = Vec::with_capacity(64) as Vec<usize>;
        let mut mask_array = Vec::with_capacity(64) as Vec<usize>;
        self.create_shift_mask_list_from_list_buf(
            target_qubit_index_list,
            target_qubit_index_count,
            &mut sort_array,
            &mut mask_array,
        );

        // matrix dim, mask, buffer
        let matrix_dim = 1 << target_qubit_index_count;
        let matrix_mask_list =
            self.create_matrix_mask_list(target_qubit_index_list, target_qubit_index_count);
        // loop variables
        let loop_dim = dim >> target_qubit_index_count;
        let mut buffer = Vec::with_capacity(matrix_dim);
        for state_index in 0..loop_dim {
            // create base index
            let mut basis_0 = state_index as usize;
            for cursor in 0..target_qubit_index_count {
                basis_0 = (basis_0 & mask_array[cursor]) + ((basis_0 & (!mask_array[cursor])) << 1);
            }
            // compute matrix-vector multiply
            for y in 0..matrix_dim {
                buffer.push(Complex64 { re: 0., im: 0. });
                //buffer[y] = 0;
                for x in 0..matrix_dim {
                    buffer[y] += matrix[y * matrix_dim + x]
                        * self.state_vector[basis_0 ^ matrix_mask_list[x]];
                }
            }
            // set result
            for y in 0..matrix_dim {
                self.state_vector[basis_0 ^ matrix_mask_list[y]] = buffer[y];
            }
        }
    }

    fn create_matrix_mask_list(
        &self,
        qubit_index_list: &Vec<usize>,
        qubit_index_count: usize,
    ) -> Vec<usize> {
        let matrix_dim = 1 << qubit_index_count;
        let mut mask_list = Vec::with_capacity(matrix_dim);
        for cursor in 0..matrix_dim {
            mask_list.push(0);
            for bit_cursor in 0..qubit_index_count {
                if (cursor >> bit_cursor) % 2 == 0 {
                    let bit_index = qubit_index_list[bit_cursor];
                    mask_list[cursor] ^= 1 << bit_index;
                }
            }
        }
        mask_list
    }

    pub fn multi_qubit_control_single_qubit_dense_matrix_gate(
        &mut self,
        control_qubit_index_list: &Vec<usize>,
        control_value_list: &Vec<usize>,
        control_qubit_index_count: usize,
        target_qubit_index: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        if control_qubit_index_count == 1 {
            self.single_qubit_control_single_qubit_dense_matrix_gate(
                control_qubit_index_list[0] as u64,
                control_value_list[0] as u64,
                target_qubit_index as u64,
                matrix,
                dim,
            );
            return;
        }

        let mut sort_array = Vec::with_capacity(64) as Vec<usize>;
        let mut mask_array = Vec::with_capacity(64) as Vec<usize>;
        self.create_shift_mask_list_from_list_and_value_buf(
            control_qubit_index_list,
            control_qubit_index_count,
            target_qubit_index,
            &mut sort_array,
            &mut mask_array,
        );
        let target_mask = 1 << target_qubit_index;
        let control_mask = self.create_control_mask(
            control_qubit_index_list,
            control_value_list,
            control_qubit_index_count,
        );
        let insert_index_list_count = control_qubit_index_count + 1;
        let loop_dim = dim >> insert_index_list_count;

        if target_qubit_index == 0 {
            for state_index in 0..loop_dim {
                // create base index
                let mut basis_0 = state_index as usize;
                for cursor in 0..insert_index_list_count {
                    basis_0 =
                        (basis_0 & mask_array[cursor]) + ((basis_0 & (!mask_array[cursor])) << 1);
                }
                basis_0 += control_mask;

                // fetch values
                let cval0 = self.state_vector[basis_0];
                let cval1 = self.state_vector[basis_0 + 1];
                // set values
                self.state_vector[basis_0] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_0 + 1] = matrix[2] * cval0 + matrix[3] * cval1;
            }
        } else if sort_array[0] == 0 {
            for state_index in 0..loop_dim {
                // create base index
                let mut basis_0 = state_index as usize;
                for cursor in 0..insert_index_list_count {
                    basis_0 =
                        (basis_0 & mask_array[cursor]) + ((basis_0 & (!mask_array[cursor])) << 1);
                }
                basis_0 += control_mask;
                let basis_1 = basis_0 + target_mask;

                // fetch values
                let cval0 = self.state_vector[basis_0];
                let cval1 = self.state_vector[basis_1];
                // set values
                self.state_vector[basis_0] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_1] = matrix[2] * cval0 + matrix[3] * cval1;
            }
        } else {
            for state_index in (0..loop_dim).step_by(2) {
                // create base index
                let mut basis_0 = state_index as usize;
                for cursor in 0..insert_index_list_count {
                    basis_0 =
                        (basis_0 & mask_array[cursor]) + ((basis_0 & (!mask_array[cursor])) << 1);
                }
                basis_0 += control_mask;
                let basis_1 = basis_0 + target_mask;

                // fetch values
                let cval0 = self.state_vector[basis_0];
                let cval1 = self.state_vector[basis_1];
                let cval2 = self.state_vector[basis_0 + 1];
                let cval3 = self.state_vector[basis_1 + 1];
                // set values
                self.state_vector[basis_0] = matrix[0] * cval0 + matrix[1] * cval1;
                self.state_vector[basis_1] = matrix[2] * cval0 + matrix[3] * cval1;
                self.state_vector[basis_0 + 1] = matrix[0] * cval2 + matrix[1] * cval3;
                self.state_vector[basis_1 + 1] = matrix[2] * cval2 + matrix[3] * cval3;
            }
        }
    }

    pub fn single_qubit_control_multi_qubit_dense_matrix_gate(
        &mut self,
        control_qubit_index: usize,
        control_value: usize,
        target_qubit_index_list: &Vec<usize>,
        target_qubit_index_count: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        // matrix dim, mask, buffer
        let matrix_dim = 1 << target_qubit_index_count;
        let matrix_mask_list =
            self.create_matrix_mask_list(target_qubit_index_list, target_qubit_index_count);
        //CTYPE* buffer = (CTYPE*)malloc((size_t)(sizeof(CTYPE) * matrix_dim));
        let mut buffer = Vec::with_capacity(matrix_dim);

        // insert list
        let insert_index_count = target_qubit_index_count + 1;
        let sorted_insert_index_list = self.create_sorted_ui_list_value(
            target_qubit_index_list,
            target_qubit_index_count,
            control_qubit_index,
        );

        // control mask
        let control_mask = (1 << control_qubit_index) * control_value;

        // loop varaibles
        let loop_dim = dim >> insert_index_count;
        for state_index in 0..loop_dim {
            // create base index
            let mut basis_0 = state_index as usize;
            for cursor in 0..insert_index_count {
                let insert_index = sorted_insert_index_list[cursor];
                basis_0 = stat_ops_probability::insert_zero_to_basis_index(
                    basis_0,
                    1 << insert_index,
                    insert_index,
                );
            }

            // flip control
            basis_0 ^= control_mask;

            // compute matrix mul
            for y in 0..matrix_dim {
                buffer.push(Complex64 { re: 0., im: 0. });
                //buffer[y] = 0;
                for x in 0..matrix_dim {
                    buffer[y] += matrix[y * matrix_dim + x]
                        * self.state_vector[basis_0 ^ matrix_mask_list[x]];
                }
            }

            // set result
            for y in 0..matrix_dim {
                self.state_vector[basis_0 ^ matrix_mask_list[y]] = buffer[y];
            }
        }
    }

    pub fn multi_qubit_control_multi_qubit_dense_matrix_gate(
        &mut self,
        control_qubit_index_list: &Vec<usize>,
        control_value_list: &Vec<usize>,
        control_qubit_index_count: usize,
        target_qubit_index_list: &Vec<usize>,
        target_qubit_index_count: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        // matrix dim, mask, buffer
        let matrix_dim = 1 << target_qubit_index_count;
        let matrix_mask_list =
            self.create_matrix_mask_list(target_qubit_index_list, target_qubit_index_count);
        //CTYPE* buffer = (CTYPE*)malloc((size_t)(sizeof(CTYPE) * matrix_dim));
        let mut buffer = Vec::with_capacity(matrix_dim);

        // insert index
        let insert_index_count = target_qubit_index_count + control_qubit_index_count;
        let sorted_insert_index_list = self.create_sorted_ui_list_list(
            target_qubit_index_list,
            target_qubit_index_count,
            control_qubit_index_list,
            control_qubit_index_count,
        );

        // control mask
        let control_mask = self.create_control_mask(
            control_qubit_index_list,
            control_value_list,
            control_qubit_index_count,
        );

        // loop varaibles
        let loop_dim = dim >> (target_qubit_index_count + control_qubit_index_count);

        for state_index in 0..loop_dim {
            // create base index
            let mut basis_0 = state_index as usize;
            for cursor in 0..insert_index_count {
                let insert_index = sorted_insert_index_list[cursor];
                basis_0 = stat_ops_probability::insert_zero_to_basis_index(
                    basis_0,
                    1 << insert_index,
                    insert_index,
                );
            }

            // flip control masks
            basis_0 ^= control_mask;

            // compute matrix mul
            for y in 0..matrix_dim {
                buffer.push(Complex64 { re: 0., im: 0. });
                //buffer[y] = 0;
                for x in 0..matrix_dim {
                    buffer[y] += matrix[y * matrix_dim + x]
                        * self.state_vector[basis_0 ^ matrix_mask_list[x]];
                }
            }

            // set result
            for y in 0..matrix_dim {
                self.state_vector[basis_0 ^ matrix_mask_list[y]] = buffer[y];
            }
        }
    }

    fn create_shift_mask_list_from_list_and_value_buf(
        &self,
        array: &Vec<usize>,
        count: usize,
        target: usize,
        dst_array: &mut Vec<usize>,
        dst_mask: &mut Vec<usize>,
    ) {
        // TODO: memcpy() to copy_from()
        // memcpy(dst_array, array, sizeof(UINT) * count);
        let size = count + 1;
        for i in 0..array.len() {
            dst_array[i] = array[i];
        }
        dst_array[count] = target;
        qsort::quick_sort(dst_array, 0, size);
        for i in 0..size {
            dst_mask[i] = (1 << dst_array[i]) - 1;
        }
    }

    fn create_shift_mask_list_from_list_buf(
        &self,
        array: &Vec<usize>,
        count: usize,
        dst_array: &mut Vec<usize>,
        dst_mask: &mut Vec<usize>,
    ) {
        // TODO: memcpy() to copy_from()
        //memcpy(dst_array, array, sizeof(UINT) * count);
        for i in 0..array.len() {
            dst_array[i] = array[i];
        }
        qsort::quick_sort(dst_array, 0, count);
        for i in 0..count {
            dst_mask[i] = (1 << dst_array[i]) - 1;
        }
    }

    fn create_sorted_ui_list_value(
        &self,
        array: &Vec<usize>,
        size: usize,
        value: usize,
    ) -> Vec<usize> {
        // TODO: memcpy() to copy_from()
        // UINT* new_array = (UINT*)calloc(size + 1, sizeof(UINT));
        // memcpy(new_array, array, size * sizeof(UINT));
        let count = size + 1;
        let mut new_array = Vec::with_capacity(count);
        for i in 0..count {
            new_array.push(0);
        }
        new_array[size] = value;
        //sort_ui(new_array, size + 1);
        qsort::quick_sort(&mut new_array, 0, count);
        new_array
    }

    fn create_sorted_ui_list_list(
        &self,
        array1: &Vec<usize>,
        size1: usize,
        array2: &Vec<usize>,
        size2: usize,
    ) -> Vec<usize> {
        // UINT* new_array = (UINT*)calloc(size1 + size2, sizeof(UINT));
        // memcpy(new_array, array1, size1 * sizeof(UINT));
        // memcpy(new_array + size1, array2, size2 * sizeof(UINT));
        let count = size1 + size2;
        let mut new_array = Vec::with_capacity(count);
        for i in 0..count {
            new_array.push(0);
        }
        qsort::quick_sort(&mut new_array, 0, count);
        new_array
    }

    fn create_control_mask(
        &self,
        qubit_index_list: &Vec<usize>,
        value_list: &Vec<usize>,
        size: usize,
    ) -> usize {
        let mut mask = 0;
        for cursor in 0..size {
            mask ^= (1 << qubit_index_list[cursor]) * value_list[cursor];
        }
        mask
    }

    pub fn multi_qubit_dense_matrix_gate(
        &mut self,
        target_qubit_index_list: &Vec<usize>,
        target_qubit_index_count: usize,
        matrix: Vec<Complex64>,
        dim: u64,
    ) {
        if target_qubit_index_count == 1 {
            self.single_qubit_dense_matrix_gate(target_qubit_index_list[0], matrix, dim);
        } else if target_qubit_index_count == 2 {
            self.double_qubit_dense_matrix_gate_c(
                target_qubit_index_list[0],
                target_qubit_index_list[1],
                matrix,
                dim,
            );
        } else {
            self.multi_qubit_dense_matrix_gate_single(
                target_qubit_index_list,
                target_qubit_index_count,
                matrix,
                dim,
            )
        }
    }

    pub fn single_qubit_Pauli_gate(
        &mut self,
        target_qubit_index: usize,
        pauli_operator_type: usize,
        dim: u64,
    ) {
        // TODO: enum化
        // match pauli_operator_type {
        //     1 => self.X_gate(target_qubit_index, dim),
        //     2 => self.Y_gate(target_qubit_index, dim),
        //     3 => self.Z_gate(target_qubit_index, dim),
        // }
        if pauli_operator_type == 1 {
            self.X_gate(target_qubit_index, dim);
        } else if pauli_operator_type == 2 {
            self.Y_gate(target_qubit_index, dim);
        } else {
            self.Z_gate(target_qubit_index, dim);
        }
    }

    pub fn X_gate(&mut self, target_qubit_index: usize, dim: u64) {
        let loop_dim = dim / 2;
        let mask = 1 << target_qubit_index;
        let mask_low = mask - 1;
        let mask_high = !mask_low;
        if target_qubit_index == 0 {
            for basis_index in (0..dim as usize).step_by(2) {
                let temp = self.state_vector[basis_index];
                self.state_vector[basis_index] = self.state_vector[basis_index + 1];
                self.state_vector[basis_index + 1] = temp;
            }
        } else {
            for state_index in (0..loop_dim as usize).step_by(2) {
                let basis_index_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);
                let basis_index_1 = basis_index_0 + mask;
                let temp0 = self.state_vector[basis_index_0];
                let temp1 = self.state_vector[basis_index_0 + 1];
                self.state_vector[basis_index_0] = self.state_vector[basis_index_1];
                self.state_vector[basis_index_0 + 1] = self.state_vector[basis_index_1 + 1];
                self.state_vector[basis_index_1] = temp0;
                self.state_vector[basis_index_1 + 1] = temp1;
            }
        }
    }

    pub fn Y_gate(&mut self, target_qubit_index: usize, dim: u64) {
        let loop_dim = dim / 2;
        let mask = 1 << target_qubit_index;
        let mask_low = mask - 1;
        let mask_high = !mask_low;
        let imag = Complex64 { re: 0., im: 1. };
        if target_qubit_index == 0 {
            for basis_index in (0..dim as usize).step_by(2) {
                let temp0 = self.state_vector[basis_index];
                self.state_vector[basis_index] = -imag * self.state_vector[basis_index + 1];
                self.state_vector[basis_index + 1] = imag * temp0;
            }
        } else {
            for state_index in (0..loop_dim as usize).step_by(2) {
                let basis_index_0 = (state_index & mask_low) + ((state_index & mask_high) << 1);
                let basis_index_1 = basis_index_0 + mask;
                let temp0 = self.state_vector[basis_index_0];
                let temp1 = self.state_vector[basis_index_0 + 1];
                self.state_vector[basis_index_0] = -imag * self.state_vector[basis_index_1];
                self.state_vector[basis_index_0 + 1] = -imag * self.state_vector[basis_index_1 + 1];
                self.state_vector[basis_index_1] = imag * temp0;
                self.state_vector[basis_index_1 + 1] = imag * temp1;
            }
        }
    }

    pub fn Z_gate(&mut self, target_qubit_index: usize, dim: u64) {
        let loop_dim = dim / 2;
        let mask = 1 << target_qubit_index;
        let mask_low = mask - 1;
        let mask_high = !mask_low;
        if target_qubit_index == 0 {
            for state_index in (1..loop_dim as usize).step_by(2) {
                self.state_vector[state_index] *= Complex64 { re: -1., im: 0. };
            }
        } else {
            for state_index in (0..loop_dim as usize).step_by(2) {
                let basis_index = (state_index & mask_low) + ((state_index & mask_high) << 1);
                self.state_vector[basis_index] *= Complex64 { re: -1., im: 0. };
                self.state_vector[basis_index + 1] *= Complex64 { re: -1., im: 0. };
            }
        }
    }

    pub fn single_qubit_Pauli_rotation_gate(
        &mut self,
        target_qubit_index: usize,
        pauli_operator_type: usize,
        angle: f64,
        dim: u64,
    ) {
        // TODO: enum化
        if pauli_operator_type == 1 {
            self.RX_gate(target_qubit_index, dim);
        } else if pauli_operator_type == 2 {
            self.RY_gate(target_qubit_index, dim);
        } else {
            self.RZ_gate(target_qubit_index, dim);
        }
    }

    // void RX_gate(UINT target_qubit_index, double angle, CTYPE* state, ITYPE dim) {
    //     CTYPE matrix[4];
    //     double c, s;
    //     c = cos(angle / 2);
    //     s = sin(angle / 2);
    //     matrix[0] = c;
    //     matrix[1] = 1.i * s;
    //     matrix[2] = 1.i * s;
    //     matrix[3] = c;
    //     single_qubit_dense_matrix_gate(target_qubit_index, matrix, state, dim);
    // }

    // void RY_gate(UINT target_qubit_index, double angle, CTYPE* state, ITYPE dim) {
    //     CTYPE matrix[4];
    //     double c, s;
    //     c = cos(angle / 2);
    //     s = sin(angle / 2);
    //     matrix[0] = c;
    //     matrix[1] = s;
    //     matrix[2] = -s;
    //     matrix[3] = c;
    //     single_qubit_dense_matrix_gate(target_qubit_index, matrix, state, dim);
    // }

    // void RZ_gate(UINT target_qubit_index, double angle, CTYPE* state, ITYPE dim) {
    //     CTYPE diagonal_matrix[2];
    //     double c, s;
    //     c = cos(angle / 2);
    //     s = sin(angle / 2);
    //     diagonal_matrix[0] = c + 1.i * s;
    //     diagonal_matrix[1] = c - 1.i * s;
    //     single_qubit_diagonal_matrix_gate(
    //         target_qubit_index, diagonal_matrix, state, dim);
    // }

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
}

impl fmt::Display for QuantumState {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // TODO: state_vector format
        let sv: Vec<String> = self.state_vector.iter().map(|x| x.to_string()).collect();
        let svs = sv.join("\n");
        let tmp = format!(
            r"* Qubit Count  : {}
* Dimension    : {}
* State vector :
{}
",
            self.qubit_count, self.dim, svs
        );
        write!(f, "{}", tmp)
    }
}
