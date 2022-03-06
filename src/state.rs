use crate::random::Xor128;
use crate::stat_ops_probability;
use num_complex::Complex64;
use std::fmt;

//use std::collections::HashMap;

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
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
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
        self.state_vector[0] = Complex64::new(0.0, 0.0);
        self.state_vector[comp_basis] = Complex64::new(1.0, 0.0);
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
            //     state[index] /= norm;
            self.state_vector[i] /= norm;
        }
    }

    pub fn get_zero_probability(&self, target_qubit_index: usize) -> f64 {
        // TODO: validate
        // if (target_qubit_index >= this->qubit_count) {
        //     throw std::invalid_argument("qubit index >= num_qubit");
        // }
        stat_ops_probability::M0_prob(target_qubit_index, &self.get_vector(), self.dim)
    }

    pub fn get_vector(&self) -> &Vec<Complex64> {
        &self.state_vector
    }

    pub fn load(&mut self, state: Vec<Complex64>) {
        // TODO: validate
        // if (_state.size() != _dim) {
        //     throw std::invalid_argument("inconsistent dim");
        // }
        self.state_vector = state
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
