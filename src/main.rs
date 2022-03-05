//use std::collections::HashMap;
use num_complex::Complex64;
use std::fmt;

pub struct QuantumState {
    dim: u64,
    qubit_count: u32,
    //classical_register: HashMap<String, i32>,
    state_vector: Vec<Complex64>,
}

impl QuantumState {
    pub fn initialize_quantum_state(&mut self, dim: u64) {
        let mut state_vector = Vec::new();
        for i in 0..dim {
            let v = if i == 0 {
                Complex64::new(1.0, 0.0)
            } else {
                Complex64::new(0.0, 0.0)
            };
            state_vector.push(v);
        }
        self.state_vector = state_vector;
    }

    pub fn new(qubit_count: u32) -> QuantumState {
        let dim = 1 << qubit_count;
        let state_vector = Vec::with_capacity(dim as usize);
        let mut this = QuantumState {
            qubit_count,
            dim,
            state_vector,
        };
        this.initialize_quantum_state(dim);
        this
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
fn main() {
    let n = 2;
    let state = QuantumState::new(n);
    println!("{}", state);
}
