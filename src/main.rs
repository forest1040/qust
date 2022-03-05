//use std::collections::HashMap;
use std::fmt;

pub struct QuantumState {
    dim: u64,
    qubit_count: u32,
    //is_state_vector: bool,
    //classical_register: HashMap<String, i32>,
    //device_number: u32,
}

impl QuantumState {
    pub fn new(qubit_count: u32) -> QuantumState {
        let dim = 1 << qubit_count;
        QuantumState {
            qubit_count,
            dim,
            //device_number: 0,
        }
    }
}

impl fmt::Display for QuantumState {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = format!(
            r"* Qubit Count : {}
* Dimension   : {}
",
            self.qubit_count, self.dim
        );
        write!(f, "{}", s)
    }
}
fn main() {
    let n = 4;
    let state = QuantumState::new(n);
    //println!("Hello, world!");
    println!("{}", state);
}
