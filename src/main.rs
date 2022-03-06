//use std::collections::HashMap;
use num_complex::Complex64;
use srand::{Rand, RngSource};
use std::fmt;

struct Xor128 {
    x: u64,
    y: u64,
    z: u64,
    w: u64,
}

impl Xor128 {
    pub fn from_seed(seed: i64) -> Xor128 {
        let mut r: Rand<_> = Rand::new(RngSource::new(seed));
        let ignore_first = 40;
        let mut res = Xor128 {
            x: r.uint64(),
            y: r.uint64(),
            z: r.uint64(),
            w: r.uint64(),
        };
        for _ in 0..ignore_first {
            res.next();
        }
        res
    }

    pub fn next(&mut self) -> u64 {
        let t = self.x ^ (self.x << 11);
        self.x = self.y;
        self.y = self.z;
        self.z = self.w;
        self.w = (self.w ^ (self.w >> 19)) ^ (t ^ (t >> 8));
        self.w
    }

    pub fn random_uniform(&mut self) -> f64 {
        self.next() as f64 / std::u64::MAX as f64
    }

    pub fn random_normal(&mut self) -> f64 {
        // return sqrt(-1.0 * log(random_uniform(state))) *
        //     sin(2.0 * M_PI * random_uniform(state));
        (-1.0 * self.random_uniform().log(std::f64::consts::E)).sqrt()
            * (2.0 * std::f64::consts::PI * self.random_uniform()).sin()
    }
}

struct QuantumState {
    dim: u64,
    qubit_count: u32,
    //classical_register: HashMap<String, i32>,
    state_vector: Vec<Complex64>,
}

impl QuantumState {
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

    pub fn new(qubit_count: u32) -> QuantumState {
        let dim = 1 << qubit_count;
        let state_vector = QuantumState::create_init_state(dim);
        QuantumState {
            qubit_count,
            dim,
            state_vector,
        }
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
        for i in 0..self.dim {
            let r1 = xor.random_normal();
            let r2 = xor.random_normal();
            //state[index] = r1 + 1.i * r2;
            self.state_vector[i as usize] = Complex64::new(r1, r2);
            norm += r1 * r1 + r2 * r2;
        }
        let norm = norm.sqrt();
        for i in 0..self.dim {
            //     state[index] /= norm;
            self.state_vector[i as usize] /= norm;
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
fn main() {
    let n = 2;
    let mut state = QuantumState::new(n);
    println!("{}", state);

    // |00>に初期化
    state.set_zero_state();
    // 基底を二進数と見た時の整数値を入れて、その状態に初期化
    state.set_computational_basis(1);
    println!("{}", state);

    // シードを指定してランダムな初期状態を生成
    // (シードを省略した場合は時刻を用いてシードを決定します。)
    let seed = 0;
    state.set_haar_random_state(seed);
    println!("{}", state);
}
