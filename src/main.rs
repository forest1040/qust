use num_complex::Complex64;
use qust::state::QuantumState;
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

    state.set_zero_state();
    let vec = state.get_vector();
    println!("{:?}", vec);

    let mut newState = Vec::new();
    newState.push(Complex64 { re: 0.5, im: 0.0 });
    newState.push(Complex64 { re: 0.5, im: 0.0 });
    newState.push(Complex64 { re: 0.5, im: 0.0 });
    newState.push(Complex64 { re: 0.5, im: 0.0 });

    state.load(newState);
    println!("{}", state);

    let n = 5;
    let mut state = QuantumState::new(n);
    state.set_haar_random_state(seed);
    let index = 3;
    let zero_probability = state.get_zero_probability(index);
    println!("prob_meas_3rd : {}", zero_probability);
}
