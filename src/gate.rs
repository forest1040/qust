//use std::collections::HashMap;
use crate::state::QuantumState;
//use crate::update_ops_matrix_dense_single;
use ndarray::prelude::*;
use num_complex::Complex64;
use std::fmt;

// TODO: constant.rs
const COMP_ZERO: Complex64 = Complex64 { re: 0., im: 0. };
const COMP_ONE: Complex64 = Complex64 { re: 1., im: 0. };
const COMP_I: Complex64 = Complex64 { re: 0., im: 1. };
const COMP_MINUS_ONE: Complex64 = Complex64 { re: -1., im: 0. };
const COMP_MINUS_I: Complex64 = Complex64 { re: 0., im: -1. };

const FLAG_COMMUTE_X: usize = 0x01;
const FLAG_COMMUTE_Y: usize = 0x02;
const FLAG_COMMUTE_Z: usize = 0x04;

// const FLAG_X_COMMUTE: usize = 0x01;
// const FLAG_Y_COMMUTE: usize = 0x02;
// const FLAG_Z_COMMUTE: usize = 0x04;

// enum MapType {
//     Basic,
//     Sequence,
//     Probabilistic,
//     CPTP,
//     Instrument,
// }

#[derive(Debug)]
enum GateMatrixType {
    DenseMatrix,
    // TODO
    // DiagonalMatrix,
    // SparseMatrix,
    PauliMatrix,
    // PermutationMatrix,
}

// enum SpecialFuncType {
//     None,
//     GateI,
//     GateX,
//     GateY,
//     GateZ,
//     GateSqrtX,
//     GateSqrtY,
//     GateSqrtXdag,
//     GateSqrtYdag,
//     GateRX,
//     GateRY,
//     GateRZ,
//     GateH,
//     GateS,
//     GateSdag,
//     GateT,
//     GateTdag,
//     GateP0,
//     GateP1,
//     GateCX,
//     GateCY,
//     GateCZ,
//     GateSWAP,
// }

// enum {
//     PAULI_ID_I = 0,
//     PAULI_ID_X = 1,
//     PAULI_ID_Y = 2,
//     PAULI_ID_Z = 3,
// }

enum PauliID {
    PAULI_ID_I,
    PAULI_ID_X,
    PAULI_ID_Y,
    PAULI_ID_Z,
}

impl PauliID {
    fn from_usize(n: usize) -> Option<PauliID> {
        match n {
            0 => Some(PauliID::PAULI_ID_I),
            1 => Some(PauliID::PAULI_ID_X),
            2 => Some(PauliID::PAULI_ID_Y),
            3 => Some(PauliID::PAULI_ID_Z),
            _ => None,
        }
    }
}

#[derive(Debug)]
pub struct QuantumGate {
    //map_type: MapType,
    matrix_type: GateMatrixType,
    name: String,
    //special_func_type: SpecialFuncType,
    target_qubit_index: Vec<usize>,
    target_qubit_commutation: Vec<usize>,
    control_qubit_index: Vec<usize>,
    control_qubit_value: Vec<usize>,
    //gate_property: usize,

    // ComplexMatrix _dense_matrix_element;
    dense_matrix_element: Array2<Complex64>,
    // ComplexVector _diagonal_matrix_element;
    // diagonal_matrix_element: Array1<Complex64>,
    // TODO: SparseComplexMatrix
    // SparseComplexMatrix _sparse_matrix_element;
    pauli_id: Vec<usize>,
    rotation_angle: f64,
    // TODO: parameter
    //parameter: HashMap<String, f64>,
}

impl QuantumGate {
    pub fn get_target_index_list(&self) -> &Vec<usize> {
        &self.target_qubit_index
    }
    pub fn get_control_index_list(&self) -> &Vec<usize> {
        &self.control_qubit_index
    }
    pub fn get_qubit_index_list(&self) -> Vec<usize> {
        // RVOになるはず
        let mut res = self.target_qubit_index.clone();
        for cqi in &self.control_qubit_index {
            res.push(*cqi);
        }
        res
    }

    // // TODO: reset_qubit_index_list
    // pub fn reset_qubit_index_list(&self, src: &Vec<usize>, dst: &Vec<usize>) {
    //     // if (src_list.size() != dst_list.size())
    //     //     throw std::invalid_argument("src.size() != dst.size()");
    //     // for (UINT index = 0; index < src_list.size(); ++index) {
    //     //     UINT src = src_list[index];
    //     //     UINT dst = dst_list[index];
    //     //     auto replace_func = [](std::vector<UINT>& vec, UINT src_ind,
    //     //                             UINT dst_ind) -> void {
    //     //         for (auto& v : vec)
    //     //             if (v == src_ind) v = dst_ind;
    //     //     };
    //     //     replace_func(_target_qubit_index, src, dst);
    //     //     replace_func(_control_qubit_index, src, dst);
    //     // }
    // }

    pub fn update_quantum_state(&self, state: &mut QuantumState) {
        // TODO: is_state_vector
        // TODO: _special_func_type
        // if (state.is_state_vector()) {
        //     if (_special_func_type == None)
        //         this->_update_state_vector_cpu_general(state);
        //     else
        //         this->_update_state_vector_cpu_special(state);
        // } else {
        //     this->_update_density_matrix_cpu_general(state);
        // }
        self.update_state_vector_cpu_general(state);
    }

    fn update_state_vector_cpu_general(&self, state: &mut QuantumState) {
        let (w, h) = {
            let shape = self.dense_matrix_element.shape();
            (shape[1], shape[0])
        };
        let vector2 = self
            .dense_matrix_element
            .clone()
            .into_shape(h * w)
            .unwrap()
            .to_vec();

        // TODO: vector2を参照渡しにする
        // TODO: Vecのsizeをわざわざ渡さなくてよい
        // TODO: usizeとu64を統一する

        match self.matrix_type {
            GateMatrixType::DenseMatrix => {
                // single qubit dense matrix gate
                if self.target_qubit_index.len() == 1 {
                    // no control qubit
                    if self.control_qubit_index.len() == 0 {
                        state.single_qubit_dense_matrix_gate(
                            self.target_qubit_index[0],
                            vector2,
                            state.get_dim(),
                        );
                    }
                    // single control qubit
                    else if self.control_qubit_index.len() == 1 {
                        state.single_qubit_control_single_qubit_dense_matrix_gate(
                            self.control_qubit_index[0] as u64,
                            self.control_qubit_value[0] as u64,
                            self.target_qubit_index[0] as u64,
                            vector2,
                            state.get_dim(),
                        );
                    }
                    // multiple control qubits
                    else {
                        if self.control_qubit_index.len() == 1 {
                            state.single_qubit_control_single_qubit_dense_matrix_gate(
                                self.control_qubit_index[0] as u64,
                                self.control_qubit_value[0] as u64,
                                self.target_qubit_index[0] as u64,
                                vector2,
                                state.get_dim(),
                            );
                        } else {
                            state.multi_qubit_control_single_qubit_dense_matrix_gate(
                                &self.control_qubit_index,
                                &self.control_qubit_value,
                                self.control_qubit_index.len(),
                                self.target_qubit_index[0],
                                vector2,
                                state.get_dim(),
                            );
                        }
                    }
                }
                // multi qubit dense matrix gate
                else {
                    // no control qubit
                    if self.control_qubit_index.len() == 0 {
                        state.multi_qubit_dense_matrix_gate(
                            &self.target_qubit_index,
                            self.target_qubit_index.len(),
                            vector2,
                            state.get_dim(),
                        );
                    }
                    // single control qubit
                    else if self.control_qubit_index.len() == 1 {
                        state.single_qubit_control_multi_qubit_dense_matrix_gate(
                            self.control_qubit_index[0],
                            self.control_qubit_value[0],
                            &self.target_qubit_index,
                            self.target_qubit_index.len(),
                            vector2,
                            state.get_dim(),
                        );
                    }
                    // multiple control qubit
                    else {
                        state.multi_qubit_control_multi_qubit_dense_matrix_gate(
                            &self.control_qubit_index,
                            &self.control_qubit_value,
                            self.control_qubit_index.len(),
                            &self.target_qubit_index,
                            self.target_qubit_index.len(),
                            vector2,
                            state.get_dim(),
                        );
                    }
                }
            }
            GateMatrixType::PauliMatrix => {
                if self.target_qubit_index.len() == 1 {
                    if self.rotation_angle.abs() < 1e-16 {
                        state.single_qubit_Pauli_gate(
                            self.target_qubit_index[0],
                            self.pauli_id[0],
                            state.get_dim(),
                        );
                    } else {
                        // invert
                        state.single_qubit_Pauli_rotation_gate(
                            self.target_qubit_index[0],
                            self.pauli_id[0],
                            -1. * self.rotation_angle,
                            state.get_dim(),
                        );
                    }
                } else {
                    if self.rotation_angle.abs() < 1e-16 {
                        state.multi_qubit_Pauli_gate_partial_list(
                            &self.target_qubit_index,
                            &self.pauli_id,
                            self.target_qubit_index.len(),
                            state.get_dim() as usize,
                        );
                    } else {
                        // invert
                        state.multi_qubit_Pauli_rotation_gate_partial_list(
                            &self.target_qubit_index,
                            &self.pauli_id,
                            self.target_qubit_index.len(),
                            self.rotation_angle,
                            state.get_dim() as usize,
                        )
                    }
                }
            }
        }
    }

    // pub fn new(matrix_type: GateMatrixType, target_qubit: &Vec<usize>, matrix: &Array2<Complex64>, target_commute: &Vec<usize>) -> QuantumGate {
    //     QuantumGate {matrix_type, }

    // }
}

impl fmt::Display for QuantumGate {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut qs = "".to_string();
        // TODO: commute
        // for qubit in &self.target_qubit_index {
        //     //qs.push_str(format!("{} : commute ", qubit));
        //     let c = qubit.to_string() + " : commute ";
        //     qs.push_str(&c);
        //     //if qubit.is_power_of_two()
        //     // << (val.is_commute_X() ? "X" : " ") << " "
        //     // << (val.is_commute_Y() ? "Y" : " ") << " "
        //     // << (val.is_commute_Z() ? "Z" : " ") << " " << std::endl;
        // }
        let msg = format!(
            r" *** gate info *** 
* gate name : {}
* target    : {}
",
            "name", "",
        );
        write!(f, "{}", msg)
    }
}

pub trait DenseMatrixGate {
    fn DenseMatrixGate(
        target_qubit_index: Vec<usize>,
        dense_matrix_element: Array2<Complex64>,
        target_qubit_commutation: Vec<usize>,
        name: &str,
    ) -> QuantumGate {
        // TODO: validation
        //let dim = 1 << target_qubit_index.len();
        // if ((unsigned)matrix.cols() != dim)
        //     throw std::invalid_argument("matrix.cols() != dim");
        // if ((unsigned)matrix.rows() != dim)
        //     throw std::invalid_argument("matrix.rows() != dim");
        QuantumGate {
            matrix_type: GateMatrixType::DenseMatrix,
            name: name.to_string(),
            target_qubit_index,
            target_qubit_commutation,
            control_qubit_index: vec![],
            control_qubit_value: vec![],
            dense_matrix_element,
            pauli_id: vec![],
            rotation_angle: 0.,
        }
    }
}
impl DenseMatrixGate for QuantumGate {}

pub trait PauliMatrixGate {
    fn PauliMatrixGate(
        target_qubit_index: Vec<usize>,
        pauli_id: Vec<usize>,
        rotation_angle: f64,
        name: &str,
    ) -> QuantumGate {
        // TODO: validation
        // if (pauli_id.size() != target_qubit.size())
        // throw std::invalid_argument(
        //     "pauli_id.size() != target_qubit.size()");
        //target_commute
        let mut target_commute = Vec::with_capacity(target_qubit_index.len());
        for i in 0..target_qubit_index.len() {
            target_commute.push(pauli_id[i]);
        }
        QuantumGate {
            matrix_type: GateMatrixType::PauliMatrix,
            name: name.to_string(),
            target_qubit_index,
            target_qubit_commutation: target_commute,
            control_qubit_index: vec![],
            control_qubit_value: vec![],
            dense_matrix_element: Array2::eye(0),
            pauli_id,
            rotation_angle,
        }
    }
}
impl PauliMatrixGate for QuantumGate {}

pub trait Identity {
    fn Identity(target_qubit: usize) -> QuantumGate {
        // let target_qubit_index = vec![target_qubit];
        // let dense_matrix_element: Array2<Complex64> = Array2::eye(2);
        // let target_qubit_commutation = vec![FLAG_COMMUTE_X | FLAG_COMMUTE_Y | FLAG_COMMUTE_Z];
        QuantumGate::DenseMatrixGate(
            vec![target_qubit],
            Array2::eye(2),
            vec![FLAG_COMMUTE_X | FLAG_COMMUTE_Y | FLAG_COMMUTE_Z],
            "Identity",
        )
    }
}
impl Identity for QuantumGate {}

pub trait X {
    fn X(target_qubit: usize) -> QuantumGate {
        QuantumGate::DenseMatrixGate(
            vec![target_qubit],
            array![[COMP_ZERO, COMP_ONE], [COMP_ONE, COMP_ZERO]],
            vec![FLAG_COMMUTE_X],
            "X",
        )
    }
}
impl X for QuantumGate {}

pub trait Y {
    fn Y(target_qubit: usize) -> QuantumGate {
        QuantumGate::DenseMatrixGate(
            vec![target_qubit],
            array![[COMP_ZERO, COMP_MINUS_I], [COMP_I, COMP_ZERO]],
            vec![FLAG_COMMUTE_Y],
            "Y",
        )
    }
}
impl Y for QuantumGate {}

pub trait Z {
    fn Z(target_qubit: usize) -> QuantumGate {
        QuantumGate::DenseMatrixGate(
            vec![target_qubit],
            array![[COMP_ONE, COMP_ZERO], [COMP_ZERO, COMP_MINUS_ONE]],
            vec![FLAG_COMMUTE_Z],
            "Z",
        )
    }
}
impl Z for QuantumGate {}
