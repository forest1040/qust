use std::collections::HashMap;

enum MapType {
    Basic,
    Sequence,
    Probabilistic,
    CPTP,
    Instrument,
}

enum GateMatrixType {
    DenseMatrix,
    DiagonalMatrix,
    SparseMatrix,
    PauliMatrix,
    PermutationMatrix,
}

enum SpecialFuncType {
    None,
    GateI,
    GateX,
    GateY,
    GateZ,
    GateSqrtX,
    GateSqrtY,
    GateSqrtXdag,
    GateSqrtYdag,
    GateRX,
    GateRY,
    GateRZ,
    GateH,
    GateS,
    GateSdag,
    GateT,
    GateTdag,
    GateP0,
    GateP1,
    GateCX,
    GateCY,
    GateCZ,
    GateSWAP,
}

struct QuantumGateBasic {
    map_type: MapType,
    matrix_type: GateMatrixType,
    special_func_type: SpecialFuncType,
    target_qubit_index: Vec<usize>,
    target_qubit_commutation: Vec<usize>,
    control_qubit_index: Vec<usize>,
    control_qubit_value: Vec<usize>,
    gate_property: usize,

    dense_matrix_element: ComplexMatrix,

    parameter: HashMap<String, f64>,
}
