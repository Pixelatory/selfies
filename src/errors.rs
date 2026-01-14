use std::fmt;

#[derive(Debug)]
pub enum GraphConstructionError {
    /// An unexpected token was used during graph construction.
    UnexpectedToken {
        smiles: String,
        message: String,
        index: usize,
    },
    /// An edge is referenced through a pair of nodes, but the edge does not exist.
    UnknownEdge {
        smiles: String,
        message: String,
        src: usize,
        dst: usize,
    },
    CannotKekulize {
        smiles: String
    },
    InvalidBondConstraints {message: String},
    EmptySMILES
}

impl fmt::Display for GraphConstructionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedToken {smiles, message, index} => {
                write!(f, "{message} Index: {index}. SMILES: {smiles}.")
            },
            Self::UnknownEdge { smiles, message, src, dst } => {
                write!(f, "{message} Bond: ({src}, {dst}). SMILES: {smiles}.")
            },
            Self::CannotKekulize { smiles } => {
                write!(f, "Could not kekulize molecule: {smiles}.")
            },
            Self::InvalidBondConstraints { message } => {
                write!(f, "{message}")
            },
            Self::EmptySMILES => write!(f, "Empty SMILES string."),
        }
    }
}

impl std::error::Error for GraphConstructionError {}

#[derive(Debug)]
pub enum EncoderError {
    UnexpectedBondOrder {
        smiles: String,
        message: String,
        src: usize,
        dst: usize,
    },
    AromaticAtom {
        smiles: String,
        message: String,
        atom_index: usize,
    }
}

impl fmt::Display for EncoderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnexpectedBondOrder { smiles, message, src, dst } => {
                write!(f, "{message} Bond: ({src}, {dst}). SMILES: {smiles}.")
            },
            Self::AromaticAtom { smiles, message, atom_index } => {
                write!(f, "{message} Atom index: {atom_index}. SMILES: {smiles}.")
            },
        }
    }
}

impl std::error::Error for EncoderError {}