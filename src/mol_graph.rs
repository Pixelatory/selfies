use std::{collections::HashMap, iter::Map};

use crate::smiles_utils::SMILESToken;


/// A molecular graph.
/// 
/// Molecules can be viewed as weighted undirected graphs. However, SMILES
/// and SELFIES strings are more naturally represented as weighted directed
/// graphs, where the direction of the edges specifies the order of atoms
/// and bonds in the string.
pub struct MolecularGraph<'atom> {
    /// Stores root atoms, where traversal begins.
    roots: Vec<usize>,
    /// Stores atoms in this graph.
    atoms: Vec<&'atom Atom>,
    /// Stores all bonds in this graph.
    bond_dict: Vec<i32>,
    /// Adjacency list, representing this graph.
    adj_list: Vec<i32>,
    /// Stores number of bonds an atom has made.
    bond_counts: Vec<i32>,
    /// Stores if an atom makes a ring bond.
    ring_bond_flags: Vec<i32>,
    /// Delocalization subgraph.
    delocal_subgraph: Vec<i32>,
    /// Attribution of each atom/bond.
    attribution: HashMap<i32, i32>,
    attributable: bool,
}

impl<'atom> MolecularGraph<'atom>{
    pub fn new(attributable: bool) -> Self {
        Self {
            roots: Vec::new(),
            atoms: Vec::new(),
            bond_dict: Vec::new(),
            adj_list: Vec::new(),
            bond_counts: Vec::new(),
            ring_bond_flags: Vec::new(),
            delocal_subgraph: Vec::new(),
            attribution: HashMap::new(),
            attributable,
        }
    }
}

#[derive(Debug, PartialEq)]
pub struct Atom {
    pub index: Option<usize>,
    pub element: String,
    pub is_aromatic: bool,
    pub isotope: Option<i32>,
    pub chirality: Option<String>,
    pub h_count: i32,
    pub charge: i32,
}

impl Atom {
    pub fn new(
        element: String,
        is_aromatic: bool,
        isotope: Option<i32>,
        chirality: Option<String>,
        h_count: i32,
        charge: i32,
    ) -> Self {
        Atom {
            index: None,
            element: element,
            is_aromatic: is_aromatic,
            isotope: isotope,
            chirality: chirality,
            h_count: h_count,
            charge: charge,
        }
    }
}