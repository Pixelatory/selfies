use std::{collections::HashMap, iter::Map};


/// A molecular graph.
/// 
/// Molecules can be viewed as weighted undirected graphs. However, SMILES
/// and SELFIES strings are more naturally represented as weighted directed
/// graphs, where the direction of the edges specifies the order of atoms
/// and bonds in the string.
pub struct MolecularGraph {
    /// Stores root atoms, where traversal begins.
    roots: Vec<i32>,
    /// Stores atoms in this graph.
    atoms: Vec<i32>,
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

impl MolecularGraph{
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
    pub element: String,
    pub is_aromatic: bool,
    pub isotope: Option<i32>,
    pub chirality: Option<String>,
    pub h_count: i32,
    pub charge: i32,
}