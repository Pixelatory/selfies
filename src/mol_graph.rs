use std::{collections::{HashMap, HashSet, VecDeque}, fmt};

use crate::{
    bond_constraints::{check_bond_constraints},
    constants::{AROMATIC_VALENCES, SMILES_BOND_ORDERS, SMILES_STEREO_BONDS, VALENCE_ELECTRONS},
    matching_utils::find_perfect_matching,
    smiles_utils::{SMILESToken, SMILESTokenType, SMILESTokenizer, smiles_to_atom},
    utilities::{get_mut_pair, last_valence, valence_any}
};

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

struct MolecularGraphContext {
    mol_graph: MolecularGraph,
    attribution: Option<Attribution>,
    delocalized_subgraph: Subgraph,
}

/// The delocalized subgraph.
struct Subgraph {
    adj_list: HashMap<usize, HashSet<usize>>
}

impl Subgraph {
    pub fn new() -> Self {
        Self {
            adj_list: HashMap::new()
        }
    }
}

struct AttributeReason {
    // TODO: fill me out!
}

struct Attribution {
    // e.g. attribute reason for each edge
    edge_attributes: HashMap<(usize, usize), AttributeReason>
}

impl Attribution {
    pub fn new() -> Self {
        Self {
            edge_attributes: HashMap::new()
        }
    }
}

/// A molecular graph.
/// 
/// Molecules can be viewed as weighted undirected graphs. However, SMILES
/// and SELFIES strings are more naturally represented as weighted directed
/// graphs, where the direction of the edges specifies the order of atoms
/// and bonds in the string.
#[derive(Debug)]
pub struct MolecularGraph {
    /// Stores atom data for this graph.
    atom_data: Vec<AtomData>,
    /// Stores all bonds in this graph.
    bond_map: HashMap<(usize, usize), DirectedBond>,
    /// Adjacency list, representing this graph. Stores indices of atoms.
    adj_list: Vec<HashSet<usize>>,
}

impl MolecularGraph{
    fn new() -> Self {
        Self {
            atom_data: Vec::new(),
            bond_map: HashMap::new(),
            adj_list: Vec::new(),
        }
    }

    fn has_bond(&self, source: usize, destination: usize) -> bool {
        self.bond_map.contains_key(&(source, destination))
    }

    /// Adds an atom to the molecular graph.
    /// 
    /// Returns the index of the inserted atom.
    fn add_atom(&mut self, atom: &Atom) -> usize {
        self.atom_data.push(AtomData::new(atom.clone()));
        self.adj_list.push(HashSet::new());

        self.adj_list.len() - 1
    }

    fn add_bond(&mut self, source: usize, destination: usize, order: f64, stereo_bond_char: Option<char>) {
        let src_dst_bond = DirectedBond::new(source, destination, order, stereo_bond_char, false);
        let dst_src_bond = DirectedBond::new(destination, source, order, stereo_bond_char, false);

        self.bond_map.insert((source, destination), src_dst_bond);
        self.bond_map.insert((destination, source), dst_src_bond);

        self.adj_list[source].insert(destination);
        self.adj_list[destination].insert(source);

        self.atom_data[source].bond_count += order;
        self.atom_data[destination].bond_count += order;
    }

    fn add_ring_bonds(&mut self, l_atom_idx: usize, r_atom_idx: usize, order: f64, l_atom_stereo: Option<char>, r_atom_stereo: Option<char>) {
        let l_bond = DirectedBond::new(l_atom_idx, r_atom_idx, order, l_atom_stereo, true);
        let r_bond = DirectedBond::new(r_atom_idx, l_atom_idx, order, r_atom_stereo, true);

        self.bond_map.insert((l_atom_idx, r_atom_idx), l_bond);
        self.bond_map.insert((r_atom_idx, l_atom_idx), r_bond);

        self.adj_list[l_atom_idx].insert(r_atom_idx);
        self.adj_list[r_atom_idx].insert(l_atom_idx);

        self.atom_data[l_atom_idx].bond_count += order;
        self.atom_data[r_atom_idx].bond_count += order;
    }

    pub fn get_atom_data(&self) -> &Vec<AtomData> {
        return &self.atom_data;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DirectedBond {
    src: usize,
    dst: usize,
    order: f64,
    stereo: Option<char>,
    ring_bond: bool,
}

impl DirectedBond {
    pub fn new(
        src: usize,
        dst: usize,
        order: f64,
        stereo: Option<char>,
        ring_bond: bool,
    ) -> Self {
        DirectedBond {src, dst, order, stereo, ring_bond}
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct Atom {
    pub element: String,
    pub is_aromatic: bool,
    pub isotope: Option<u32>,
    pub chirality: Option<String>,
    pub h_count: Option<u32>,
    pub charge: i32,
}

impl Atom {
    pub fn new(
        element: String,
        is_aromatic: bool,
        isotope: Option<u32>,
        chirality: Option<String>,
        h_count: Option<u32>,
        charge: i32,
    ) -> Self {
        Atom {
            element: element,
            is_aromatic: is_aromatic,
            isotope: isotope,
            chirality: chirality,
            h_count: h_count,
            charge: charge,
        }
    }
}

#[derive(Debug)]
pub struct AtomData {
    pub atom: Atom,
    pub bond_count: f64,
}

impl AtomData {
    pub fn new(atom: Atom) -> Self {
        Self {
            atom: atom,
            bond_count: 0.0,
        }
    }
}

/// Reads a molecular graph from a SMILES string.
pub fn create_mol_graph(smiles: &str, kekulize: bool, strict: bool) -> Result<MolecularGraph, GraphConstructionError> {
    if smiles.is_empty() {
        return Err(GraphConstructionError::EmptySMILES);
    }
    let mut mol = MolecularGraph::new();
    let tokens = SMILESTokenizer::new(smiles).into_iter().collect::<Result<Vec<SMILESToken>, GraphConstructionError>>()?;
    let mut tokens_queue = VecDeque::from(tokens);

    derive_mol_from_tokens(&mut mol, &mut tokens_queue, smiles)?;

    if kekulize {
        let subgraph = create_delocalized_subgraph(&mol);
        let kekulized = kekulize_mol(smiles, &mut mol, subgraph)?;
        if !kekulized {
            return Err(GraphConstructionError::CannotKekulize { smiles: smiles.to_string() });
        }
    }

    if strict {
        check_bond_constraints(&mol, smiles)?;
    }

    return Ok(mol)
}

/// Dearomatizes atoms in the molecular graph according to the subgraph.
/// 
/// Atoms that have a corresponding sub-node have their aromatic flag set to false.
/// Bonds that have a corresponding sub-edge have their order set to 1.0.
fn dearomatize_atoms(smiles: &str, mol: &mut MolecularGraph, subgraph: &Subgraph) -> Result<(), GraphConstructionError> {
    for (src, dst_vec) in subgraph.adj_list.iter() {
        for dst in dst_vec {
            let src = *src;
            let dst = *dst;

            let (src_atom, dst_atom) = get_mut_pair(&mut mol.atom_data, src, dst);
            
            // Reset the aromatic bond information. Use single bonds for the cycle.
            let [maybe_src_dst_bond, maybe_dst_src_bond] = mol.bond_map.get_disjoint_mut([&(src, dst), &(dst, src)]);

            match (maybe_src_dst_bond, maybe_dst_src_bond) {
                (Some(src_dst_bond), Some(dst_src_bond)) => {
                    let old_order = src_dst_bond.order;
                    src_atom.bond_count += 1.0 - old_order;
                    dst_atom.bond_count += 1.0 - old_order;

                    src_dst_bond.order = 1.0;
                    dst_src_bond.order = 1.0;
                },
                _ => {
                    return Err(GraphConstructionError::UnknownEdge {
                        smiles: smiles.to_string(),
                        message: format!("An edge in the subgraph does not have a bond in molecular graph.").to_string(),
                        src: src,
                        dst: dst,
                    });
                }
            }
            src_atom.atom.is_aromatic = false;
            dst_atom.atom.is_aromatic = false;
        }
    }

    Ok(())
}

/// Kekulize a molecule.
/// 
/// This removes a molecule's aromatic data. Aromatic atoms/bonds will be set
/// with either a single or double bond.
fn kekulize_mol(smiles: &str, mol: &mut MolecularGraph, subgraph: Subgraph) -> Result<bool, GraphConstructionError> {
    if let Some(matching) = find_perfect_matching(&subgraph.adj_list) {
        dearomatize_atoms(smiles, mol, &subgraph)?;

        // Edges in the matching have their order set to 2.0 for a double bond.
        for (src, dst) in matching.iter() {
            let maybe_src_dst_bond = mol.bond_map.get_mut(&(*src, *dst));
            if let Some(src_dst_bond) = maybe_src_dst_bond {
                let (src_atom, dst_atom) = get_mut_pair(&mut mol.atom_data, *src, *dst);
                src_atom.bond_count += 2.0 - src_dst_bond.order;
                dst_atom.bond_count += 2.0 - src_dst_bond.order;
                src_dst_bond.order = 2.0;
            } else {
                return Err(GraphConstructionError::UnknownEdge {
                    smiles: smiles.to_string(),
                    message: format!("An edge in the perfect matching does not have an edge in the subgraph.").to_string(),
                    src: *src,
                    dst: *dst,
                });
            }
        }
        return Ok(true);
    }

    Ok(false)
}

/// Creates a delocalized subgraph from a constructed MolecularGraph.
/// 
/// The result is an already pruned subgraph.
fn create_delocalized_subgraph(mol: &MolecularGraph) -> Subgraph {
    let mut subgraph = Subgraph::new();
    let mut should_prune_map: HashMap<usize, bool> = HashMap::new();

    for bond in mol.bond_map.values() {
        if bond.order == 1.5 {
            let should_prune_src = memoize_should_prune(&mut should_prune_map, mol, bond.src);
            let should_prune_dst = memoize_should_prune(&mut should_prune_map, mol, bond.dst);
            if should_prune_src || should_prune_dst {
                continue;
            }
            subgraph.adj_list.entry(bond.src).or_default().insert(bond.dst);
            subgraph.adj_list.entry(bond.dst).or_default().insert(bond.src);
        }
    }
    subgraph
}

/// Memoizes the output of should_prune().
/// 
/// Memoized outputs are stured in `should_prune_map`.
fn memoize_should_prune(should_prune_map: &mut HashMap<usize, bool>, mol: &MolecularGraph, node: usize) -> bool {
    match should_prune_map.get(&node) {
        Some(result) => return *result,
        None => {
            let result = should_prune(mol, node);
            should_prune_map.insert(node, result);
            return result;
        }
    }
}

/// Whether the node (atom index in MolecularGraph) should be pruned from the delocalized subgraph.
fn should_prune(mol: &MolecularGraph, node: usize) -> bool {
    let adj_nodes = &mol.adj_list[node];
    if adj_nodes.is_empty() {
        return true;
    }

    let atom = &mol.atom_data[node].atom;
    let valences = &AROMATIC_VALENCES[&atom.element.as_str()];

    // Each bond in delocalized subgraph has order 1.5 - treat them as single bonds.
    let mut used_electrons = adj_nodes.len() as i32;

    // Account for implicit Hs.
    return match atom.h_count {
        None => {
            assert!(atom.charge == 0);
            valence_any(&valences, used_electrons)
        },
        Some(h_count) => {
            let bond_count = adj_nodes.iter().map(|x| mol.bond_map[&(node, *x)].order).sum::<f64>() as i32;
            let valence = last_valence(&valences) - atom.charge;
            used_electrons += h_count as i32;

            let bound_electrons = atom.charge.max(0) + h_count as i32 + bond_count;
            let radical_electrons = (VALENCE_ELECTRONS[&atom.element.as_str()] - bound_electrons).max(0) % 2;
            let free_electrons = valence - used_electrons - radical_electrons;

            if valence_any(valences, used_electrons + atom.charge) {
                return true;
            } else {
                return !((free_electrons >= 0) && (free_electrons % 2 != 0));
            }
        }
    };
}

fn handle_atom_token(
    smiles: &str,
    token: &SMILESToken,
    mol: &mut MolecularGraph,
    prev_stack: &mut Vec<usize>,
) -> Result<(), GraphConstructionError>{
    let curr = match smiles_to_atom(&token.token) {
        Some(atom) => atom,
        None => {
            return Err(
                GraphConstructionError::UnexpectedToken {
                    smiles: smiles.to_string(),
                    message: format!("Invalid atom symbol {}", token.token),
                    index: token.start_idx,
                }
            )
        },
    };
    
    // Add the atom to the graph.
    let inserted_idx = mol.add_atom(&curr);
    
    // Add bonds to the graph if there's a previous atom.
    if let Some(prev_atom_idx) = prev_stack.last() {
        let (order, stereo_bond_char) = smiles_to_bond(
            smiles,
            token.start_idx,
            token.bond_token,
            &curr,
            &mol.atom_data[*prev_atom_idx].atom,
        )?;
        mol.add_bond(*prev_atom_idx, inserted_idx, order, stereo_bond_char);
    }

    // Remove the previous atom index, and add the current inserted index as the previous index.
    prev_stack.pop();
    prev_stack.push(inserted_idx);

    Ok(())
}

fn handle_branch_token(
    smiles: &str,
    token: &SMILESToken,
    prev_stack: &mut Vec<usize>,
    branch_stack: &mut Vec<(String, usize)>,
) -> Result<(), GraphConstructionError> {
    if &token.token == "(" {
        if let Some(prev_atom_idx) = prev_stack.last() {
            prev_stack.push(*prev_atom_idx);
            branch_stack.push((token.token.clone(), token.start_idx));
        } else {
            return Err(GraphConstructionError::UnexpectedToken {
                smiles: smiles.to_string(),
                message: "Branch has no previous atom.".to_string(),
                index: token.start_idx,
            });
        }
    } else {
        if branch_stack.is_empty() {
            return Err(GraphConstructionError::UnexpectedToken {
                smiles: smiles.to_string(),
                message: "Hanging ')' bracket.".to_string(),
                index: token.start_idx,
            });
        }
        branch_stack.pop();
        prev_stack.pop();
    }

    Ok(())
}

fn handle_ring_token(
    smiles: &str,
    token: &SMILESToken,
    mol: &mut MolecularGraph,
    prev_stack: &Vec<usize>,
    ring_log: &mut HashMap<String, (Option<char>, usize)>,
) -> Result<(), GraphConstructionError> {
    if let Some((maybe_latom_bond_char, latom_idx)) = ring_log.remove(&token.token) {
        // The ending ring bond token.
        if let Some(ratom_idx) = prev_stack.last() {
            let ratom_idx = *ratom_idx;
            // Validate whether a ring bond should be added.
            if mol.has_bond(latom_idx, ratom_idx) {
                return Err(GraphConstructionError::UnexpectedToken {
                    smiles: smiles.to_string(),
                    message: "Attempted to make a ring bond between already-bonded atoms.".to_string(),
                    index: token.start_idx,
                });
            }
            
            match (maybe_latom_bond_char, token.bond_token) {
                (Some(latom_bond_char), Some(ratom_bond_char)) => {
                    if latom_bond_char != ratom_bond_char && 
                    (!SMILES_STEREO_BONDS.contains(&latom_bond_char) || !SMILES_STEREO_BONDS.contains(&ratom_bond_char)) {
                        return Err(GraphConstructionError::UnexpectedToken {
                            smiles: smiles.to_string(),
                            message: "A ring bond is specified at both ends, but the bond type does not match.".to_string(),
                            index: token.start_idx,
                        });
                    }
                }
                _ => {},
            };
            
            // Attempt to include a ring bond.
            let latom = &mol.atom_data[latom_idx].atom;
            let ratom = &mol.atom_data[ratom_idx].atom;

            let (l_order, l_stereo) = smiles_to_bond(smiles, token.start_idx, maybe_latom_bond_char, ratom, latom)?;
            let (r_order, r_stereo) = smiles_to_bond(smiles, token.start_idx, token.bond_token, latom, ratom)?;
            
            let order: f64;
            if latom.is_aromatic && ratom.is_aromatic && maybe_latom_bond_char.is_none() && token.bond_token.is_none() {
                order = 1.5;
            } else {
                order = l_order.max(r_order);
            }

            mol.add_ring_bonds(latom_idx, ratom_idx, order, l_stereo, r_stereo);
        } else {
            return Err(GraphConstructionError::UnexpectedToken {
                smiles: smiles.to_string(),
                message: "Ending ring token has no previous atom.".to_string(),
                index: token.start_idx,
            });
        }
    } else {
        // The starting ring bond token.
        if let Some(prev_atom_idx) = prev_stack.last() {
            ring_log.insert(token.token.clone(), (token.bond_token, *prev_atom_idx));
        } else {
            return Err(GraphConstructionError::UnexpectedToken {
                smiles: smiles.to_string(),
                message: "Starting ring token has no previous atom.".to_string(),
                index: token.start_idx,
            });
        }
    }
    Ok(())
}

fn derive_mol_from_tokens(
    mol: &mut MolecularGraph,
    tokens: &mut VecDeque<SMILESToken>,
    smiles: &str,
) -> Result<(), GraphConstructionError> {
    // Holds indexes of previous atoms in the graph.
    let mut prev_stack = Vec::new();
    // Holds (branch token char, branch token index).
    let mut branch_stack = Vec::new();
    // Maps a ring token string to tuple of bond char and previous atom index.
    let mut ring_log = HashMap::new();

    while !tokens.is_empty() {
        let token = tokens.pop_front().unwrap();
        match token.token_type {
            SMILESTokenType::Dot => {
                // Clear all containers, as if starting from scratch with a new molecule.
                prev_stack.clear();
                branch_stack.clear();
                ring_log.clear();
            },
            SMILESTokenType::Atom => handle_atom_token(smiles, &token, mol, &mut prev_stack)?,
            SMILESTokenType::Branch => handle_branch_token(smiles, &token, &mut prev_stack, &mut branch_stack)?,
            SMILESTokenType::Ring => handle_ring_token(smiles, &token, mol, &prev_stack, &mut ring_log)?,
        }
    }
    Ok(())
}

fn smiles_to_bond(smiles: &str, token_idx: usize, maybe_bond_char: Option<char>, curr_atom: &Atom, prev_atom: &Atom) -> Result<(f64, Option<char>), GraphConstructionError> {
    match maybe_bond_char {
        Some(bond_char) => {
            if let Some(bond_order) = SMILES_BOND_ORDERS.get(&bond_char) {
                let stereo_bond_char = SMILES_STEREO_BONDS.get(&bond_char).copied();
                return Ok((*bond_order, stereo_bond_char));
            } else {
                // All bond types must be defined in SMILES_BOND_ORDERS.
                return Err(GraphConstructionError::UnexpectedToken {
                    smiles: smiles.to_string(),
                    message: format!("Unknown bond character: '{}'.", bond_char),
                    index: token_idx,
                });
            }
        }
        None => {
            if curr_atom.is_aromatic && prev_atom.is_aromatic {
                // The same bond order as a directly specified aromatic bond (':').
                return Ok((1.5, None));
            }
            // Assumes single bond ('-') if not specified.
            return Ok((1.0, None));
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_string() {
        let x = create_mol_graph("", false, false);
        assert!(matches!(x, Err(GraphConstructionError::EmptySMILES)));
    }

    #[test]
    fn test_create_mol_graph_stereo() {
        let x = create_mol_graph("F/C=C\\F", false, false);
        assert!(x.is_ok());
        let x = x.unwrap();

        let expected_atoms = vec![
            Atom::new("F".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("F".to_string(), false, None, None, None, 0),
        ];
        let expected_adj_list = vec![
            HashSet::from([1]),
            HashSet::from([0, 2]),
            HashSet::from([1, 3]),
            HashSet::from([2]),
        ];
        let expected_bond_map = HashMap::from(
            [
                ((0, 1), DirectedBond::new(0, 1, 1.0, Some('/'), false)),
                ((1, 0), DirectedBond::new(1, 0, 1.0, Some('/'), false)),
                ((1, 2), DirectedBond::new(1, 2, 2.0, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 2.0, None, false)),
                ((2, 3), DirectedBond::new(2, 3, 1.0, Some('\\'), false)),
                ((3, 2), DirectedBond::new(3, 2, 1.0, Some('\\'), false)),
            ]
        );
        let expected_bond_count = vec![1.0, 3.0, 3.0, 1.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert_eq!(x.bond_map, expected_bond_map);
        assert_eq!(bond_count, expected_bond_count);
    }

    #[test]
    fn test_create_mol_graph_dot() {
        let x = create_mol_graph("[Cu+2].[O-]S(=O)(=O)[O-]", false, false);
        assert!(x.is_ok());
        let x = x.unwrap();

        let expected_atoms = vec![
            Atom::new("Cu".to_string(), false, None, None, Some(0), 2),
            Atom::new("O".to_string(), false, None, None, Some(0), -1),
            Atom::new("S".to_string(), false, None, None, None, 0),
            Atom::new("O".to_string(), false, None, None, None, 0),
            Atom::new("O".to_string(), false, None, None, None, 0),
            Atom::new("O".to_string(), false, None, None, Some(0), -1),
        ];
        let expected_adj_list = vec![
            HashSet::new(),  // Cu+2 has no bonds with any other atom.
            HashSet::from([2]),
            HashSet::from([1, 3, 4, 5]),
            HashSet::from([2]),
            HashSet::from([2]),
            HashSet::from([2]),
        ];
        let expected_bond_map: HashMap<(usize, usize), DirectedBond> = HashMap::from(
            [
                ((1, 2), DirectedBond::new(1, 2, 1.0, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 1.0, None, false)),
                ((2, 3), DirectedBond::new(2, 3, 2.0, None, false)),
                ((3, 2), DirectedBond::new(3, 2, 2.0, None, false)),
                ((2, 4), DirectedBond::new(2, 4, 2.0, None, false)),
                ((4, 2), DirectedBond::new(4, 2, 2.0, None, false)),
                ((2, 5), DirectedBond::new(2, 5, 1.0, None, false)),
                ((5, 2), DirectedBond::new(5, 2, 1.0, None, false)),
            ]
        );
        let expected_bond_count = vec![0.0, 1.0, 6.0, 2.0, 2.0, 1.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert_eq!(x.bond_map, expected_bond_map);
        assert_eq!(bond_count, expected_bond_count);
    }

    #[test]
    fn test_create_mol_graph_ring_bond() {
        let x = create_mol_graph("O1C(CCl)=CCN=1", false, false);
        assert!(x.is_ok());
        let x = x.unwrap();

        let expected_atoms = vec![
            Atom::new("O".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("Cl".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("N".to_string(), false, None, None, None, 0),
        ];
        let expected_adj_list = vec![
            HashSet::from([1, 6]),
            HashSet::from([0, 2, 4]),
            HashSet::from([1, 3]),
            HashSet::from([2]),
            HashSet::from([1, 5]),
            HashSet::from([4, 6]),
            HashSet::from([5, 0]),
        ];
        let expected_bond_map: HashMap<(usize, usize), DirectedBond> = HashMap::from(
            [
                ((0, 1), DirectedBond::new(0, 1, 1.0, None, false)),
                ((1, 0), DirectedBond::new(1, 0, 1.0, None, false)),

                ((1, 2), DirectedBond::new(1, 2, 1.0, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 1.0, None, false)),

                ((2, 3), DirectedBond::new(2, 3, 1.0, None, false)),
                ((3, 2), DirectedBond::new(3, 2, 1.0, None, false)),

                ((4, 1), DirectedBond::new(4, 1, 2.0, None, false)),
                ((1, 4), DirectedBond::new(1, 4, 2.0, None, false)),

                ((4, 5), DirectedBond::new(4, 5, 1.0, None, false)),
                ((5, 4), DirectedBond::new(5, 4, 1.0, None, false)),

                ((5, 6), DirectedBond::new(5, 6, 1.0, None, false)),
                ((6, 5), DirectedBond::new(6, 5, 1.0, None, false)),

                ((6, 0), DirectedBond::new(6, 0, 2.0, None, true)),
                ((0, 6), DirectedBond::new(0, 6, 2.0, None, true)),
            ]
        );
        let expected_bond_count = vec![3.0, 4.0, 2.0, 1.0, 3.0, 2.0, 3.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert_eq!(x.bond_map, expected_bond_map);
        assert_eq!(bond_count, expected_bond_count);
    }

    /// A molecule cannot be kekulized if there is no perfect matching.
    #[test]
    fn test_unkekulized_no_perfect_matching() {
        let x = create_mol_graph("n1c[nH]cc1", false, false);
        assert!(x.is_ok());
        let x = x.unwrap();
        // [nH] is pruned from the subgraph.
        assert_eq!(should_prune(&x, 2), true);

        let expected_atoms = vec![
            Atom::new("N".to_string(), true, None, None, None, 0),
            Atom::new("C".to_string(), true, None, None, None, 0),
            Atom::new("N".to_string(), true, None, None, Some(1), 0),
            Atom::new("C".to_string(), true, None, None, None, 0),
            Atom::new("C".to_string(), true, None, None, None, 0),
        ];
        let expected_adj_list = vec![
            HashSet::from([1, 4]),
            HashSet::from([0, 2]),
            HashSet::from([1, 3]),
            HashSet::from([2, 4]),
            HashSet::from([3, 0]),
        ];
        let expected_bond_map: HashMap<(usize, usize), DirectedBond> = HashMap::from(
            [
                ((0, 1), DirectedBond::new(0, 1, 1.5, None, false)),
                ((1, 0), DirectedBond::new(1, 0, 1.5, None, false)),

                ((1, 2), DirectedBond::new(1, 2, 1.5, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 1.5, None, false)),

                ((2, 3), DirectedBond::new(2, 3, 1.5, None, false)),
                ((3, 2), DirectedBond::new(3, 2, 1.5, None, false)),

                ((4, 3), DirectedBond::new(4, 3, 1.5, None, false)),
                ((3, 4), DirectedBond::new(3, 4, 1.5, None, false)),

                ((0, 4), DirectedBond::new(0, 4, 1.5, None, true)),
                ((4, 0), DirectedBond::new(4, 0, 1.5, None, true)),
            ]
        );
        let expected_bond_count = vec![3.0, 3.0, 3.0, 3.0, 3.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert_eq!(x.bond_map, expected_bond_map);
        assert_eq!(bond_count, expected_bond_count);

        // A CannotKekulize error occurs if we attempt to kekulize and fail.
        let x = create_mol_graph("n1c[nH]cc1", true, false);
        assert!(matches!(x.err().unwrap(), GraphConstructionError::CannotKekulize { smiles: _ }))
    }

    #[test]
    fn test_kekulized_mol_defined_bond() {
        // A single bond is defined in what would otherwise be an aromatic ring.
        let x = create_mol_graph("c1ccc#cc1", true, false);
        assert!(x.is_ok());
        let x = x.unwrap();
        let expected_atoms = vec![
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
        ];
        let expected_adj_list = vec![
            HashSet::from([1, 5]),
            HashSet::from([0, 2]),
            HashSet::from([1, 3]),
            HashSet::from([2, 4]),
            HashSet::from([3, 5]),
            HashSet::from([4, 0]),
        ];
        let expected_bond_map: HashMap<(usize, usize), DirectedBond> = HashMap::from(
            [
                ((0, 1), DirectedBond::new(0, 1, 2.0, None, false)),
                ((1, 0), DirectedBond::new(1, 0, 2.0, None, false)),

                ((1, 2), DirectedBond::new(1, 2, 1.0, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 1.0, None, false)),

                ((2, 3), DirectedBond::new(2, 3, 2.0, None, false)),
                ((3, 2), DirectedBond::new(3, 2, 2.0, None, false)),

                ((4, 3), DirectedBond::new(4, 3, 3.0, None, false)),
                ((3, 4), DirectedBond::new(3, 4, 3.0, None, false)),

                ((4, 5), DirectedBond::new(4, 5, 2.0, None, false)),
                ((5, 4), DirectedBond::new(5, 4, 2.0, None, false)),

                ((5, 0), DirectedBond::new(5, 0, 1.0, None, true)),
                ((0, 5), DirectedBond::new(0, 5, 1.0, None, true)),
            ]
        );
        let expected_bond_count = vec![4.0, 4.0, 4.0, 6.0, 6.0, 4.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert_eq!(x.bond_map, expected_bond_map);
        assert_eq!(bond_count, expected_bond_count);
    }

    /// A kekulized molecule removes the aromatic flags.
    #[test]
    fn test_kekulized_mol() {
        let x = create_mol_graph("c1cc(ccc1)C", true, false);
        assert!(x.is_ok());
        let x = x.unwrap();
        let expected_atoms = vec![
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
            Atom::new("C".to_string(), false, None, None, None, 0),
        ];
        let expected_adj_list = vec![
            HashSet::from([1, 5]),
            HashSet::from([0, 2]),
            HashSet::from([1, 3, 6]),
            HashSet::from([2, 4]),
            HashSet::from([3, 5]),
            HashSet::from([4, 0]),
            HashSet::from([2]),
        ];
        let expected_bond_map_1: HashMap<(usize, usize), DirectedBond> = HashMap::from(
            [
                ((0, 1), DirectedBond::new(0, 1, 1.0, None, false)),
                ((1, 0), DirectedBond::new(1, 0, 1.0, None, false)),

                ((1, 2), DirectedBond::new(1, 2, 2.0, None, false)),
                ((2, 1), DirectedBond::new(2, 1, 2.0, None, false)),

                ((2, 3), DirectedBond::new(2, 3, 1.0, None, false)),
                ((3, 2), DirectedBond::new(3, 2, 1.0, None, false)),

                ((4, 3), DirectedBond::new(4, 3, 2.0, None, false)),
                ((3, 4), DirectedBond::new(3, 4, 2.0, None, false)),

                ((4, 5), DirectedBond::new(4, 5, 1.0, None, false)),
                ((5, 4), DirectedBond::new(5, 4, 1.0, None, false)),

                ((5, 0), DirectedBond::new(5, 0, 2.0, None, true)),
                ((0, 5), DirectedBond::new(0, 5, 2.0, None, true)),

                ((2, 6), DirectedBond::new(2, 6, 1.0, None, false)),
                ((6, 2), DirectedBond::new(6, 2, 1.0, None, false)),
            ]
        );
        // Kekulization could return one of two bond formations, depending on the matching.
        // expected_bond_map_2 reverses the single and double bonds in the aromatic ring.
        let mut expected_bond_map_2 = expected_bond_map_1.clone();
        for ((src, dst), bond) in expected_bond_map_2.iter_mut() {
            if *src != 6 && *dst != 6 {
                if bond.order == 2.0 {
                    bond.order = 1.0;
                } else {
                    bond.order = 2.0;
                }
            }
        }
        let expected_bond_count = vec![4.0, 4.0, 5.0, 4.0, 4.0, 4.0, 1.0];
        let atoms: Vec<Atom> = x.atom_data.iter().map(|x| x.atom.clone()).collect();
        let bond_count: Vec<f64> = x.atom_data.iter().map(|x| x.bond_count).collect();
        assert_eq!(atoms, expected_atoms);
        assert_eq!(x.adj_list, expected_adj_list);
        assert!(x.bond_map == expected_bond_map_1 || x.bond_map == expected_bond_map_2);
        assert_eq!(bond_count, expected_bond_count);
    }

    #[test]
    fn test_strict_check() {
        let x = create_mol_graph("c1cc(ccc1)C", true, true);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), GraphConstructionError::InvalidBondConstraints {message: _}));

        let x = create_mol_graph("c1ccc#cc1", true, true);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), GraphConstructionError::InvalidBondConstraints {message: _}));

        let x = create_mol_graph("F#C", true, true);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), GraphConstructionError::InvalidBondConstraints {message: _}));
    }
}