use std::collections::{HashMap, HashSet, VecDeque};

use crate::{constants::{AROMATIC_VALENCES, VALENCE_ELECTRONS, ValenceTuple}, smiles_utils::{SMILESParserError, SMILESToken, SMILESTokenizer}, utilities::last_valence};


/// A molecular graph.
/// 
/// Molecules can be viewed as weighted undirected graphs. However, SMILES
/// and SELFIES strings are more naturally represented as weighted directed
/// graphs, where the direction of the edges specifies the order of atoms
/// and bonds in the string.
pub struct MolecularGraph {
    /// Stores root atoms, where traversal begins.
    roots: Vec<usize>,
    /// Stores atoms in this graph.
    atoms: Vec<Atom>,
    /// Stores all bonds in this graph.
    bond_dict: HashMap<(usize, usize), DirectedBond>,
    /// Adjacency list, representing this graph.
    adj_list: Vec<Vec<Option<DirectedBond>>>,
    /// Stores number of bonds an atom has made.
    bond_counts: Vec<f64>,
    /// Stores if an atom makes a ring bond.
    ring_bond_flags: Vec<bool>,
    /// Delocalization subgraph.
    /// A mapping from 
    delocal_subgraph: HashMap<usize, Vec<usize>>,
    /// Attribution of each atom/bond.
    attribution: HashMap<i32, i32>,
    attributable: bool,
}

impl MolecularGraph{
    fn new(attributable: bool) -> Self {
        Self {
            roots: Vec::new(),
            atoms: Vec::new(),
            bond_dict: HashMap::new(),
            adj_list: Vec::new(),
            bond_counts: Vec::new(),
            ring_bond_flags: Vec::new(),
            delocal_subgraph: HashMap::new(),
            attribution: HashMap::new(),
            attributable,
        }
    }

    /// Adds an atom to the molecular graph.
    /// 
    /// Returns the index of the inserted atom.
    fn add_atom(mut self, atom: Atom, mark_root: bool) -> usize {
        let num_atoms = self.atoms.len();
        let mut atom = atom;
        atom.index = Some(num_atoms);
        
        if mark_root {
            self.roots.push(num_atoms);
        }

        self.atoms.push(atom);
        self.adj_list.push(Vec::new());
        self.bond_counts.push(0.);
        self.ring_bond_flags.push(false);

        if self.atoms[num_atoms].is_aromatic {
            self.delocal_subgraph.insert(num_atoms, Vec::new());
        }

        num_atoms
    }

    fn add_bond(mut self, src: usize, dst: usize, order: f64, stereo: String) {
        if src >= dst {
            panic!("Source must be less than destination.")
        }

        let bond = DirectedBond::new(src, dst, order, Some(stereo), false);
        self.add_bond_at_loc(bond, None);
        self.bond_counts[src] += order;
        self.bond_counts[dst] += order;

        if order == 1.5 {
            // Insert vectors if the keys do not have entries.
            self.delocal_subgraph.entry(src).or_insert(Vec::new());
            self.delocal_subgraph.entry(dst).or_insert(Vec::new());
            self.delocal_subgraph.get_mut(&src).unwrap().push(dst);
            self.delocal_subgraph.get_mut(&dst).unwrap().push(src);
        }
    }

    fn add_placeholder_bond(mut self, src: usize) -> usize {
        let out_edges = &mut self.adj_list[src];
        out_edges.push(None);
        out_edges.len() - 1
    }

    fn add_ring_bond(
        mut self,
        a: usize,
        b: usize,
        order: f64,
        a_stereo: Option<String>,
        b_stereo: Option<String>,
        maybe_a_pos: Option<usize>,
        maybe_b_pos: Option<usize>,
    ) {
        let a_bond = DirectedBond::new(a, b, order, a_stereo, true);
        let b_bond = DirectedBond::new(a, b, order, b_stereo, true);
        self.add_bond_at_loc(a_bond, maybe_a_pos);
        self.add_bond_at_loc(b_bond, maybe_b_pos);
        self.bond_counts[a] += order;
        self.bond_counts[b] += order;
        self.ring_bond_flags[a] = true;
        self.ring_bond_flags[b] = true;

        if order == 1.5 {
            self.delocal_subgraph.entry(a).or_insert(Vec::new());
            self.delocal_subgraph.entry(b).or_insert(Vec::new());
            self.delocal_subgraph.get_mut(&a).unwrap().push(b);
            self.delocal_subgraph.get_mut(&b).unwrap().push(a);
        }
    }

    fn update_bond_order(mut self, a: usize, b: usize, new_order: f64) {
        if new_order < 1.0 || new_order > 3.0 {
            panic!("new_order must be within [1.0, 3.0]");
        }
        
        let mut a = a;
        let mut b = b;
        if a > b {
            let tmp = a;
            a = b;
            b = tmp;
        }
        
        let a_to_b = self.bond_dict.get(&(a, b)).unwrap();

        if new_order == a_to_b.order {
            return
        }
        
        if a_to_b.ring_bond {
            let b_to_a = self.bond_dict.get_mut(&(b, a)).unwrap();
            b_to_a.order = new_order;
        }

        let a_to_b = self.bond_dict.get_mut(&(a, b)).unwrap();

        let old_order = a_to_b.order;
        a_to_b.order = new_order;
        self.bond_counts[a] += new_order - old_order;
        self.bond_counts[b] += new_order - old_order;
    }

    fn add_bond_at_loc(&mut self, bond: DirectedBond, maybe_pos: Option<usize>) {
        let bond_src = bond.src;
        self.bond_dict.insert((bond.src, bond.dst), bond.clone());

        let out_edges = &mut self.adj_list[bond_src];
        if (maybe_pos.is_none()) || (maybe_pos.unwrap() == out_edges.len()) {
            out_edges.push(Some(bond));
        } else if let Some(pos) = maybe_pos {
            out_edges.insert(pos, Some(bond));
        }
    }

    fn is_kekulized(&self) -> bool {
        return !self.delocal_subgraph.is_empty();
    }

    /// Algorithm based on Depth-First article by Richard L. Apodaca
    /// Reference:
    /// https://depth-first.com/articles/2020/02/10/a-comprehensive-treatment-of-aromaticity-in-the-smiles-language/
    fn kekulize(self) -> bool {
        if self.is_kekulized() {
            return true;
        }

        let kept_nodes: HashSet<usize> = self.delocal_subgraph.keys().filter(|node| !self.prune_from_ds(**node)).cloned().collect();

        // relabel kept DS nodes to be 0, 1, 2, ...
        let mut labels: Vec<usize> = kept_nodes.iter().cloned().collect();
        labels.sort();
        let node_to_label: HashMap<usize, usize> = labels.iter().enumerate().map(|(idx, label)| (label.clone(), idx)).collect::<HashMap<_, _>>();

        // pruned and relabelled DS
        let mut pruned_ds: Vec<Vec<usize>> = vec![Vec::new(); labels.len()];
        for node in &kept_nodes {
            let label = node_to_label[&node];
            let adj_nodes: Vec<usize> = self.delocal_subgraph[&node].iter().filter(|v| kept_nodes.contains(v)).cloned().collect();
            for adj in adj_nodes {
                pruned_ds[label].push(node_to_label[&adj]);
            }
        }
        return false;
    }

    fn prune_from_ds(&self, node: usize) -> bool {
        if let Some(adj_nodes) = self.delocal_subgraph.get(&node) {
            let atom = &self.atoms[node];
            let valences = AROMATIC_VALENCES.get(&atom.element.as_str()).unwrap();
            
            // each bond in DS has order 1.5 - we treat them as single bonds
            let mut used_electrons = (self.bond_counts[node] - 0.5 * (adj_nodes.len() as f64)) as i32;

            let valence = last_valence(valences) - atom.charge;
            used_electrons += atom.h_count;
            
            // count the total number of bound electrons of each atom
            let bond_count_int = self.bond_counts[node] as i32;
            let bound_electrons = atom.charge.max(0) + atom.h_count + bond_count_int + (2 * (bond_count_int % 1));
            
            // calculate the number of unpaired electrons of each atom
            let radical_electrons = VALENCE_ELECTRONS[&atom.element.as_str()].max(0) - bound_electrons;
            
            // unpaired electrons do not contribute to the aromatic system
            let free_electrons = valence - used_electrons - radical_electrons;

            let used_electrons_any = match valences {
                ValenceTuple::One(x) => used_electrons == x - atom.charge,
                ValenceTuple::Two(x, y) => used_electrons == x - atom.charge || used_electrons == y - atom.charge,
            };

            if used_electrons_any {
                return true;
            } else {
                return !((free_electrons >= 0) && (free_electrons % 2 != 0))
            }
        } else {
            true
        }
    }
}

#[derive(Debug, Clone)]
pub struct DirectedBond {
    src: usize,
    dst: usize,
    order: f64,
    stereo: Option<String>,
    ring_bond: bool,
}

impl DirectedBond {
    pub fn new(
        src: usize,
        dst: usize,
        order: f64,
        stereo: Option<String>,
        ring_bond: bool,
    ) -> Self {
        DirectedBond {src, dst, order, stereo, ring_bond}
    }
}

#[derive(Debug, PartialEq, Clone)]
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

/// Reads a molecular graph from a SMILES string.
pub fn create_mol_graph(attributable: bool, smiles: &str) -> Result<MolecularGraph, SMILESParserError> {
    if smiles.is_empty() {
        return Err(SMILESParserError::new(smiles.to_string(), "Empty SMILES".to_string(), 0));
    }
    let mut mol = MolecularGraph::new(attributable);
    let mut tokens: Vec<SMILESToken> = SMILESTokenizer::new(smiles).into_iter().collect::<Result<_, _>>()?;

    let q = VecDeque::from(tokens);

    return Ok(mol)
}

fn derive_mol_from_tokens(mol: &mut MolecularGraph, smiles: &str, tokens: &mut Vec<SMILESToken>) {
    
}