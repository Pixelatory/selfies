use std::collections::HashMap;

use once_cell::sync::Lazy;

use crate::{mol_graph::{GraphConstructionError, MolecularGraph}, smiles_utils::atom_to_smiles};

static UNKNOWN_CONSTRAINT: u32 = 8;
static CURRENT_CONSTRAINTS: Lazy<HashMap<String, u32>> = Lazy::new(|| HashMap::from(
    [
        ("H".to_string(), 1),
        ("F".to_string(), 1),
        ("Cl".to_string(), 1),
        ("Br".to_string(), 1),
        ("I".to_string(), 1),
        ("B".to_string(), 3),
        ("B+1".to_string(), 2),
        ("B-1".to_string(), 4),
        ("O".to_string(), 2),
        ("O+1".to_string(), 3),
        ("O-1".to_string(), 1),
        ("N".to_string(), 3),
        ("N+1".to_string(), 4),
        ("N-1".to_string(), 2),
        ("C".to_string(), 4),
        ("C+1".to_string(), 3),
        ("C-1".to_string(), 3),
        ("P".to_string(), 5),
        ("P+1".to_string(), 4),
        ("P-1".to_string(), 6),
        ("S".to_string(), 6),
        ("S+1".to_string(), 5),
        ("S-1".to_string(), 5),
    ]
));

fn get_bonding_capacity(element: &str, charge: i32) -> u32 {
    let mut key = String::from(element);
    if charge > 0 {
        key += &format!("+{}", charge);
    } else if charge < 0 {
        key += &charge.to_string();
    }

    if let Some(current_constraint) = CURRENT_CONSTRAINTS.get(&key) {
        return *current_constraint;
    }
    return UNKNOWN_CONSTRAINT;
}

/// Check that each atom in the molecular graph has its bond constraints met.
/// 
/// An atom's bond count cannot exceed its bonding capacity.
pub fn check_bond_constraints(mol: &MolecularGraph, smiles: &str) -> Result<(), GraphConstructionError> {
    let mut troubled_smiles = Vec::new();
    for (atom_index, atom) in mol.get_atoms().iter().enumerate() {
        let bond_cap = get_bonding_capacity(&atom.element, atom.charge) as f64;
        let bond_count = mol.get_bond_count(atom_index);
        if bond_count > bond_cap {
            troubled_smiles.push((atom_to_smiles(atom, true), bond_count, bond_cap));
        }
    }

    if !troubled_smiles.is_empty() {
        let mut error_message = String::from("Input violates the currently set semantic constraints.\n");
        error_message += &format!("\tSMILES: {smiles}\n");
        error_message += &String::from("\tErrors:\n");

        for (smiles, bond_count, bond_cap) in troubled_smiles {
            error_message += &format!(
                "\t[{} with {} bond(s) - a max of {} bond(s) is allowed]",
                smiles, bond_count, bond_cap,
            );
        }
        return Err(GraphConstructionError::InvalidBondConstraints { message: error_message });
    }

    Ok(())
}