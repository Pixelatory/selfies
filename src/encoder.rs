use crate::{errors::EncoderError, mol_graph::{Atom, DirectedBond}, smiles_utils::{atom_to_smiles, bond_to_smiles}};

fn bond_to_selfies(smiles: &str, bond: &DirectedBond, show_stereo: bool) -> Result<String, EncoderError> {
    if !show_stereo && bond.order == 1.0 {
        return Ok(String::new());
    }
    return bond_to_smiles(smiles, bond);
}

fn get_stereo_or_single(bond: &DirectedBond) -> char {
    if let Some(stereo) = bond.stereo {
        return stereo;
    } else {
        return '-';
    }
}

/// Convert a ring bond to a SELFIES string.
fn ring_bonds_to_selfies(smiles: &str, l_bond: &DirectedBond, r_bond: &DirectedBond) -> Result<String, EncoderError> {
    if !l_bond.ring_bond {
        return Err(
            EncoderError::UnexpectedBondOrder {
                smiles: smiles.to_string(),
                message: String::from("Attempting to convert a ring bond to SELFIES which is not a ring bond."),
                src: l_bond.src,
                dst: l_bond.dst,
            }
        )
    }
    if !r_bond.ring_bond {
        return Err(
            EncoderError::UnexpectedBondOrder {
                smiles: smiles.to_string(),
                message: String::from("Attempting to convert a ring bond to SELFIES which is not a ring bond."),
                src: r_bond.src,
                dst: r_bond.dst,
            }
        )
    }
    if l_bond.order != r_bond.order {
        return Err(
            EncoderError::UnexpectedBondOrder {
                smiles: smiles.to_string(),
                message: format!("Left and right bonds do not have matching orders: ({}, {})", l_bond.order, r_bond.order),
                src: l_bond.src,
                dst: l_bond.dst,
            }
        );
    }

    if l_bond.order != 1.0 || (l_bond.stereo.is_none() && r_bond.stereo.is_none()) {
        return bond_to_selfies(smiles, l_bond, false);
    } else {
        let mut bond_char = get_stereo_or_single(l_bond).to_string();
        bond_char += &get_stereo_or_single(r_bond).to_string();
        return Ok(bond_char);
    }
}

/// Convert an Atom and Bond to a SELFIES string.
fn atom_to_selfies(smiles: &str, atom_index: usize, maybe_bond: Option<&DirectedBond>, atom: &Atom) -> Result<String, EncoderError> {
    if atom.is_aromatic {
        return Err(
            EncoderError::AromaticAtom {
                smiles: smiles.to_string(),
                message: String::from(""),
                atom_index: atom_index,
            }
        );
    }

    let bond_str: String;
    if let Some(bond) = maybe_bond {
        bond_str = bond_to_selfies(smiles, bond, true)?;
    } else {
        bond_str = String::from("");
    }

    let atom_str = atom_to_smiles(atom, false);

    Ok(format!("[{}{}]", bond_str, atom_str))
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_to_selfies() {
        // Aromatic atoms are not allowed. Kekulization must occur beforehand.
        let atom = Atom::new(String::from("C"), true, None, None, None, 0);
        let x = atom_to_selfies("c", 0, None, &atom);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::AromaticAtom {smiles: _, message: _, atom_index: _}));

        // Bond has an order that marks aromaticness. This is not allowed.
        // Typically the atom would be marked as aromatic too, but this is just for testing.
        let atom = Atom::new(String::from("C"), false, None, None, None, 0);
        let bond = DirectedBond::new(0, 1, 1.5, None, false);
        let x = atom_to_selfies("cc", 0, Some(&bond), &atom);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::UnexpectedBondOrder {smiles: _, message: _, src: _, dst: _}));

        // Single bond with carbon.
        let atom = Atom::new(String::from("C"), false, None, None, None, 0);
        let bond = DirectedBond::new(0, 1, 1.0, None, false);
        let x = atom_to_selfies("CC", 0, Some(&bond), &atom);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("[C]"));

        // Double bond with nitrogen.
        let atom = Atom::new(String::from("N"), false, None, None, None, 0);
        let bond = DirectedBond::new(0, 1, 2.0, None, false);
        let x = atom_to_selfies("N=N", 0, Some(&bond), &atom);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("[=N]"));

        // Triple bond with flourine.
        let atom = Atom::new(String::from("F"), false, None, None, None, 0);
        let bond = DirectedBond::new(0, 1, 3.0, None, false);
        let x = atom_to_selfies("F#F", 0, Some(&bond), &atom);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("[#F]"));
    }

    #[test]
    fn test_ring_bonds_to_selfies() {
        // l_bond is has ring_bond set to false.
        let l_bond = DirectedBond::new(0, 5, 1.0, None, false);
        let r_bond = DirectedBond::new(5, 0, 1.0, None, true);
        let x = ring_bonds_to_selfies("c1ccccc1", &l_bond, &r_bond);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::UnexpectedBondOrder { smiles: _, message: _, src: 0, dst: 5 }));

        // r_bond is has ring_bond set to false.
        let l_bond = DirectedBond::new(0, 5, 1.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 1.0, None, false);
        let x = ring_bonds_to_selfies("c1ccccc1", &l_bond, &r_bond);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::UnexpectedBondOrder { smiles: _, message: _, src: 5, dst: 0 }));

        // l_bond and r_bond have mismatching orders.
        let l_bond = DirectedBond::new(0, 5, 1.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 2.0, None, true);
        let x = ring_bonds_to_selfies("c1ccccc1", &l_bond, &r_bond);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::UnexpectedBondOrder { smiles: _, message: _, src: 0, dst: 5 }));
        
        // Both bonds have no stereo, and order 1.0.
        let l_bond = DirectedBond::new(0, 5, 1.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 1.0, None, true);
        let x = ring_bonds_to_selfies("C1CCCCC1", &l_bond, &r_bond);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from(""));

        // l_bond has stereo.
        let l_bond = DirectedBond::new(0, 5, 1.0, Some('/'), true);
        let r_bond = DirectedBond::new(5, 0, 1.0, None, true);
        let x = ring_bonds_to_selfies("C1CCCCC1", &l_bond, &r_bond);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("/-"));

        // r_bond has stereo.
        let l_bond = DirectedBond::new(0, 5, 1.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 1.0, Some('\\'), true);
        let x = ring_bonds_to_selfies("C1CCCCC1", &l_bond, &r_bond);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("-\\"));

        // Both bonds have order 2.0.
        let l_bond = DirectedBond::new(0, 5, 2.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 2.0, None, true);
        let x = ring_bonds_to_selfies("C1CCCCC1", &l_bond, &r_bond);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("="));

        // Both bonds have order 5.0, which is unknown.
        let l_bond = DirectedBond::new(0, 5, 5.0, None, true);
        let r_bond = DirectedBond::new(5, 0, 5.0, None, true);
        let x = ring_bonds_to_selfies("C1CCCCC1", &l_bond, &r_bond);
        assert!(x.is_err());
        assert!(matches!(x.unwrap_err(), EncoderError::UnexpectedBondOrder { smiles: _, message: _, src: _, dst: _ }))
    }

    #[test]
    fn test_bond_to_selfies() {
        // Single bond, and do not show stereo.
        let bond = DirectedBond::new(0, 1, 1.0, Some('/'), false);
        let x = bond_to_selfies("C", &bond, false);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from(""));

        // Single bond, and show stereo.
        let bond = DirectedBond::new(0, 1, 1.0, Some('/'), false);
        let x = bond_to_selfies("C", &bond, true);
        assert!(x.is_ok());
        assert_eq!(x.unwrap(), String::from("/"));
    }
}