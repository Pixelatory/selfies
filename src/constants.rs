use std::{collections::{HashMap, HashSet}};

use once_cell::sync::Lazy;

#[derive(Debug)]
pub enum ValenceTuple {
    One(i32),
    Two(i32, i32),
}

pub static ELEMENTS: Lazy<HashSet<&'static str>> = Lazy::new(|| HashSet::from(
    [
        "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
        "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
        "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
        "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
        "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "Hf",
        "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
        "Po", "At", "Rn", "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
        "Ds", "Rg", "Cn", "Fl", "Lv", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
        "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Ac", "Th",
        "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
        "No", "Lr"
    ]
));


pub static ORGANIC_SUBSET: Lazy<HashSet<&'static str>> = Lazy::new(|| HashSet::from(["B", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]));

pub static AROMATIC_VALENCES: Lazy<HashMap<&'static str, ValenceTuple>> = Lazy::new(|| HashMap::from(
    [
        ("B", ValenceTuple::One(3)),
        ("Al", ValenceTuple::One(3)),
        ("C", ValenceTuple::One(4)),
        ("Si", ValenceTuple::One(4)),
        ("N", ValenceTuple::Two(3, 5)),
        ("P", ValenceTuple::Two(3, 5)),
        ("As", ValenceTuple::Two(3, 5)),
        ("O", ValenceTuple::Two(2, 4)),
        ("S", ValenceTuple::Two(2, 4)),
        ("Se", ValenceTuple::Two(2, 4)),
        ("Te", ValenceTuple::Two(2, 4)),
    ]
));

pub static VALENCE_ELECTRONS: Lazy<HashMap<&'static str, i32>> = Lazy::new(|| HashMap::from(
    [
        ("B", 3),
        ("Al", 3),
        ("C", 3),
        ("Si", 4),
        ("N", 5),
        ("P", 5),
        ("As", 5),
        ("O", 6),
        ("S", 6),
        ("Se", 6),
        ("Te", 6),
    ]
));

// AROMATIC_SUBSET should have the same keys as AROMATIC_VALENCES, but be lowercase.
pub static AROMATIC_SUBSET: Lazy<HashSet<&'static str>> = Lazy::new(|| HashSet::from(
    ["b", "al", "c", "si", "n", "p", "as", "o", "s", "se", "te"]
));