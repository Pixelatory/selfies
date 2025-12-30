use std::{fmt, str::Chars};

use once_cell::sync::Lazy;
use regex::Regex;

use crate::mol_graph::{Atom, GraphConstructionError};
use crate::constants::{AROMATIC_SUBSET, ELEMENTS, ORGANIC_SUBSET, SMILES_BOND_ORDERS};
use crate::utilities::{capitalize_first};

static SMILES_BRACKETED_ATOM_PATTERN: Lazy<Regex> = Lazy::new(|| Regex::new(
r"(?x)
    ^[\[]                     # opening square bracket [
    (\d*)                     # isotope number (optional, e.g. 123, 26)
    ([A-Za-z][a-z]?)          # element symbol
    ([@]{0,2})                # chiral_tag (optional, only @ and @@ supported)
    ((?:[H]\d?)?)             # H count (optional, e.g. H, H0, H3)
    ((?:[+]+|[-]+|[+-]\d+)?)  # charge (optional, e.g. ---, +1, ++)
    ((?:[:]\d+)?)             # atom class (optional, e.g. :12, :1)
    []]                       # closing square bracket ]
").unwrap());

// ---------- Enums ----------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SMILESTokenType {
    Atom,
    Branch,
    Ring,
    Dot,
}

// ---------- Token Struct ----------

#[derive(Debug, Clone, PartialEq)]
pub struct SMILESToken {
    pub bond_token: Option<char>,
    pub start_idx: usize,
    pub end_idx: usize,
    pub token_type: SMILESTokenType,
    pub token: String,
}

impl fmt::Display for SMILESToken {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.token)
    }
}

impl SMILESToken {
    pub fn new(
        bond_token: Option<char>,
        start_idx: usize,
        end_idx: usize,
        token_type: SMILESTokenType,
        token: String,
    ) -> Self {
        Self {
            bond_token,
            start_idx,
            end_idx,
            token_type,
            token,
        }
    }
}

// ---------- Tokenizer Iterator ----------

pub struct SMILESTokenizer<'a> {
    smiles: &'a str,
    chars: Chars<'a>,
    i: usize,
}

impl<'a> SMILESTokenizer<'a> {
    pub fn new(smiles: &'a str) -> Self {
        Self {
            smiles,
            chars: smiles.chars(),
            i: 0,
        }
    }
}

impl<'a> Iterator for SMILESTokenizer<'a> {
    type Item = Result<SMILESToken, GraphConstructionError>;

    fn next(&mut self) -> Option<Self::Item> {
        let smiles = self.smiles;
        let chars: Vec<char> = smiles.chars().collect();
        let mut i = self.i;

        while i < chars.len() {
            let ch = chars[i];

            // DOT.
            if ch == '.' {
                let token = SMILESToken {
                    bond_token: None,
                    start_idx: i,
                    end_idx: i + 1,
                    token_type: SMILESTokenType::Dot,
                    token: ch.to_string(),
                };
                self.i = i + 1;
                return Some(Ok(token));
            }

            // BOND.
            let (bond_token, bond_idx) = if SMILES_BOND_ORDERS.contains_key(&ch) {
                i += 1;
                (Some(ch), i-1)
            } else {
                (None, i)
            };

            if i >= chars.len() {
                return Some(Err(GraphConstructionError::UnexpectedToken {
                    smiles: smiles.to_string(),
                    message: "Hanging bond.".to_string(),
                    index: i.saturating_sub(1),
                }));
            }

            let ch = chars[i];
            let token: SMILESToken;

            // ATOMS (one or two letters).
            if ch.is_ascii_alphabetic() {
                let next_two = if i + 2 <= smiles.len() {
                    &smiles[i..i + 2]
                } else {
                    ""
                };
                if next_two == "Br" || next_two == "Cl" {
                    token = SMILESToken {
                        bond_token,
                        start_idx: i,
                        end_idx: i + 2,
                        token_type: SMILESTokenType::Atom,
                        token: next_two.to_string(),
                    };
                } else {
                    token = SMILESToken {
                        bond_token,
                        start_idx: i,
                        end_idx: i + 1,
                        token_type: SMILESTokenType::Atom,
                        token: ch.to_string(),
                    };
                }

            // BRACKETED ATOMS.
            } else if ch == '[' {
                if let Some(r_idx) = smiles[i + 1..].find(']') {
                    let end_idx = i + r_idx + 2;
                    token = SMILESToken {
                        bond_token,
                        start_idx: i,
                        end_idx,
                        token_type: SMILESTokenType::Atom,
                        token: smiles[i..end_idx].to_string(),
                    };
                } else {
                    return Some(Err(GraphConstructionError::UnexpectedToken {
                        smiles: smiles.to_string(),
                        message: "Hanging bracket '['.".to_string(),
                        index: i,
                    }));
                }

            // BRANCHES.
            } else if ch == '(' || ch == ')' {
                if bond_token.is_some() {
                    return Some(Err(GraphConstructionError::UnexpectedToken {
                        smiles: smiles.to_string(),
                        message: "Hanging bond.".to_string(),
                        index: bond_idx,
                    }));
                }
                token = SMILESToken {
                    bond_token: None,
                    start_idx: i,
                    end_idx: i + 1,
                    token_type: SMILESTokenType::Branch,
                    token: ch.to_string(),
                };

            // RINGS.
            } else if ch.is_ascii_digit() {
                token = SMILESToken {
                    bond_token,
                    start_idx: i,
                    end_idx: i + 1,
                    token_type: SMILESTokenType::Ring,
                    token: ch.to_string(),
                };

            } else if ch == '%' {
                if i + 3 <= smiles.len() {
                    let rnum = &smiles[i + 1..i + 3];
                    if rnum.chars().all(|c| c.is_ascii_digit()) {
                        token = SMILESToken {
                            bond_token,
                            start_idx: i,
                            end_idx: i + 3,
                            token_type: SMILESTokenType::Ring,
                            token: smiles[i..i + 3].to_string(),
                        };
                    } else {
                        return Some(Err(GraphConstructionError::UnexpectedToken {
                            smiles: smiles.to_string(),
                            message: format!("Invalid ring number '%{}'.", rnum),
                            index: i,
                        }));
                    }
                } else {
                    return Some(Err(GraphConstructionError::UnexpectedToken {
                        smiles: smiles.to_string(),
                        message: "Incomplete ring number.".to_string(),
                        index: i,
                    }));
                }

            // UNKNOWN SYMBOL.
            } else {
                return Some(Err(GraphConstructionError::UnexpectedToken {
                    smiles: smiles.to_string(),
                    message: format!("Unrecognized symbol '{}'.", ch),
                    index: i,
                }));
            }

            self.i = token.end_idx;
            return Some(Ok(token));
        }

        None
    }
}


/// Reads an atom from its SMILES representation.
pub fn smiles_to_atom(atom_symbol: &str) -> Option<Atom> {
    let chars: Vec<char> = atom_symbol.chars().collect();

    if ORGANIC_SUBSET.contains(atom_symbol) {
        return Some(Atom::new(atom_symbol.to_string(),  false, None, None, None, 0));
    } else if AROMATIC_SUBSET.contains(atom_symbol) {
        return Some(Atom::new(atom_symbol.to_string().to_uppercase(), true, None, None, None, 0));
    } else if !(*chars.first()? == '[' && *chars.last()? == ']') {
        return None;
    }

    let matches: Vec<(&str, &str, &str, &str, &str, &str)> = SMILES_BRACKETED_ATOM_PATTERN
        .captures_iter(atom_symbol)
        .map(|caps| 
        (
            caps.get(1).map_or("", |m| m.as_str()), // isotope
            caps.get(2).map_or("", |m| m.as_str()), // element
            caps.get(3).map_or("", |m| m.as_str()), // chiral tag
            caps.get(4).map_or("", |m| m.as_str()), // hydrogens
            caps.get(5).map_or("", |m| m.as_str()), // charge
            caps.get(6).map_or("", |m| m.as_str()), // atom class
        )
    )
        .collect();
    
    if matches.is_empty() {
        return None;
    }
    let (
        isotope_str,
        element_match,
        chirality_match,
        h_count_match,
        charge_match,
        _,
    ) = matches.first()?;
    

    let isotope = match isotope_str.is_empty() {
        false => Some(isotope_str.parse::<u32>().ok()?),
        true => None,
    };
    let is_aromatic = element_match.chars().all(|c| c.is_ascii_lowercase())
        && AROMATIC_SUBSET.contains(element_match);
    let element = capitalize_first(element_match);
    if !ELEMENTS.contains(&element.as_str()) {
        return None;
    }

    let chirality = match chirality_match.is_empty() {
        false => Some(chirality_match.to_string()),
        true => None,
    };

    let h_count = match h_count_match.is_empty() {
        false => {
            let h_count_suffix = h_count_match.strip_prefix("H")?;
            match h_count_suffix.is_empty() {
                false => Some(h_count_suffix.parse::<u32>().ok()?),
                true => Some(1),
            }
        },
        true => Some(0),
    };

    let charge = match charge_match.is_empty() {
        false => {
            let mut tmp_charge;
            if charge_match.ends_with(|ch: char| ch.is_ascii_digit()) {
                tmp_charge = charge_match
                    .strip_prefix(|ch: char| ch == '+' || ch == '-')?
                    .parse::<i32>()
                    .ok()?;
            } else {
                tmp_charge = charge_match.len().try_into().ok()?;
            }
            if charge_match.starts_with('-') {
                tmp_charge *= -1;
            }
            tmp_charge
        },
        true => 0,
    };
    
    Some(
        Atom::new(element, is_aromatic, isotope, chirality, h_count, charge)
    )
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_smiles_to_atom() {
        // Element, h_count, and charge are filled out.
        let x = smiles_to_atom("[OH3+]");
        assert!(x.is_some());
        assert_eq!(
            x.unwrap(),
            Atom::new(String::from("O"), false, None, None, Some(3), 1)
        );
        
        // Charge counts the many negative symbols.
        let x = smiles_to_atom("[Co---]");
        assert!(x.is_some());
        assert_eq!(
            x.unwrap(),
            Atom::new(String::from("Co"), false, None, None, Some(0), -3)
        );

        // Isotope and is_aromatic are filled out.
        let x = smiles_to_atom("[14cH]");
        assert!(x.is_some());
        assert_eq!(
            x.unwrap(),
            Atom::new(String::from("C"), true, Some(14), None, Some(1), 0)
        );

        // Chirality is filled out.
        let x = smiles_to_atom("[C@@H]");
        assert!(x.is_some());
        assert_eq!(
            x.unwrap(),
            Atom::new(String::from("C"), false, None, Some(String::from("@@")), Some(1), 0)
        );

        // Unknown element.
        let x = smiles_to_atom("[Zh++]");
        assert!(x.is_none());
        
        // Aromatic element.
        let x = smiles_to_atom("c");
        assert!(x.is_some());
        assert_eq!(
            x.unwrap(),
            Atom::new(String::from("C"), true, None, None, None, 0)
        );
    }

    #[test]
    fn test_smiles_tokenizer() {
        // The molecule is not structurally sound, but let's care more about testing the tokenization logic
        // than the molecular validity.
        let x = SMILESTokenizer::new("C(Br)C=O.N#C1CC1Cl[2FeH+++]%10CC%10");
        let tokens: Vec<SMILESToken> = x.into_iter().collect::<Result<Vec<_>, _>>().unwrap();
        let expected_tokens = Vec::from(
            [
                SMILESToken{bond_token: None, start_idx: 0, end_idx: 1, token_type: SMILESTokenType::Atom, token: String::from("C")},
                SMILESToken{bond_token: None, start_idx: 1, end_idx: 2, token_type: SMILESTokenType::Branch, token: String::from("(")},
                SMILESToken{bond_token: None, start_idx: 2, end_idx: 4, token_type: SMILESTokenType::Atom, token: String::from("Br")},
                SMILESToken{bond_token: None, start_idx: 4, end_idx: 5, token_type: SMILESTokenType::Branch, token: String::from(")")},
                SMILESToken{bond_token: None, start_idx: 5, end_idx: 6, token_type: SMILESTokenType::Atom, token: String::from("C")},
                SMILESToken{bond_token: Some('='), start_idx: 7, end_idx: 8, token_type: SMILESTokenType::Atom, token: String::from("O")},
                SMILESToken{bond_token: None, start_idx: 8, end_idx: 9, token_type: SMILESTokenType::Dot, token: String::from(".")},
                SMILESToken{bond_token: None, start_idx: 9, end_idx: 10, token_type: SMILESTokenType::Atom, token: String::from("N")},
                SMILESToken{bond_token: Some('#'), start_idx: 11, end_idx: 12, token_type: SMILESTokenType::Atom, token: String::from("C")},
                SMILESToken{bond_token: None, start_idx: 12, end_idx: 13, token_type: SMILESTokenType::Ring, token: String::from("1")},
                SMILESToken{bond_token: None, start_idx: 13, end_idx: 14, token_type: SMILESTokenType::Atom, token: String::from("C")},
                SMILESToken{bond_token: None, start_idx: 14, end_idx: 15, token_type: SMILESTokenType::Atom, token: String::from("C")},
                SMILESToken{bond_token: None, start_idx: 15, end_idx: 16, token_type: SMILESTokenType::Ring, token: String::from("1") },
                SMILESToken{bond_token: None, start_idx: 16, end_idx: 18, token_type: SMILESTokenType::Atom, token: String::from("Cl") },
                SMILESToken{bond_token: None, start_idx: 18, end_idx: 27, token_type: SMILESTokenType::Atom, token: String::from("[2FeH+++]") },
                SMILESToken{bond_token: None, start_idx: 27, end_idx: 30, token_type: SMILESTokenType::Ring, token: String::from("%10") },
                SMILESToken{bond_token: None, start_idx: 30, end_idx: 31, token_type: SMILESTokenType::Atom, token: String::from("C") },
                SMILESToken{bond_token: None, start_idx: 31, end_idx: 32, token_type: SMILESTokenType::Atom, token: String::from("C") },
                SMILESToken{bond_token: None, start_idx: 32, end_idx: 35, token_type: SMILESTokenType::Ring, token: String::from("%10") },
            ]
        );
        assert_eq!(tokens, expected_tokens);
    }

    #[test]
    fn test_smiles_tokenizer_2() {
        let x = SMILESTokenizer::new("F/C=C/F");
        let tokens: Vec<SMILESToken> = x.into_iter().collect::<Result<Vec<_>, _>>().unwrap();
        let expected_tokens = Vec::from(
            [
                SMILESToken::new(None, 0, 1, SMILESTokenType::Atom, "F".to_string()),
                SMILESToken::new(Some('/'), 2, 3, SMILESTokenType::Atom, "C".to_string()),
                SMILESToken::new(Some('='), 4, 5, SMILESTokenType::Atom, "C".to_string()),
                SMILESToken::new(Some('/'), 6, 7, SMILESTokenType::Atom, "F".to_string()),
            ]
        );
        assert_eq!(tokens, expected_tokens);
    }
}