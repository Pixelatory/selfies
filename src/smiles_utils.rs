use std::{fmt, str::Chars};

// ---------- Enums ----------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SMILESTokenType {
    Atom,
    Branch,
    Ring,
    Dot,
}

// ---------- Token Struct ----------

#[derive(Debug, Clone)]
pub struct SMILESToken {
    pub bond_idx: Option<usize>,
    pub start_idx: usize,
    pub end_idx: usize,
    pub token_type: SMILESTokenType,
    pub token: String,
}

impl SMILESToken {
    pub fn extract_bond_char(&self, smiles: &str) -> Option<char> {
        self.bond_idx.and_then(|i| smiles.chars().nth(i))
    }

    pub fn extract_symbol<'a>(&self, smiles: &'a str) -> &'a str {
        &smiles[self.start_idx..self.end_idx]
    }
}

impl fmt::Display for SMILESToken {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.token)
    }
}

// ---------- Error Handling ----------

#[derive(Debug)]
pub struct SMILESParserError {
    pub smiles: String,
    pub message: String,
    pub index: usize,
}

impl fmt::Display for SMILESParserError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SMILESParserError at {}: {} (input: '{}')",
            self.index, self.message, self.smiles
        )
    }
}

impl std::error::Error for SMILESParserError {}

// ---------- Constants ----------

const SMILES_BOND_ORDERS: &[char] = &['-', '=', '#', ':', '/', '\\'];

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
    type Item = Result<SMILESToken, SMILESParserError>;

    fn next(&mut self) -> Option<Self::Item> {
        let smiles = self.smiles;
        let chars: Vec<char> = smiles.chars().collect();
        let mut i = self.i;

        while i < chars.len() {
            let ch = chars[i];

            // DOT.
            if ch == '.' {
                let token = SMILESToken {
                    bond_idx: None,
                    start_idx: i,
                    end_idx: i + 1,
                    token_type: SMILESTokenType::Dot,
                    token: ch.to_string(),
                };
                self.i = i + 1;
                return Some(Ok(token));
            }

            // BOND.
            let bond_idx = if SMILES_BOND_ORDERS.contains(&ch) {
                i += 1;
                Some(i - 1)
            } else {
                None
            };

            if i >= chars.len() {
                return Some(Err(SMILESParserError {
                    smiles: smiles.to_string(),
                    message: "hanging bond".to_string(),
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
                        bond_idx,
                        start_idx: i,
                        end_idx: i + 2,
                        token_type: SMILESTokenType::Atom,
                        token: next_two.to_string(),
                    };
                } else {
                    token = SMILESToken {
                        bond_idx,
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
                        bond_idx,
                        start_idx: i,
                        end_idx,
                        token_type: SMILESTokenType::Atom,
                        token: smiles[i..end_idx].to_string(),
                    };
                } else {
                    return Some(Err(SMILESParserError {
                        smiles: smiles.to_string(),
                        message: "hanging bracket [".to_string(),
                        index: i,
                    }));
                }

            // BRANCHES.
            } else if ch == '(' || ch == ')' {
                if bond_idx.is_some() {
                    return Some(Err(SMILESParserError {
                        smiles: smiles.to_string(),
                        message: "hanging bond".to_string(),
                        index: bond_idx.unwrap(),
                    }));
                }
                token = SMILESToken {
                    bond_idx: None,
                    start_idx: i,
                    end_idx: i + 1,
                    token_type: SMILESTokenType::Branch,
                    token: ch.to_string(),
                };

            // RINGS.
            } else if ch.is_ascii_digit() {
                token = SMILESToken {
                    bond_idx,
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
                            bond_idx,
                            start_idx: i,
                            end_idx: i + 3,
                            token_type: SMILESTokenType::Ring,
                            token: smiles[i..i + 3].to_string(),
                        };
                    } else {
                        return Some(Err(SMILESParserError {
                            smiles: smiles.to_string(),
                            message: format!("invalid ring number '%{}'", rnum),
                            index: i,
                        }));
                    }
                } else {
                    return Some(Err(SMILESParserError {
                        smiles: smiles.to_string(),
                        message: "incomplete ring number".to_string(),
                        index: i,
                    }));
                }

            // UNKNOWN SYMBOL.
            } else {
                return Some(Err(SMILESParserError {
                    smiles: smiles.to_string(),
                    message: format!("unrecognized symbol '{}'", ch),
                    index: i,
                }));
            }

            self.i = token.end_idx;
            return Some(Ok(token));
        }

        None
    }
}