mod smiles_utils;
mod mol_graph;
mod constants;
mod utilities;

#[pyo3::pymodule]
mod selfies_rust {
    use pyo3::prelude::*;

    use crate::smiles_utils::SMILESTokenizer;

    #[pyfunction]
    #[pyo3(text_signature = "(smiles: str strict: bool attribute: bool) -> list[str]")]
    fn encoder(smiles: String, strict: Option<bool>, attribute: Option<bool>) -> PyResult<Vec<String>> {
        let strict = strict.unwrap_or(true);
        let tokenizer = SMILESTokenizer::new(&smiles);
        
        let mut tokens = Vec::new();
        for token in tokenizer {
            match token {
                Ok(tok) => tokens.push(tok.token),
                Err(e) => eprintln!("{}", e),
            }
        }
        Ok(tokens)
    }
}