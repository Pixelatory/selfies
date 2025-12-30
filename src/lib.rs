mod smiles_utils;
mod mol_graph;
mod constants;
mod utilities;
mod matching_utils;

#[pyo3::pymodule]
mod selfies_rust {
    use pyo3::prelude::*;

    use crate::{mol_graph::create_mol_graph};

    #[pyfunction]
    #[pyo3(text_signature = "(smiles: str strict: bool attribute: bool) -> list[str]")]
    fn encoder(smiles: String, strict: Option<bool>, attribute: Option<bool>) -> PyResult<Vec<String>> {
        let mol = create_mol_graph(&smiles, true).unwrap();

        return Ok(Vec::new());
    }
}