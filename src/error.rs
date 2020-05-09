use thiserror::Error;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not parse motif positions")]
    IncorrectMotifPositions(#[from] std::num::ParseIntError),
    #[error("Could not load any kir ligand information")]
    KirLigandMapError,
}
