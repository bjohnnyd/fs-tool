use thiserror::Error;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not parse motif positions")]
    IncorrectMotifPositions(#[from] std::num::ParseIntError),
}
