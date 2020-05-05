use thiserror::Error;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error {
    #[error("Attempted to access sequence at positon '{0}' when the protein sequence is shorter than this (current length {1})")]
    ProteinTooShort(usize, usize),
}
