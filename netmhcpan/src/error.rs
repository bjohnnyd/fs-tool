use thiserror::Error;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error<I>
where
    I: std::fmt::Debug,
{
    #[error("Attempted to access sequence at positon '{0}' when the protein sequence is shorter than this (current length {1})")]
    ProteinTooShort(usize, usize),
    #[error("Could not open file containing netmhcpan binding data")]
    CouldNotOpenFile(#[from] std::io::Error),
    #[error("Bytes not valid UTF-8.")]
    CouldNotCreateString(#[from] std::string::FromUtf8Error),
    #[error("Could not parse binding data in line: {0:?}. Issue with  parser {1:?}")]
    ParseError(I, nom::error::ErrorKind),
}

