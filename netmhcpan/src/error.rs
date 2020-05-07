use thiserror::Error;
use csv::FromUtf8Error;
use nom::error::{ParseError, ErrorKind};


/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error {
    #[error("Attempted to access sequence at positon '{0}' when the protein sequence is shorter than this (current length {1})")]
    ProteinTooShort(usize, usize),
    #[error("Could not open file containing netmhcpan binding data")]
    CouldNotOpenFile(#[from] std::io::Error),
    #[error("Bytes not valid UTF-8.")]
    CouldNotCreateString(#[from] std::string::FromUtf8Error),
}


pub enum ResultParsingError {

}