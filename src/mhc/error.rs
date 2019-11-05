use crate::mhc::hla::HLA;
use std::fmt::{self, Display};

type Result<T> = std::result::Result<T, HLAErr>;

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum HLAErr {
    ParseError(String),
    GeneNameTooShort,
}

impl fmt::Display for HLAErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            HLAErr::ParseError(s) => write!(f, "HLA allele specified incorrectly: {}", s),
            HLAErr::GeneNameTooShort => write!(f, "The name of the hla gene name was too short"),
        }
    }
}

impl std::error::Error for HLAErr {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            HLAErr::ParseError(_) => None,
            HLAErr::GeneNameTooShort => None,
        }
    }
}
//impl From<std::option::NoneError> for HLAErr {
//        fn from(err: std::option::NoneError) -> Self {
//            HLAErr::GeneNameTooShort
//        }
//    }
