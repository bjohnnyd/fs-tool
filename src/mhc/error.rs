use std::fmt::{self,Display};

type Result<T> = std::result::Result<T, HLAErr::HLAError>;

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub enum HLAErr {
    ParseError(String)
}

impl fmt::Display for HLAErr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            HLAErr::ParseError(s) => {write!(f, "HLA allele specified incorrectly: {}", s)}
        }
    }
}

impl std::error::Error for HLAErr{
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            HLAErr::ParseError(_) => None,
        }
    }
}
//impl From<std::io::Error> for Fast5Error {
//        fn from(err: std::io::Error) -> Self {
//            Fast5Error::PathError(err)
//        }
//
//    }

