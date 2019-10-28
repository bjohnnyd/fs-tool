use std::error;
use std::fmt;

type Result<T> = std::result::Result<T, HLAErr::HLAError>;


pub mod HLAErr  {
use super::*;

#[derive(Debug, Clone, Eq, PartialEq, Copy)]
pub (crate) struct HLAError;

impl fmt::Display for HLAError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Invalid HLA locus/gene specified the allowd options are A, B, C, DP, DR, DQ, DN")
    }
}

impl error::Error for HLAError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        None
    }
}
}


#[cfg(test)]
mod tests {
    #[test]
    fn check_hla_error() {
    }
}
    