use std::fmt::{Error, Formatter};

#[derive(Debug)]
pub enum RetrieveLigandError {
    InvalidHLA(String),
    RequestError(reqwest::Error),
}

impl std::fmt::Display for RetrieveLigandError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        match self {
            RetrieveLigandError::InvalidHLA(hla) =>
                {writeln!(f, "Please ensure that HLA alleles are of right format (e.g. C*01:102 or C*01:02), the passed HLA was {}", hla)},
            RetrieveLigandError::RequestError(err) => std::fmt::Display::fmt(&err, f)
        }
    }
}

impl From<reqwest::Error> for RetrieveLigandError {
    fn from(e: reqwest::Error) -> Self {
        RetrieveLigandError::RequestError(e)
    }
}
