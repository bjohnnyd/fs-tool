use std::fmt::{Error, Formatter};

#[derive(Debug)]
pub enum RetrieveLigandError {
    InvalidHLA(String),
    RequestError(reqwest::Error),
    NoLigandTableFound(String),
    CSSParseError(u32, u32),
}

impl std::fmt::Display for RetrieveLigandError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        match self {
            RetrieveLigandError::InvalidHLA(hla) =>
                {writeln!(f, "Please ensure that HLA alleles are of right format (e.g. C*01:102 or C*01:02), the passed HLA was {}", hla)},
            RetrieveLigandError::RequestError(err) => std::fmt::Display::fmt(&err, f),
            RetrieveLigandError::NoLigandTableFound(url) =>
                {writeln!(f, "No ligand table found at URL:\n{}", url)},
            RetrieveLigandError::CSSParseError(line, column) =>
            {writeln!(f, "Website parsing error encountered at line {} and character number {}", line, column)},
        }
    }
}

impl From<reqwest::Error> for RetrieveLigandError {
    fn from(e: reqwest::Error) -> Self {
        RetrieveLigandError::RequestError(e)
    }
}

impl std::error::Error for RetrieveLigandError {}
