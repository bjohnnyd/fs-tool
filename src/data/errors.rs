use std::fmt::{Error, Formatter};

#[derive(Debug)]
pub enum RetrieveLigandError {
    InvalidHLA(String),
    WebsiteError(attohttpc::Error),
    NoLigandTableFound(String),
    CSSParseError(u32, u32),
    FileError(std::io::Error),
    ErrorWithIPDWebsite(String),
    CouldNotAccessData,
}

impl std::fmt::Display for RetrieveLigandError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        match self {
            RetrieveLigandError::InvalidHLA(hla) =>
                {writeln!(f, "Please ensure that HLA alleles are of right format (e.g. C*01:102 or C*01:02), the passed HLA was {}", hla)},
            RetrieveLigandError::WebsiteError(err) => std::fmt::Display::fmt(&err, f),
            RetrieveLigandError::NoLigandTableFound(url) =>
                {writeln!(f, "No ligand table found at URL:\n{}", url)},
            RetrieveLigandError::CSSParseError(line, column) =>
            {writeln!(f, "Website parsing error encountered at line {} and character number {}", line, column)},
            RetrieveLigandError::FileError(err) => std::fmt::Display::fmt(&err, f),
            RetrieveLigandError::ErrorWithIPDWebsite(url) =>
                {writeln!(f, "Error with accessing IPD Website at {}", url)},
            RetrieveLigandError::CouldNotAccessData =>
                {writeln!(f, "Error with accessing local resources")},
        }
    }
}

impl From<std::io::Error> for RetrieveLigandError {
    fn from(e: std::io::Error) -> Self {
        RetrieveLigandError::FileError(e)
    }
}

impl From<attohttpc::Error> for RetrieveLigandError {
    fn from(e: attohttpc::Error) -> Self {
        RetrieveLigandError::WebsiteError(e)
    }
}
impl std::error::Error for RetrieveLigandError {}
