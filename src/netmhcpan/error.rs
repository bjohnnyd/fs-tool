use std::error;
use std::fmt;

type Result<T> = std::result::Result<T, HLAErr::HLAError>;

#[derive(Default)]
pub struct ParseError;

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "A parsing error occurred.")
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn check_netmhcpan_error() {}
}
