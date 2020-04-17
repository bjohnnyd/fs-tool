use std::iter::FromIterator;
use std::str::FromStr;

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum Gene {
    A,
    B,
    C,
    DP,
    DM,
    DO,
    DQ,
    DR,
    Unknown,
}

impl Gene {
    pub fn is_unknown(&self) -> bool {
        self.eq(&Gene::Unknown)
    }
}

impl From<&Gene> for &str {
    fn from(g: &Gene) -> Self {
        match g {
            Gene::A => "A",
            Gene::B => "B",
            Gene::C => "C",
            Gene::DP => "DP",
            Gene::DM => "DM",
            Gene::DO => "DO",
            Gene::DQ => "DQ",
            Gene::DR => "DR",
            Gene::Unknown => "",
        }
    }
}

impl FromIterator<char> for Gene {
    fn from_iter<T: IntoIterator<Item = char>>(iter: T) -> Self {
        let mut iter = iter.into_iter();
        match iter.next() {
            Some('A') => Gene::A,
            Some('B') => Gene::B,
            Some('C') => Gene::C,
            Some('D') => match iter.next() {
                Some('P') => Gene::DP,
                Some('M') => Gene::DM,
                Some('O') => Gene::DO,
                Some('Q') => Gene::DQ,
                Some('R') => Gene::DR,
                _ => Gene::Unknown,
            },
            _ => Gene::Unknown,
        }
    }
}

impl std::fmt::Display for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use Gene::*;

        let s = match self {
            A => "A",
            B => "B",
            C => "C",
            DP => "DP",
            DM => "DM",
            DO => "DO",
            DQ => "DQ",
            DR => "DR",
            Unknown => "Unknown",
        }
            .to_string();

        write!(f, "{}", s)
    }
}


