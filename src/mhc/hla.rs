use crate::mhc::errors::HLAError;
use std::iter::FromIterator;
use std::str::FromStr;

type Result<T> = std::result::Result<T, HLAError>;
#[derive(Debug, Eq, PartialEq)]
pub struct HLA {
    pub gene: Gene,
    pub allele_group: String,
    pub hla_protein: Option<String>,
    pub cds_synonymous_sub: Option<String>,
    pub non_coding_diff: Option<String>,
    pub expression_change: ExpressionChange,
}

#[derive(Debug, Eq, PartialEq)]
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
    fn is_unknown(&self) -> bool {
        self.eq(&Gene::Unknown)
    }
}
#[derive(Debug, Eq, PartialEq)]
pub enum ExpressionChange {
    N,
    L,
    S,
    C,
    A,
    Q,
    Unknown,
}

impl From<char> for ExpressionChange {
    fn from(c: char) -> Self {
        match c {
            'N' => ExpressionChange::N,
            'L' => ExpressionChange::L,
            'S' => ExpressionChange::S,
            'C' => ExpressionChange::C,
            'A' => ExpressionChange::A,
            'Q' => ExpressionChange::Q,
            _ => ExpressionChange::Unknown,
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

impl HLA {
    pub fn new<T: AsRef<str>>(s: T) -> Result<Self> {
        let mut hla_name = s.as_ref().trim_start_matches("HLA-");
        hla_name.parse::<HLA>()
    }
}

macro_rules! to_option {
    ($s:ident) => {{
        Some($s).filter(|s| !s.is_empty())
    }};
}

impl std::str::FromStr for HLA {
    type Err = HLAError;

    fn from_str(s: &str) -> Result<Self> {
        let hla = s
            .trim_start_matches("HLA-")
            .replace("*", "")
            .replace(":", "");

        let gene: Gene = hla.chars().take(2).collect::<Gene>();

        if gene.is_unknown() {
            return Err(HLAError::IncorrectGeneLocus(
                hla.chars().take(2).collect::<String>(),
            ));
        }

        let expression_change =
            ExpressionChange::from(s.chars().last().ok_or(HLAError::GeneNameTooShort)?);

        let mut nomenclature_digits = hla.chars().filter(|c| c.is_numeric());

        let second_field_size = match nomenclature_digits.clone().count() % 2 {
            0 => 2usize,
            _ => 3usize,
        };

        Ok(Self {
            gene,
            allele_group: extract(&mut nomenclature_digits, 2),
            hla_protein: string_to_option(extract(&mut nomenclature_digits, second_field_size)),
            cds_synonymous_sub: string_to_option(extract(&mut nomenclature_digits, 2)),
            non_coding_diff: string_to_option(extract(&mut nomenclature_digits, 2)),
            expression_change,
        })
    }
}

fn string_to_option<T: AsRef<str>>(s: T) -> Option<String> {
    Some(s.as_ref().to_string()).filter(|s| !s.is_empty())
}

fn extract<I>(it: I, count: usize) -> String
where
    I: IntoIterator<Item = char>,
{
    it.into_iter().take(count).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn expression_change_from_char() {
        let c = 'C';

        assert_eq!(ExpressionChange::from(c), ExpressionChange::C);
    }

    #[test]
    fn gene_from_char_iter() {
        let chars = "DP".chars();
        let gene: Gene = chars.collect::<Gene>();

        assert_eq!(gene, Gene::DP);
    }

    #[test]
    fn to_hla() {
        let hla = HLA::new("HLA-C*01:101").unwrap();
        let expected = HLA {
            gene: Gene::C,
            allele_group: "01".to_string(),
            hla_protein: Some("101".to_string()),
            cds_synonymous_sub: None,
            non_coding_diff: None,
            expression_change: ExpressionChange::Unknown,
        };

        assert_eq!(hla, expected);
    }
}
