use crate::prelude::error::HLAError;
use crate::prelude::fs_tool::LigandInfo;
use crate::prelude::traits::*;

trait ToDisplay {
    type CanDisplay: std::fmt::Display;
    fn to_display(&self) -> Self::CanDisplay;
}

impl<T> ToDisplay for Option<T>
where
    T: std::fmt::Display,
{
    type CanDisplay = String;

    fn to_display(&self) -> String {
        match self {
            Some(val) => format!("{}", val),
            None => String::from(""),
        }
    }
}

impl ToDisplay for ExpressionChange {
    type CanDisplay = &'static str;

    fn to_display(&self) -> &'static str {
        self.into()
    }
}

impl ToDisplay for Gene {
    type CanDisplay = &'static str;

    fn to_display(&self) -> &'static str {
        self.into()
    }
}
impl ToDisplay for HLA {
    type CanDisplay = String;

    fn to_display(&self) -> String {
        format!(
            "{}*{}{}{}{}{}",
            self.gene.to_display(),
            self.allele_group,
            self.hla_protein.to_display(),
            self.cds_synonymous_sub.to_display(),
            self.non_coding_diff.to_display(),
            self.expression_change.to_display()
        )
    }
}

impl std::fmt::Display for HLA {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.to_display())
    }
}

type HLAResult<T> = std::result::Result<T, HLAError>;
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct HLA {
    pub gene: Gene,
    pub allele_group: String,
    pub hla_protein: Option<String>,
    pub cds_synonymous_sub: Option<String>,
    pub non_coding_diff: Option<String>,
    pub expression_change: ExpressionChange,
    pub ligand_group: Option<LigandGroup>,
    pub ipd_frequency: Option<IPDFrequency>,
}

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

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum LigandGroup {
    A11,
    A3,
    Bw4_80T,
    Bw4_80I,
    C1,
    C2,
    Unclassified,
}

impl FromStr for LigandGroup {
    type Err = HLAError;

    fn from_str(s: &str) -> HLAResult<Self> {
        let s = s
            .chars()
            .filter(|c| !c.is_ascii_whitespace())
            .collect::<String>();
        match s.as_str() {
            "A11" => Ok(LigandGroup::A11),
            "A3" => Ok(LigandGroup::A3),
            "Bw4-80T" => Ok(LigandGroup::Bw4_80T),
            "Bw4-80I" => Ok(LigandGroup::Bw4_80I),
            "C1" => Ok(LigandGroup::C1),
            "C2" => Ok(LigandGroup::C2),
            "Unclassified" => Ok(LigandGroup::Unclassified),
            _ => Err(HLAError::UnknownLigandGroup(s.to_string())),
        }
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum IPDFrequency {
    Rare,
    Common,
    Unknown,
}

impl FromStr for IPDFrequency {
    type Err = HLAError;

    fn from_str(s: &str) -> HLAResult<Self> {
        match s.trim() {
            "Common or Well Defined" => Ok(IPDFrequency::Common),
            "Rare" => Ok(IPDFrequency::Rare),
            "Unknown" => Ok(IPDFrequency::Unknown),
            _ => Err(HLAError::IPDFrequencyIncorrect(s.to_string())),
        }
    }
}

impl Gene {
    fn is_unknown(&self) -> bool {
        self.eq(&Gene::Unknown)
    }
}
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
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

impl From<&ExpressionChange> for &str {
    fn from(e: &ExpressionChange) -> &'static str {
        match e {
            ExpressionChange::N => "N",
            ExpressionChange::L => "L",
            ExpressionChange::S => "S",
            ExpressionChange::C => "C",
            ExpressionChange::A => "A",
            ExpressionChange::Q => "Q",
            ExpressionChange::Unknown => "",
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
    pub fn new<T: AsRef<str>>(s: T) -> HLAResult<Self> {
        let hla_name = s.as_ref().trim_start_matches("HLA-");
        hla_name.parse::<HLA>()
    }

    pub fn set_ligand_group(&mut self, lg: LigandGroup) {
        self.ligand_group = Some(lg)
    }
}

macro_rules! to_option {
    ($s:ident) => {{
        Some($s).filter(|s| !s.is_empty())
    }};
}

impl std::str::FromStr for HLA {
    type Err = HLAError;

    fn from_str(s: &str) -> HLAResult<Self> {
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
            ligand_group: None,
            ipd_frequency: None,
        })
    }
}

impl TryFrom<LigandInfo> for HLA {
    type Error = HLAError;

    fn try_from(ligand_info: LigandInfo) -> HLAResult<Self> {
        let mut hla = ligand_info.0.parse::<HLA>()?;
        let lg = ligand_info.1.parse::<LigandGroup>()?;
        let freq = ligand_info.2.parse::<IPDFrequency>()?;
        hla.ligand_group = Some(lg);
        hla.ipd_frequency = Some(freq);
        Ok(hla)
    }
}

#[inline]
fn string_to_option<T: AsRef<str>>(s: T) -> Option<String> {
    Some(s.as_ref().to_string()).filter(|s| !s.is_empty())
}

#[inline]
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
            ligand_group: None,
            ipd_frequency: None,
        };

        assert_eq!(hla, expected);
    }
    #[test]
    fn to_hla_incomplete() {
        let hla = HLA::new("A01").unwrap();
        let expected = HLA {
            gene: Gene::A,
            allele_group: "01".to_string(),
            hla_protein: None,
            cds_synonymous_sub: None,
            non_coding_diff: None,
            expression_change: ExpressionChange::Unknown,
            ligand_group: None,
            ipd_frequency: None,
        };

        assert_eq!(hla, expected);
    }

    #[test]
    fn test_display_hla() {
        let hla1 = HLA::new("A01").unwrap();
        let hla2 = HLA::new("A01:102").unwrap();
        println!("{}{}", hla1, hla2);
    }
}
