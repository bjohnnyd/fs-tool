use log::{debug, error, info, warn};

use crate::error::HLAError;
use crate::mhc::hla::Gene;
use std::fmt::Formatter;
use std::str::FromStr;
use std::process::id;

type Result<T> = std::result::Result<T, HLAError>;

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
        use ExpressionChange::*;
        match c {
            'N' => N,
            'L' => L,
            'S' => S,
            'C' => C,
            'A' => A,
            'Q' => Q,
            _ => Unknown,
        }
    }
}

impl std::fmt::Display for ExpressionChange {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use ExpressionChange::*;

        let s = match self {
            N => "N",
            L => "L",
            S => "S",
            C => "C",
            A => "A",
            Q => "Q",
            Unknown => "",
        }
        .to_string();

        write!(f, "{}", s)
    }
}

impl FromStr for ExpressionChange {
    type Err = HLAError;

    fn from_str(s: &str) -> Result<Self> {
        use ExpressionChange::*;
        match s {
            "N" => Ok(N),
            "L" => Ok(L),
            "S" => Ok(S),
            "C" => Ok(C),
            "A" => Ok(A),
            "Q" => Ok(Q),
            "" => Ok(Unknown),
            s => Err(HLAError::UnknownExpressionChangeTag(s.to_string())),
        }
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum LigandGroup {
    A11,
    A3,
    Bw4_80T,
    Bw4_80I,
    Bw6,
    C1,
    C2,
    Unclassified,
}

impl FromStr for LigandGroup {
    type Err = HLAError;

    fn from_str(s: &str) -> Result<Self> {
        let s = s.chars().filter(|c| !c.is_whitespace()).collect::<String>();

        match s.as_str() {
            "A11" => Ok(LigandGroup::A11),
            "A3" | "A03" => Ok(LigandGroup::A3),
            "Bw4-80T" => Ok(LigandGroup::Bw4_80T),
            "Bw4-80I" => Ok(LigandGroup::Bw4_80I),
            "Bw6" => Ok(LigandGroup::Bw6),
            "C1" | "C01" => Ok(LigandGroup::C1),
            "C2" | "C02" => Ok(LigandGroup::C2),
            "Unclassified" => Ok(LigandGroup::Unclassified),
            s => Err(HLAError::UnknownLigandGroup(s.to_string())),
        }
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct MHCI {
    pub gene: Gene,
    pub allele_group: String,
    pub hla_protein: Option<String>,
    pub cds_synonymous: Option<String>,
    pub non_coding: Option<String>,
    pub expression_change: ExpressionChange,
    pub ligand_group: Option<LigandGroup>,
}

impl std::str::FromStr for MHCI {
    type Err = HLAError;

    fn from_str(s: &str) -> Result<Self> {
        let hla = s
            .trim_start_matches("HLA-")
            .replace("*", "");

        if !hla.contains(':') && hla.len() > 3 {
            error!("MHCI allele {} is longer than 3 characters but does not contain colons which are required", s);
            Err(HLAError::CouldNotParseMHCI(s.to_string()))
        } else {

            let mut hla_parts =  hla
                        .split(':');

            let gene = hla_parts
                .next()
                .map(|gene| gene.chars().collect::<Gene>())
                .ok_or_else(|| HLAError::GeneUnknown(s.to_string()))?;

            let allele_group = hla_parts
                .next()
                .ok_or_else(|| HLAError::NoAlleleGroup(s.to_string()))?
                .to_string();

            let expression_change = hla
                .chars()
                .last()
                .map(ExpressionChange::from)
                .unwrap_or_else(|| ExpressionChange::Unknown);

            Ok(
                Self {
                    gene,
                    allele_group,
                    hla_protein: hla_parts.next().map(String::from),
                    cds_synonymous: hla_parts.next().map(String::from),
                    non_coding: hla_parts.next().map(String::from),
                    expression_change,
                    ligand_group: None,
                }
            )

        }
    }
}

#[cfg(test)]
mod tests {
    use crate::mhc::mhc_I::{ExpressionChange, LigandGroup, MHCI};

    #[test]
    fn test_expression_change() {
        let expression_change = "Q".parse().unwrap();

        assert_eq!(ExpressionChange::Q, expression_change);
        assert_eq!(String::from("Q"), expression_change.to_string());
    }

    #[test]
    fn test_known_ligands() {
        let lg_group = "A03".parse::<LigandGroup>().unwrap();
        assert_eq!(LigandGroup::A3, lg_group)
    }
    #[test]
    fn test_mhcI_parse() {
        let mhc_I = "A03:02".parse::<MHCI>().unwrap();
        dbg!(mhc_I);
        let mhc_I = "A03:02:101".parse::<MHCI>().unwrap();
        dbg!(mhc_I);
    }
}
