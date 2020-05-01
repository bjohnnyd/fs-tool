use std::iter::FromIterator;
use std::str::FromStr;

use crate::error::NomenclatureError;
use crate::ig_like::kir_ligand::LigandMotif;

use log::{debug, error, info, warn};

type Result<T> = std::result::Result<T, NomenclatureError>;

/* Nomenclature */

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
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
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
    type Err = NomenclatureError;

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
            s => Err(NomenclatureError::UnknownExpressionChangeTag(s.to_string())),
        }
    }
}


/* HLA Class I */

#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub struct ClassI {
    pub gene: Gene,
    pub allele_group: String,
    pub hla_protein: Option<String>,
    pub cds_synonymous: Option<String>,
    pub non_coding: Option<String>,
    pub expression_change: ExpressionChange,
    pub ligand_motif: Option<LigandMotif>,
}

impl std::str::FromStr for ClassI {
    type Err = NomenclatureError;

    fn from_str(s: &str) -> Result<Self> {
        let hla = s
            .trim_start_matches("HLA-")
            .replace("*", "");

        if !hla.contains(':') && hla.len() > 3 {
            error!("HLA-I allele {} is longer than 3 characters but does not contain colons which are required", s);
            Err(NomenclatureError::CouldNotParseClassI(s.to_string()))
        } else {

            let mut hla_parts =  hla
                        .split(':');

            let (gene, allele_group) = hla_parts
                .next()
                .map(|required_fields|  {

                    let gene = required_fields[0..1].chars().collect::<Gene>();
                    let allele_group = required_fields[1..].to_string();

                    (gene, allele_group)
                })
                .ok_or_else(|| NomenclatureError::GeneUnknown(s.to_string()))?;

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
                    ligand_motif: None,
                }
            )
        }
    }
}

impl std::fmt::Display for ClassI {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {


        let s  = format!(
            "{}{}{}{}{}{}",
            self.gene,
            self.allele_group,
            self.hla_protein.clone().unwrap_or_else(|| "".to_string()),
            self.cds_synonymous.clone().unwrap_or_else(|| "".to_string()),
            self.non_coding.clone().unwrap_or_else(|| "".to_string()),
            self.expression_change
        );
        write!(f, "{}", s)
    }
}


impl ClassI {
    pub fn to_nomenclature_string(&self) -> String {
        format!(
            "HLA-{}*{}{}{}{}{}",
            self.gene,
            self.allele_group,
            self.hla_protein.clone().map(|protein| format!(":{}", protein)).unwrap_or_else(|| "".to_string()),
            self.cds_synonymous.clone().map(|synonymous| format!(":{}", synonymous)).unwrap_or_else(|| "".to_string()),
            self.non_coding.clone().map(|non_coding| format!(":{}", non_coding)).unwrap_or_else(|| "".to_string()),
            self.expression_change
        )

    }
}

#[cfg(test)]
mod tests {
    use crate::mhc::mhc_I::{ExpressionChange , ClassI};
    use crate::ig_like::kir_ligand::LigandMotif;
    use crate::mhc::hla::Gene;

    #[test]
    fn test_expression_change() {
        let expression_change = "Q".parse().unwrap();

        assert_eq!(ExpressionChange::Q, expression_change);
        assert_eq!(String::from("Q"), expression_change.to_string());
    }

    #[test]
    fn test_known_ligands() {
        let ligand_motif = "A03".parse::<LigandMotif>().unwrap();
        assert_eq!(LigandMotif::A3, ligand_motif)
    }
    #[test]
    fn test_mhcI_parse() {
        let mhc_I = "A03:02:101".parse::<ClassI>().unwrap();
        assert_eq!(
            ClassI {
                gene: Gene::A,
                allele_group: "03".to_string(),
                hla_protein: Some("02".to_string()),
                cds_synonymous: Some("101".to_string()),
                non_coding: None,
                expression_change: ExpressionChange::Unknown,
                ligand_motif: None
            },
            mhc_I);
    }

    #[test]
    fn test_mhcI_into_str() {
        let mhc_I = ClassI {
                gene: Gene::A,
                allele_group: "03".to_string(),
                hla_protein: Some("02".to_string()),
                cds_synonymous: Some("101".to_string()),
                non_coding: None,
                expression_change: ExpressionChange::Unknown,
                ligand_motif: None
            };


        assert_eq!(
            mhc_I.to_string(),
            "A0302101".to_string()
        );

        assert_eq!(
            mhc_I.to_nomenclature_string(),
            "HLA-A*03:02:101".to_string()
        );
    }
}
