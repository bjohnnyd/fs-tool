// nomenclature URL: https://www.ebi.ac.uk/ipd/kir/alleles.html
// allele info example URL: https://www.ebi.ac.uk/cgi-bin/ipd/kir/allele.cgi?KIR2DL1*

use std::fmt::Formatter;
use std::str::FromStr;

use crate::error::NomenclatureError;

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum Tail {
    Long,  // Inhibitory
    Short, // Activating
    Pseudo,
}

impl std::fmt::Display for Tail {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use Tail::*;
        let s = match self {
            Long => "L",
            Short => "S",
            Pseudo => "P",
        };
        write!(f, "{}", s)
    }
}

impl FromStr for Tail {
    type Err = NomenclatureError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Tail::*;
        match s {
            "L" => Ok(Long),
            "S" => Ok(Short),
            "P" => Ok(Pseudo),
            s => Err(NomenclatureError::IncorrectKirTail(s.to_string())),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub enum Domain {
    Two,
    Three,
}

impl std::fmt::Display for Domain {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use Domain::*;
        let s = match self {
            Two => "2D",
            Three => "3D",
        };
        write!(f, "{}", s)
    }
}

impl FromStr for Domain {
    type Err = NomenclatureError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use Domain::*;
        match s {
            "2D" => Ok(Two),
            "3D" => Ok(Three),
            s => Err(NomenclatureError::IncorrectKirDomain(s.to_string())),
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub enum KirProtein {
    KP1,
    KP2,
    KP3,
    KP4,
    KP5,
    KP5A,
    KP5B,
}

impl std::fmt::Display for KirProtein {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        use KirProtein::*;
        let s = match self {
            KP1 => "1",
            KP2 => "2",
            KP3 => "3",
            KP4 => "4",
            KP5 => "5",
            KP5A => "5A",
            KP5B => "5B",
        };
        write!(f, "{}", s)
    }
}

impl FromStr for KirProtein {
    type Err = NomenclatureError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        use KirProtein::*;
        match s {
            "1" => Ok(KP1),
            "2" => Ok(KP2),
            "3" => Ok(KP3),
            "4" => Ok(KP4),
            "5" => Ok(KP5),
            "5A" => Ok(KP5A),
            "5B" => Ok(KP5B),
            s => Err(NomenclatureError::IncorrectKirProtein(s.to_string())),
        }
    }
}
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct KirAllele {
    series: String,
    cds_syn_sub: Option<String>,
    non_coding_sub: Option<String>,
}

impl std::fmt::Display for KirAllele {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}",
            self.series,
            self.cds_syn_sub.clone().unwrap_or_else(|| "".to_string()),
            self.non_coding_sub
                .clone()
                .unwrap_or_else(|| "".to_string()),
        )
    }
}

impl FromStr for KirAllele {
    type Err = NomenclatureError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let series = s.chars().take(3).collect::<String>();
        let cds_syn_sub = s.chars().skip(3).take(2).collect::<String>();
        let non_coding_sub = s.chars().skip(5).take(2).collect::<String>();

        if series.len() != 3 {
            Err(NomenclatureError::IncorrectKirAllele(s.to_string()))
        } else {
            Ok(Self {
                series,
                cds_syn_sub: if cds_syn_sub.is_empty() {
                    None
                } else {
                    Some(cds_syn_sub)
                },
                non_coding_sub: if non_coding_sub.is_empty() {
                    None
                } else {
                    Some(non_coding_sub)
                },
            })
        }
    }
}

#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Kir {
    ig_like_domain: Domain,
    cytoplasmic_tail: Tail,
    protein: KirProtein,
    allele: Option<KirAllele>,
}

#[cfg(test)]
mod tests {
    use crate::ig_like::kir::Domain::{self, *};
    use crate::ig_like::kir::Tail::{self, *};

    #[test]
    fn test_kir_tail_naming_correct() {
        let correct = "L";
        let incorrect = "Z";

        assert_eq!(correct.parse::<Tail>().unwrap(), Long);
    }
    #[test]
    #[should_panic]
    fn test_kir_tail_naming_incorrect() {
        let incorrect = "Z";
        incorrect.parse::<Tail>().unwrap();
    }
}
