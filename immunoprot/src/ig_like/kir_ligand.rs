use std::str::FromStr;
use crate::error::HLAError;
use std::collections::HashSet;

type Result<T> = std::result::Result<T, HLAError>;

const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";
const GENE_LOCI: [&str; 3] = ["A", "B", "C"];


// #[derive(Debug, Eq, PartialEq, Hash, Clone)]
// pub struct KirLigandMap {
//     alleles: HashSet,
//     :
//
//
// }
#[derive(Debug, Eq, PartialEq, Hash, Clone)]
pub enum LigandMotif {
    A11,
    A3,
    Bw4_80T,
    Bw4_80I,
    Bw6,
    C1,
    C2,
    Unclassified,
}

impl FromStr for LigandMotif {
    type Err = HLAError;

    fn from_str(s: &str) -> Result<Self> {
        let s = s.chars().filter(|c| !c.is_whitespace()).collect::<String>();

        match s.as_str() {
            "A11" => Ok(LigandMotif::A11),
            "A3" | "A03" => Ok(LigandMotif::A3),
            "Bw4-80T" => Ok(LigandMotif::Bw4_80T),
            "Bw4-80I" => Ok(LigandMotif::Bw4_80I),
            "Bw6" => Ok(LigandMotif::Bw6),
            "C1" | "C01" => Ok(LigandMotif::C1),
            "C2" | "C02" => Ok(LigandMotif::C2),
            "Unclassified" => Ok(LigandMotif::Unclassified),
            s => Err(HLAError::UnknownLigandMotif(s.to_string())),
        }
    }
}

// pub struct IPDFrequency {
//     Rare,
//
// }

#[cfg(test)]
mod tests {
    use crate::mhc::mhc_I::{ExpressionChange, MHCI};
    use crate::ig_like::kir_ligand::LigandMotif;

    #[test]
    fn test_known_ligands() {
        let lg_group = "A03".parse::<LigandMotif>().unwrap();
        assert_eq!(LigandMotif::A3, lg_group)
    }
    #[test]
    fn test_ligand_info() {
        let lg_info = include_str!("../resources/2019-12-29_lg.tsv");

        lg_info
            .lines()
            .for_each(|l|{
                dbg!(l);
            });
    }
}
