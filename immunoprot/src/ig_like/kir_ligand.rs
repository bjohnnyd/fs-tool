use std::str::FromStr;
use std::collections::HashSet;

use crate::error::HLAError;
use crate::mhc::mhc_I::MHCI;
use scraper::{Html, Selector};

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

// TODO: Create struct for IPD information storage that implements can be create from vec of strings,
// written to file etc.
// TODO: Need to create appropriate errors and split function into one that gets HTML and one that
// parses the tables and creates IPD info struct
pub fn connect_to_ipd<T>(gene_locus: T)
where T: AsRef<str> + std::fmt::Display {
    let resp = attohttpc::get(format!("https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?{}", gene_locus));
    let text = resp.send().unwrap().text().unwrap();

    let page = Html::parse_document(&text);
    let selector = Selector::parse("tr").unwrap();

    for row in page.select(&selector).skip(1) {
        let row_text = row.text().collect::<Vec<&str>>();
        let ligand_motif: LigandMotif = row_text[1].parse().unwrap();
    }
}

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

impl std::fmt::Display for LigandMotif {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use LigandMotif::*;
        let motif = match self {
            A3 => "A3",
            A11 => "A11",
            Bw4_80I => "Bw4-80I",
            Bw4_80T => "Bw4-80T",
            Bw6 => "Bw6",
            C1 => "C1",
            C2 => "C2",
            Unclassified => "Unclassified"
        };
        write!(f, "{}", motif)
    }
}


#[cfg(test)]
mod tests {
    use crate::mhc::mhc_I::{ExpressionChange, MHCI};
    use crate::ig_like::kir_ligand::{LigandMotif, connect_to_ipd};

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

    #[test]
    fn test_connect_to_ipd() {
        // let lg_info = include_str!("../resources/2019-12-29_lg.tsv");
        connect_to_ipd("C*01:02");
    }
}
