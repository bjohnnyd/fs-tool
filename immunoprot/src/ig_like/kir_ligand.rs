use std::str::FromStr;
use std::collections::HashSet;

use crate::error::NomenclatureError;
use crate::mhc::hla::ClassI;

type Result<T> = std::result::Result<T, NomenclatureError>;

const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";
const GENE_LOCI: [&str; 3] = ["A", "B", "C"];


#[derive(Debug, Eq, PartialEq)]
pub struct LigandInfo(pub ClassI, pub LigandMotif, pub AlleleFreq);

#[derive(Debug, Eq, PartialEq)]
pub enum AlleleFreq {
    Rare,
    Common,
    Unknown
}

impl<T> From<T> for AlleleFreq
where T: AsRef<str>  {
    fn from(s: T) -> Self {
        use AlleleFreq::*;
        match s.as_ref().trim() {
            _match if _match.contains("Common") => Common,
            "Rare" => Rare,
            _ => Unknown
        }
    }
}

// TODO: Create struct for IPD information storage that implements can be create from vec of strings,
// written to file etc.
// TODO: Need to create appropriate errors and split function into one that gets HTML and one that
// parses the tables and creates IPD info struct
mod reader {
    use scraper::{Html, Selector};
    use crate::ig_like::kir_ligand::{LigandInfo, AlleleFreq, LigandMotif};
    use crate::mhc::hla::ClassI;

    pub fn get_ipd_html<T>(gene_locus: T) -> Html
        where T: AsRef<str> + std::fmt::Display {
        let resp = attohttpc::get(format!("https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?{}", gene_locus));
        let text = resp.send().unwrap().text().unwrap();

        Html::parse_document(&text)
    }

    pub fn read_table(html: &Html, skip_rows: usize) -> Vec<LigandInfo> {
        let mut result = Vec::<LigandInfo>::new();

        let selector = Selector::parse("tr").unwrap();

        for row in html.select(&selector).skip(skip_rows) {
            let table_row = row.text().collect::<Vec<&str>>();
            let ligand_info = LigandInfo(
                table_row[0].parse::<ClassI>().unwrap(),
                table_row[1].parse::<LigandMotif>().unwrap(),
                table_row[2].into()
            );

            result.push(ligand_info);
        }

        result
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
    type Err = NomenclatureError;

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
            s => Err(NomenclatureError::UnknownLigandMotif(s.to_string())),
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
    use crate::mhc::hla::{ExpressionChange, ClassI};
    use crate::ig_like::kir_ligand::{LigandMotif, reader::get_ipd_html, reader::read_table};

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
        let html = get_ipd_html("C*01:02");
        let ligand_info = read_table(&html, 1);
        dbg!(ligand_info);
    }
}
