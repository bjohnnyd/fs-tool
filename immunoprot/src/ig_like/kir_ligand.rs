// TODO: Curently motifs are assigned by sorting by allele field (e.g. protein number), Ligand motif, then frequency
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::error::{HtmlParseError, NomenclatureError};
use crate::mhc::hla::ClassI;
use log::info;
use scraper::{Html, Selector};

type Result<T> = std::result::Result<T, NomenclatureError>;

pub const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";
pub const GENE_LOCI: [&str; 3] = ["A", "B", "C"];
pub const SKIP_ROWS: usize = 1;

#[derive(Debug, Eq, PartialEq, Hash, Clone, Ord, PartialOrd)]
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
            Unclassified => "Unclassified",
        };
        write!(f, "{}", motif)
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone, Ord, PartialOrd)]
pub enum AlleleFreq {
    Common,
    Rare,
    Unknown,
}

impl<T> From<T> for AlleleFreq
where
    T: AsRef<str>,
{
    fn from(s: T) -> Self {
        use AlleleFreq::*;
        match s.as_ref().trim() {
            _match if _match.contains("Common") => Common,
            "Rare" => Rare,
            _ => Unknown,
        }
    }
}

#[derive(Debug, Eq, PartialEq, Hash, Clone, PartialOrd, Ord)]
pub struct KirLigandInfo(ClassI, LigandMotif, AlleleFreq);

impl KirLigandInfo {
    pub fn new(hla: ClassI, motif: LigandMotif, freq: AlleleFreq) -> Self {
        Self(hla, motif, freq)
    }

    pub fn allele(&self) -> &ClassI {
        &self.0
    }

    pub fn motif(&self) -> &LigandMotif {
        &self.1
    }

    pub fn freq(&self) -> &AlleleFreq {
        &self.2
    }
}

#[derive(Debug, Eq, PartialEq)]
pub struct KirLigandMap {
    pub alleles: HashSet<ClassI>,
    pub cache: HashMap<ClassI, KirLigandInfo>,
}

impl KirLigandMap {
    fn new(loci: &[&str]) -> std::result::Result<Self, HtmlParseError> {
        let mut alleles = HashSet::<ClassI>::new();
        let mut cache = HashMap::<ClassI, KirLigandInfo>::new();

        let results: std::result::Result<Vec<_>, HtmlParseError> = loci
            .iter()
            .map(|locus| {
                let raw_html = get_ipd_html(locus)?;
                let allele_infos = read_table(&raw_html, SKIP_ROWS)?;

                for allele_info in allele_infos {
                    alleles.insert(allele_info.0.clone());
                    cache.insert(allele_info.0.clone(), allele_info);
                }

                Ok(())
            })
            .collect();

        results?;

        Ok(Self { alleles, cache })
    }

    // TODO: Messy and inefficient as it iterates over all backwards removing fields
    // cannot think of case where more than once is really necessary but there might be A*02 returning A*02:01:01
    // might be better to return multiple removed at once so cases like A02:07 returning A*02:07:01:01 and A*02:07:02
    // will both have A*02:07:01 and A*02:07:02 present?
    fn get_allele_info(&self, allele: &ClassI) -> Vec<&KirLigandInfo> {
        let mut kir_ligand_info = Vec::<&KirLigandInfo>::new();

        if let Some(allele_info) = self.cache.get(allele) {
            kir_ligand_info.push(allele_info)
        } else {
            let mut present_alleles = self
                .alleles
                .iter()
                .map(|allele_with_info| (allele_with_info.clone(), allele_with_info.generalize()))
                .filter(|(_, general)| general.is_some())
                .collect::<HashMap<ClassI, Option<ClassI>>>();

            let mut no_more_alleles = true;

            'find_alleles: loop {
                for (original_allele, generalized_allele) in present_alleles.iter() {
                    if let Some(generalized_allele) = generalized_allele {
                        no_more_alleles = false;
                        if generalized_allele == allele {
                            kir_ligand_info.push(self.cache.get(original_allele).unwrap());
                        }
                    }
                }
                if !kir_ligand_info.is_empty() || no_more_alleles {
                    break 'find_alleles;
                } else {
                    no_more_alleles = true;
                    present_alleles = present_alleles
                        .into_iter()
                        .map(|(original, general)| {
                            if let Some(general) = general {
                                let general = general.generalize();
                                (original, general)
                            } else {
                                (original, None)
                            }
                        })
                        .filter(|(_, general)| general.is_some())
                        .collect::<HashMap<ClassI, Option<ClassI>>>();
                }
            }
        }
        kir_ligand_info
    }
}

/// Obtains raw HTL from the EBI website
pub fn get_ipd_html<T>(gene_locus: T) -> std::result::Result<Html, HtmlParseError>
where
    T: AsRef<str> + std::fmt::Display,
{
    let url = format!("{}{}", IPD_KIR_URL, &gene_locus);
    info!("Connecting to {}...", &url);
    let request = attohttpc::get(&url);
    let response = request.send()?;
    let text = response
        .text()
        .or_else(|err| Err(HtmlParseError::CouldNotReadResponse(err)))?;

    info!(
        "Obtained response, looking for allele table for locus '{}'",
        gene_locus
    );

    Ok(Html::parse_document(&text))
}

// TODO: Need to deal with cases where supplied HLA allele is incorrect format
/// Find the first HTML table and can skip a desired set of rows
pub fn read_table(
    html: &Html,
    skip_rows: usize,
) -> std::result::Result<Vec<KirLigandInfo>, HtmlParseError> {
    let mut result = Vec::<KirLigandInfo>::new();

    // TODO: Deal with error
    let selector = Selector::parse("tr").unwrap();
    info!("Found HLA allele table! Reading rows...");

    for row in html.select(&selector).skip(skip_rows) {
        let table_row = row.text().collect::<Vec<&str>>();

        if table_row.len() == 3 || table_row.len() == 2 {
            let allele = table_row[0]
                .parse::<ClassI>()
                .or_else(|_| Err(HtmlParseError::CouldNotReadClassI(table_row[0].to_string())))?;
            let motif = table_row[1]
                .parse::<LigandMotif>()
                .or_else(|_| Err(HtmlParseError::CouldNotReadClassI(table_row[1].to_string())))?;
            let freq: AlleleFreq = if table_row.len() == 3 {
                table_row[2].into()
            } else {
                AlleleFreq::Unknown
            };

            let ligand_info = KirLigandInfo::new(allele, motif, freq);
            result.push(ligand_info);
        } else {
            return Err(HtmlParseError::IncorrectNumberOfColumns(
                table_row.len(),
                table_row.join(""),
            ));
        }
    }
    info!("Finished reading table");
    Ok(result)
}

#[cfg(test)]
mod tests {
    use crate::ig_like::kir_ligand::{
        get_ipd_html, read_table, AlleleFreq, KirLigandInfo, KirLigandMap, LigandMotif,
    };
    use crate::mhc::hla::ClassI;

    #[test]
    fn test_known_ligands() {
        let lg_group = "A03".parse::<LigandMotif>().unwrap();
        assert_eq!(LigandMotif::A3, lg_group)
    }

    #[test]
    fn test_ligand_info() {
        let lg_info = include_str!("../resources/2019-12-29_lg.tsv");

        lg_info.lines().for_each(|l| {
            dbg!(l);
        });
    }

    #[test]
    fn test_connect_to_ipd() {
        let html = get_ipd_html("C*01:02").unwrap();
        let ligand_info = read_table(&html, 1).unwrap();
        let expected = KirLigandInfo::new(
            "C01:02:01:01".parse::<ClassI>().unwrap(),
            LigandMotif::C1,
            AlleleFreq::Common,
        );

        assert_eq!(expected, ligand_info[0]);
    }

    #[test]
    fn test_create_ligand_map() {
        let loci = ["A*02:07", "B*57:01", "C0*01:02"];
        let ligand_map = KirLigandMap::new(&loci).unwrap();

        let mut motifs = Vec::<LigandMotif>::new();
        let mut expected: Vec<LigandMotif> = vec![
            "Unclassified".parse().unwrap(),
            "Bw4-80I".parse().unwrap(),
            "Bw4-80I".parse().unwrap(),
        ];

        for (k, v) in ligand_map.cache.into_iter() {
            motifs.push(v.1)
        }

        motifs.sort();
        motifs.dedup();

        expected.sort();
        expected.dedup();

        assert_eq!(expected, motifs);
    }

    // Website has bugs as for example https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?A*02:15
    // A*02:15 should not be resolved to any allele (might loop long!!)
    #[test]
    fn test_query_ligand_map() {
        stderrlog::new()
            .module("immunoprot")
            .verbosity(3)
            .init()
            .unwrap();
        let loci = ["A*02:15", "A*02:16", "A*02:07"];
        let ligand_map = KirLigandMap::new(&loci).unwrap();

        let query_allele_missing = "A*02:15".parse::<ClassI>().unwrap();
        let query_allele_singly_matched = "A*02:16".parse::<ClassI>().unwrap();
        let query_allele_variable_output = "A*02:07".parse::<ClassI>().unwrap();

        let mut matching_none = ligand_map.get_allele_info(&query_allele_missing);

        let mut matching_once = ligand_map.get_allele_info(&query_allele_singly_matched);

        let mut matching_second_option = ligand_map.get_allele_info(&query_allele_variable_output);

        assert!(matching_none.is_empty());
        assert_eq!(1, matching_once.len())
    }
}
