// TODO: Curently motifs are assigned by sorting by allele field (e.g. protein number), Ligand motif, then frequency
use std::collections::{HashMap, HashSet};
use std::str::FromStr;

use crate::error::{HtmlParseError, IoError, NomenclatureError};
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

impl std::fmt::Display for AlleleFreq {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use AlleleFreq::*;

        let s = match self {
            Common => "Common or Well Defined",
            Rare => "Rare",
            Unknown => "",
        };

        write!(f, "{}", s)
    }
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

// NOTE: Could cause errors if the ligand map file is wrong
impl Default for KirLigandMap {
    fn default() -> Self {
        let alleles = HashSet::<ClassI>::new();
        let cache = HashMap::<ClassI, KirLigandInfo>::new();

        Self { alleles, cache }
    }
}

impl KirLigandMap {
    pub fn new() -> Self {
        KirLigandMap::default()
    }

    pub fn updated() -> std::result::Result<Self, HtmlParseError> {
        KirLigandMap::from_loci(&GENE_LOCI)
    }

    pub fn init() -> std::result::Result<Self, IoError> {
        let mut alleles = HashSet::<ClassI>::new();
        let mut cache = HashMap::<ClassI, KirLigandInfo>::new();

        for (row, line) in crate::LIGAND_MAP_DEF.lines().enumerate() {
            let entry = line.split('\t').collect::<Vec<&str>>();

            let allele = entry[0]
                .parse::<ClassI>()
                .or_else(|_| Err(IoError::CouldNotReadAllele(row + 1)))?;

            let motif = entry[1]
                .parse::<LigandMotif>()
                .or_else(|_| Err(IoError::CouldNotReadMotif(row + 1)))?;

            let freq: AlleleFreq = if entry.len() == 3 {
                entry[2].into()
            } else {
                AlleleFreq::Unknown
            };

            let info = KirLigandInfo::new(allele.clone(), motif, freq);

            alleles.insert(allele.clone());
            cache.insert(allele, info);
        }

        Ok(Self { alleles, cache })
    }

    pub fn insert_info(&mut self, info: KirLigandInfo) {
        self.alleles.insert(info.0.clone());
        self.cache.insert(info.0.clone(), info);
    }

    pub fn from_loci(loci: &[&str]) -> std::result::Result<Self, HtmlParseError> {
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

    pub fn from_path<T>(p: T) -> std::result::Result<Self, IoError>
    where
        T: AsRef<std::path::Path>,
    {
        let mut map = KirLigandMap::new();

        let mut rdr = csv::ReaderBuilder::new()
            .has_headers(false)
            .delimiter(b'\t')
            .comment(Some(b'#'))
            .from_path(&p)?;

        for (row, result) in rdr.records().enumerate() {
            let entry = result?;

            let allele = entry[0]
                .parse::<ClassI>()
                .or_else(|_| Err(IoError::CouldNotReadAllele(row + 1)))?;

            let motif = entry[1]
                .parse::<LigandMotif>()
                .or_else(|_| Err(IoError::CouldNotReadMotif(row + 1)))?;

            let freq: AlleleFreq = if entry.len() == 3 {
                entry[2].into()
            } else {
                AlleleFreq::Unknown
            };

            let info = KirLigandInfo::new(allele, motif, freq);
            map.insert_info(info);
        }

        Ok(map)
    }

    // TODO: Look at ebd55784dc5a2d11645c0a754650b705b7407705 for previous messier implementation
    // cannot think of case where more than twice is really necessary but there might be.
    // Decided to for now to only generalise the allele twice. Cases like A02:07,
    // returning A*02:07:01:01 and A*02:07:02 will both have A*02:07:01 and A*02:07:02 present.
    fn get_allele_info(&self, allele: &ClassI) -> Vec<&KirLigandInfo> {
        let mut kir_ligand_info = Vec::<&KirLigandInfo>::new();

        if let Some(allele_info) = self.cache.get(allele) {
            kir_ligand_info.push(allele_info)
        } else {
            self.alleles.iter().for_each(|cached_allele| {
                if let Some(generalised_once) = cached_allele.generalize() {
                    if generalised_once == *allele {
                        kir_ligand_info.push(self.cache.get(cached_allele).unwrap());
                    } else if let Some(generalised_twice) = generalised_once.generalize() {
                        if generalised_twice == *allele {
                            kir_ligand_info.push(self.cache.get(cached_allele).unwrap());
                        }
                    }
                }
            });
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
/// Finds the first HTML table with the possibility of skipping a set of rows
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

        assert_eq!(ligand_info[0], expected);
    }

    #[test]
    fn test_create_ligand_map() {
        let loci = ["A*02:07", "B*57:01", "C*01:02"];
        let ligand_map = KirLigandMap::from_loci(&loci).unwrap();

        let mut motifs = Vec::<LigandMotif>::new();
        let mut expected: Vec<LigandMotif> = vec![
            "Unclassified".parse().unwrap(),
            "C1".parse().unwrap(),
            "Bw4-80I".parse().unwrap(),
        ];

        for (_, v) in ligand_map.cache.into_iter() {
            motifs.push(v.1)
        }

        motifs.sort();
        motifs.dedup();

        expected.sort();
        expected.dedup();

        assert_eq!(motifs, expected);
    }

    // Website has bugs as for example https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?A*02:15
    // A*02:15 should not be resolved to any allele (might loop long!!)
    #[test]
    fn test_query_ligand_map() {
        let loci = ["A*02:15", "A*02:16", "A*02:07"];
        let ligand_map = KirLigandMap::from_loci(&loci).unwrap();

        let query_allele_missing = "A*02:15".parse::<ClassI>().unwrap();
        let query_allele_singly_matched = "A*02:16".parse::<ClassI>().unwrap();
        let query_allele_variable_output = "A*02:07".parse::<ClassI>().unwrap();

        let matching_none = ligand_map.get_allele_info(&query_allele_missing);
        let matching_once = ligand_map.get_allele_info(&query_allele_singly_matched);

        let mut matching_second_option = ligand_map.get_allele_info(&query_allele_variable_output);
        let expected = "A*02:07:01:01".parse::<ClassI>().unwrap();

        assert!(matching_none.is_empty());
        assert_eq!(matching_once.len(), 1);

        matching_second_option.sort();
        assert_eq!(*matching_second_option[0].allele(), expected)
    }

    #[test]
    fn test_map_from_file() {
        let mut map = KirLigandMap::from_path("tests/lg.tsv").unwrap();

        let expected_allele = "A*01:01:01:01".parse::<ClassI>().unwrap();
        let expected_motif = "Unclassified".parse::<LigandMotif>().unwrap();
        let expected_freq = AlleleFreq::Common;

        let expected_ligand_info =
            KirLigandInfo::new(expected_allele.clone(), expected_motif, expected_freq);

        let first_allele = map.alleles.iter().next().unwrap();
        let ligand_info = map.cache.remove(first_allele).unwrap();

        assert_eq!(*first_allele, expected_allele);
        assert_eq!(ligand_info, expected_ligand_info);
    }
}
