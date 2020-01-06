use crate::error::*;
use crate::prelude::external::{Html, Selector};
use crate::prelude::io::*;
use rayon::prelude::*;

const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";
const GENE_LOCI: [&str; 3] = ["A", "B", "C"];

type Result<T, E = Error> = std::result::Result<T, E>;

/// Converts `cssparser::ParseError` to a RetrieveLigandError
macro_rules! selector_error_convert {
    ($e:expr) => {{
        match Selector::parse($e) {
            Ok(selector) => Ok(selector),
            Err(err) => {
                let line = err.location.line + 1;
                let column = err.location.column;
                Err(Error::WebsiteParsingError { line, column })
            }
        }
    }};
}

#[derive(Debug, Eq, PartialEq)]
pub struct LigandInfo(pub String, pub String, pub String);

impl From<Vec<&str>> for LigandInfo {
    fn from(info: Vec<&str>) -> Self {
        let mut frequency = String::new();
        if info.len() == 2 {
            frequency.push_str("Unknown");
        } else {
            frequency = info[2].to_string();
        }
        LigandInfo(info[0].to_string(), info[1].to_string(), frequency)
    }
}

impl From<&str> for LigandInfo {
    fn from(s: &str) -> Self {
        LigandInfo::from(s.split_ascii_whitespace().collect::<Vec<&str>>())
    }
}

impl std::fmt::Display for LigandInfo {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        writeln!(
            f,
            "HLA: {}, LigandGroup: {}, Frequency: {}",
            self.0, self.1, self.2
        )
    }
}

fn is_ipd_search_safe<T: AsRef<str>>(hla: T) -> bool {
    let query = hla.as_ref();

    if query.len() == 1 {
        match query {
            "A" | "B" | "C" => true,
            _ => false,
        }
    } else {
        query.contains('*') && query.len() <= 4 || query.contains('*') && query.contains(':')
    }
}

//fn clean_hla<T: AsRef<str>>(hla: &T) -> Result<&str, RetrieveLigandError> {
fn clean_hla<T: AsRef<str>>(hla: &T) -> Result<&str> {
    let non_prefix_hla = hla.as_ref().trim_start_matches("HLA-");

    if is_ipd_search_safe(non_prefix_hla) {
        Ok(non_prefix_hla)
    } else {
        Err(Error::IncorrectHLA {
            hla: hla.as_ref().to_string(),
        })
    }
}

fn parse_ipd_response(response_html: String) -> Result<Vec<LigandInfo>> {
    let mut result = Vec::<LigandInfo>::new();
    let table_selector = selector_error_convert!("table")?;
    let row_selector = selector_error_convert!("tr")?;
    let document = Html::parse_document(&response_html);

    if let Some(table) = document.select(&table_selector).next() {
        let rows = table.select(&row_selector).skip(1);

        for row in rows {
            let hla_info: LigandInfo = row.text().collect::<Vec<&str>>().into();
            result.push(hla_info);
        }
        Ok(result)
    } else {
        let response = response_html.to_string();
        Err(Error::NoLigandTableFound { response })
    }
}

fn get_ipd_website(url: &str) -> Result<String> {
    let res = attohttpc::get(url).send().context(WebsiteAccessError {
        url: url.to_string(),
    })?;
    Ok(res.text().context(CouldNotReadResponse {
        url: url.to_string(),
    })?)
}

pub fn retrieve_ligand_group<T>(hla: &T) -> Result<Vec<LigandInfo>>
where
    T: AsRef<str>,
{
    let url = format!("{}{}", IPD_KIR_URL, &clean_hla(&hla)?);
    let resp = get_ipd_website(url.as_ref())?;
    Ok(parse_ipd_response(resp)?)
}

fn get_ligand_db_file() -> Result<PathBuf> {
    if let Some(proj_dirs) = directories::ProjectDirs::from("", "", "fs-tool") {
        let out_dir = proj_dirs.data_local_dir();
        if !out_dir.exists() {
            fs::create_dir_all(&out_dir).context(CouldNotCreateDirectory {
                out_dir: out_dir.to_string_lossy().to_string(),
            });
        }
        Ok(out_dir.join("ligand_groups.tsv"))
    } else {
        Err(Error::NoLocalData {})
    }
}

fn save_ligand_groups(p: &PathBuf) -> Result<()> {
    let ligands = GENE_LOCI
        .par_iter()
        .filter_map(|gene| retrieve_ligand_group(gene).ok())
        .flatten()
        .collect::<Vec<LigandInfo>>();

    if ligands.is_empty() {
        return Err(Error::NoLigandInformation {
            url: IPD_KIR_URL.to_string(),
        });
    }
    {
        let mut f = File::create(p).context(CouldNotOpenFile {
            f_path: p.clone().to_owned(),
        })?;
        for ligand_info in &ligands {
            f.write_all(
                format!("{}\t{}\t{}\n", ligand_info.0, ligand_info.1, ligand_info.2).as_ref(),
            )
            .context(WriteLigandInformationToFile {
                p: p.clone().to_owned(),
                info_1: ligand_info.0.clone(),
                info_2: ligand_info.1.clone(),
                info_3: ligand_info.2.clone(),
            });
        }
    }
    Ok(())
}

pub fn get_ligand_table(update_ligand_groups: bool) -> Result<impl AsRef<str>> {
    let mut ligand_table = String::new();
    let ligand_db_file = get_ligand_db_file()?;
    if update_ligand_groups || !ligand_db_file.exists() {
        save_ligand_groups(&ligand_db_file)?;
    }
    let mut f = File::open(&ligand_db_file).context(CouldNotOpenFile {
        f_path: ligand_db_file.clone(),
    })?;
    f.read_to_string(&mut ligand_table)
        .context(CouldNotReadLigandFile {
            ligand_db_file: ligand_db_file.clone(),
        })?;

    Ok(ligand_table)
}

pub fn parse_ligand_table<T: AsRef<str>>(text: T) -> Vec<LigandInfo> {
    text.as_ref().lines().map(LigandInfo::from).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_valid_hla() {
        let hla = "HLA-A*01";
        let cleaned_hla = clean_hla(&hla).unwrap();

        assert_eq!(cleaned_hla, "A*01");
    }

    #[test]
    fn test_reqwest() {
        let hla = "C*01:102";
        let result = retrieve_ligand_group(&hla).unwrap();

        assert_eq!(
            result[0],
            LigandInfo(
                "C*01:102".to_string(),
                "C1".to_string(),
                "Unknown".to_string()
            )
        );
    }
}
