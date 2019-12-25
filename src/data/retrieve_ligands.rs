use crate::data::errors::RetrieveLigandError::{self, InvalidHLA};
use scraper::{Html, Selector};
use std::fmt::{Error, Formatter};
use std::path::Path;
//use tokio::fs::File;
use tokio::prelude::*;
use std::fs::File;
use std::io::Write;

const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";

/// Converts `cssparser::ParseError` to a RetrieveLigandError
macro_rules! selector_error_convert {
    ($e:expr) => {{
        match Selector::parse($e) {
            Ok(selector) => Ok(selector),
            Err(err) => {
                let line = err.location.line + 1;
                Err(RetrieveLigandError::CSSParseError(
                    line,
                    err.location.column,
                ))
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

impl std::fmt::Display for LigandInfo {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        writeln!(
            f,
            "HLA: {}, LigandGroup: {}, Frequency: {}",
            self.0, self.1, self.2
        )
    }
}

fn is_ipd_search_safe<T>(hla: T) -> bool
where
    T: AsRef<str>,
{
    let query = hla.as_ref();

    if query.len() == 1 {
        match query {
            "A" | "B" | "C" => true,
            _ => false,
        }
    } else {
        query.contains('*') && query.len() <= 3 || query.contains('*') && query.contains(':')
    }
}
fn clean_hla<T>(hla: &T) -> Result<&str, RetrieveLigandError>
where
    T: AsRef<str>,
{
    let non_prefix_hla = hla.as_ref().trim_start_matches("HLA-");

    if is_ipd_search_safe(non_prefix_hla) {
        Ok(non_prefix_hla)
    } else {
        Err(InvalidHLA(hla.as_ref().to_string()))
    }
}

pub async fn retrieve_ligand_group<T>(hla: &T) -> Result<Vec<LigandInfo>, RetrieveLigandError>
where
    T: AsRef<str>,
{
    let mut result = Vec::<LigandInfo>::new();
    let url = format!("{}{}", IPD_KIR_URL, &clean_hla(&hla)?);
//    let mut resp = reqwest::get(&url).await?;
    let mut resp = surf::get(&url).await?;

    if let Ok(resp_html) = resp.body_string().await {
        let table_selector = selector_error_convert!("table")?;
        let row_selector = selector_error_convert!("tr")?;
        let document = Html::parse_document(&resp_html);

        if let Some(table) = document.select(&table_selector).next() {
            let mut rows = table.select(&row_selector).skip(1);

            for row in rows {
                let mut hla_info: LigandInfo = row.text().collect::<Vec<&str>>().into();
                result.push(hla_info);
            }
            Ok(result)
        } else {
            Err(RetrieveLigandError::NoLigandTableFound(url))
        }
    } else {
        Err(RetrieveLigandError::InvalidHLA(hla.as_ref().to_string()))
    }
}

pub async fn obtain_hla_ligand_groups<T>(ligand_file: T) -> Result<(), Box<dyn std::error::Error>>
where
    T: AsRef<Path>,
{
    let mut lg = retrieve_ligand_group(&"A").await?;
    let mut b_lg = retrieve_ligand_group(&"B").await?;
    let mut c_lg = retrieve_ligand_group(&"C").await?;

    lg.append(&mut b_lg);
    lg.append(&mut c_lg);

    {
        let mut f = File::create(ligand_file)?;
        for ligand_info in lg {
            f.write(format!("{}\t{}\t{}\n", ligand_info.0, ligand_info.1, ligand_info.2).as_ref());
//                .await?;
        }
//        f.flush().await?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_hla() {
        let hla = "C*01:02";
        clean_hla(&hla).unwrap();
    }

    #[tokio::test]
    async fn test_reqwest() {
        let hla = "C*01:102";
        let result = retrieve_ligand_group(&hla).await.unwrap();

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
