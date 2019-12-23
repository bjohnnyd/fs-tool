use crate::data::errors::RetrieveLigandError::{self, InvalidHLA};
use scraper::{Html, Selector};

const IPD_KIR_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";

fn is_ipd_search_safe<T>(hla: T) -> bool
where
    T: AsRef<str>,
{
    if hla.as_ref().len() == 1 {
        match hla.as_ref() {
            "A" | "B" | "C" => true,
            _ => false,
        }
    } else {
        hla.as_ref().contains('*')
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

pub fn retrieve_ligand_group<T>(hla: &T) -> Result<Vec<String>, RetrieveLigandError>
where
    T: AsRef<str>,
{
    let url = format!("{}{}", IPD_KIR_URL, &clean_hla(&hla)?);

    let mut resp = reqwest::get(&url)?;

    if let Ok(resp_html) = resp.text() {
        let table_selector = Selector::parse("table").unwrap();
        let row_selector = Selector::parse("tr").unwrap();
        let document = Html::parse_document(&resp_html);

        let table = document.select(&table_selector).next().unwrap();
        let mut rows = table.select(&row_selector).skip(1);

        //        while let Some(row) = rows.next() {
        for row in rows {
            let hla_info = row.text().collect::<Vec<&str>>();
            //            println!("HLA: {}, LigandGroup: {}, Frequency: {}", hla_info[0], hla_info[1], hla_info[2]);
            println!("{:#?}", hla_info);
        }
        Ok(table
            .text()
            .map(|val| val.to_string())
            .collect::<Vec<String>>())
    } else {
        Err(RetrieveLigandError::InvalidHLA(hla.as_ref().to_string()))
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn invalid_hla() {
        let hla = "C*01:02";
        clean_hla(&hla).unwrap();
    }

    #[test]
    fn test_reqwest() {
        let hla = "C*01:102";
        let result = retrieve_ligand_group(&hla).unwrap();
    }
}
