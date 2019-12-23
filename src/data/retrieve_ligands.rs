use crate::mhc::ligand_group::LigandGroup;

const IPD_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";


/// Makes GET requests to EBI to retrieve ligand information. A URL requests are of the form
/// `https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?<hla allele>` where `<hla allele>`
/// lacks the prefix (e.g. for `HLA-C*01:02` it would be `C*01:02`)

pub fn get_ligand_table(hla: &str) -> String {
    let sanitized_hla = hla.trim_start_matches("HLA-");
    let url = format!("{}{}", IPD_URL, sanitized_hla);

    let mut resp = reqwest::get(&url).unwrap();
//    println!("{}", resp.text().unwrap());
    if resp.status().is_success() {
        Document::from_read(resp)
            .unwrap()
            .find(Class("contenttable_lmenu"))
            .filter_map(|l| l.attr("tr"))
            .for_each(|x| println!("{}", x));
    }
    return "".to_string();
}

#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_request() {
        let hla = "C*01:02";

        let response = get_ligand_table(hla);

//        dbg!(response);

    }
}
