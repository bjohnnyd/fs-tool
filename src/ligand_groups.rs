const IPD_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";

pub mod ligands_groups {
    use super::*;
    use select::document::Document;
    use select::predicate::Name;

    enum P80 {
        I,
        T,
        N,
        K,
    }

    enum Bw4 {
        P80,
    }

    enum C1 {
        P80,
    }

    enum C2 {
        P80,
    }

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) enum LigandGroup {
        A3,
        A11,
        Bw4,
        Bw6,
        C1,
        C2,
        Unclassified,
        Unknown,
    }

    pub fn get_ligand_group(hla: &str) -> String {
        let url = format!("{}{}", IPD_URL, hla);
        println!("URL = {}", url);

        let mut resp = reqwest::get(&url).unwrap();
        println!("{}", resp.text().unwrap());
        if resp.status().is_success() {
            Document::from_read(resp)
                .unwrap()
                .find(Name("table"))
                .filter_map(|l| l.attr("tr"))
                .for_each(|x| println!("{}", x));
        }
        return "".to_string();
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_ligand_group() {
        assert_eq!(ligands_groups::get_ligand_group("C*01:02"), "TEST");
    }
}
