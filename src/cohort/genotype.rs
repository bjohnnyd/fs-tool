use crate::prelude::fs_tool::{LigandGroup, HLA};
use crate::prelude::fs_trait::FindSimilarHLA;
use crate::prelude::traits::TryFrom;

//use serde::{Serialize, Deserialize};
use serde_derive::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Genotype {
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "A1",
        alias = "A.1",
        alias = "a.1"
    )]
    pub a1: HLA,
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "A2",
        alias = "A.2",
        alias = "a.2"
    )]
    pub a2: HLA,
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "B1",
        alias = "B.1",
        alias = "b.1"
    )]
    pub b1: HLA,
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "B2",
        alias = "B.2",
        alias = "b.2"
    )]
    pub b2: HLA,
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "C1",
        alias = "C.1",
        alias = "c.1"
    )]
    pub c1: HLA,
    #[serde(
        with = "serde_with::rust::display_fromstr",
        alias = "C2",
        alias = "C.2",
        alias = "c.2"
    )]
    pub c2: HLA,
    #[serde(skip)]
    pub ligand_groups: Option<LigandGroup>,
}

macro_rules! update_if_present {
   ($self:ident, $($hla:ident) +, $lg_info:ident) => {{
    $(
    if let Some(lg_info_hla) = $lg_info.get_hla(&$self.$hla) {
    $self.$hla.set_ligand_info(&lg_info_hla);
    })+
   }};
}

impl Genotype {
    fn update_ligand_groups(&mut self, ligand_info: &Vec<HLA>) {
        update_if_present!(self, a1 a2 b1 b2 c1 c2, ligand_info);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_reading_in() {
        let mut rdr = csv::Reader::from_path("examples/example_cohort.csv").unwrap();
        for result in rdr.deserialize() {
            let genotype: Genotype = result.unwrap();
            println!("{:#?}", genotype);
        }
    }

    #[test]
    fn test_update_meta() {
        let LIGAND_TABLE: &str = include_str!("../resources/2019-12-29_lg.tsv");
        let ligand_data = crate::data::retrieve_ligands::parse_ligand_table(LIGAND_TABLE);
        let ligand_data = ligand_data
            .into_iter()
            .map(|lg| {
                if let Ok(hla) = HLA::try_from(lg.clone()) {
                    Some(hla)
                } else {
                    None
                }
            })
            .filter_map(|hla| hla)
            .collect();

        let mut rdr = csv::Reader::from_path("examples/example_cohort.csv").unwrap();
        for result in rdr.deserialize() {
            let mut genotype: Genotype = result.unwrap();
            genotype.update_ligand_groups(&ligand_data);
            println!("{:#?}", genotype);
        }
    }
}
