use crate::prelude::fs_tool::{LigandGroup, HLA};
use crate::prelude::fs_trait::FindSimilarHLA;
use crate::prelude::traits::{FromStr, TryFrom};

//use serde::{Serialize, Deserialize};
use serde_derive::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct Genotype<'a> {
    #[serde(
        alias = "ID",
        alias = "Individual",
        alias = "Genotype",
        alias = "id",
        borrow
    )]
    pub id: Option<&'a str>,
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
    /* Need to implement desearializer for KIR genotype */
    #[serde(skip)]
    pub kir: Vec<KIR<'a>>,
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

impl<'a> Genotype<'a> {
    fn update_ligand_groups(&mut self, ligand_info: &Vec<HLA>) {
        update_if_present!(self, a1 a2 b1 b2 c1 c2, ligand_info);
    }
}

#[derive(Debug)]
pub enum Tail<'a> {
    Long(&'a str),
    Short(&'a str),
    Pseudo(&'a str),
    Unknown(&'a str),
}

impl<'a> From<&'a str> for Tail<'a> {
    fn from(s: &'a str) -> Self {
        match s.chars().next() {
            Some('L') => Tail::Long(s[1..].as_ref()),
            Some('S') => Tail::Short(s[1..].as_ref()),
            Some('P') => Tail::Pseudo(s[1..].as_ref()),
            _ => Tail::Unknown(s),
        }
    }
}

#[derive(Debug)]
pub enum Domain {
    Two,
    Three,
}

impl FromStr for Domain {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.chars().next() {
            Some('2') => Ok(Domain::Two),
            Some('3') => Ok(Domain::Three),
            None => Err("KIR name not complete missing domain information"),
            _ => Err(format!("KIR domain type not known {} expected 2/3", s).as_ref()),
        }
    }
}

#[derive(Debug)]
pub struct KIR<'a>(pub Domain, pub Tail<'a>);

impl<'a> FromStr for KIR<'a> {
    type Err = &'static str;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim_start_matches("KIR").trim_start_matches("kir");
        let domain = s[0..2].parse::<Domain>()?;
        let tail = Tail::from(s[2..].as_ref());
        Ok(KIR(domain, tail))
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
