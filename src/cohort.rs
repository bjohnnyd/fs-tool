use immunoprot::ig_like::kir::Kir;
use immunoprot::ig_like::kir_ligand::{KirLigandMap, LigandMotif};
use immunoprot::mhc::hla::ClassI;
use serde::Deserialize;
use std::collections::HashMap;

// TODO: Need to implement a way to deal with cases where an allele in the genotype is missing
// TODO: Need to create an alternative cohort representation
// TODO: implement a macro for `CohortTemp` definition
#[derive(Debug, Deserialize)]
pub struct CohortTemp {
    #[serde(
        alias = "ID",
        alias = "sample",
        alias = "rowid",
        alias = "id",
        with = "serde_with::rust::display_fromstr"
    )]
    pub id: String,
    #[serde(
        alias = "A1",
        alias = "A.1",
        with = "serde_with::rust::display_fromstr"
    )]
    pub a1: ClassI,
    #[serde(
        alias = "A2",
        alias = "A.2",
        with = "serde_with::rust::display_fromstr"
    )]
    pub a2: ClassI,
    #[serde(
        alias = "B1",
        alias = "B.1",
        with = "serde_with::rust::display_fromstr"
    )]
    pub b1: ClassI,
    #[serde(
        alias = "B2",
        alias = "B.2",
        with = "serde_with::rust::display_fromstr"
    )]
    pub b2: ClassI,
    #[serde(
        alias = "C1",
        alias = "C.1",
        with = "serde_with::rust::display_fromstr"
    )]
    pub c1: ClassI,
    #[serde(
        alias = "C2",
        alias = "C.2",
        with = "serde_with::rust::display_fromstr"
    )]
    pub c2: ClassI,
    #[serde(
        alias = "KIR2DL1",
        alias = "2DL1",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2dl1: Option<bool>,
    #[serde(
        alias = "KIR2DL2",
        alias = "2DL2",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2dl2: Option<bool>,
    #[serde(
        alias = "KIR2DL3",
        alias = "2DL3",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2dl3: Option<bool>,
    #[serde(
        alias = "KIR2DL4",
        alias = "2DL4",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2dl4: Option<bool>,
    #[serde(
        alias = "KIR2DL5",
        alias = "2DL5",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2dl5: Option<bool>,
    #[serde(
        alias = "KIR2DS1",
        alias = "2DS1",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2ds1: Option<bool>,
    #[serde(
        alias = "KIR2DS2",
        alias = "2DS2",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2ds2: Option<bool>,
    #[serde(
        alias = "KIR2DS3",
        alias = "2DS3",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2ds3: Option<bool>,
    #[serde(
        alias = "KIR2DS4",
        alias = "2DS4",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2ds4: Option<bool>,
    #[serde(
        alias = "KIR2DS5",
        alias = "2DS5",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir2ds5: Option<bool>,
    #[serde(
        alias = "KIR3DS1",
        alias = "3DS1",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir3ds1: Option<bool>,
    #[serde(
        alias = "KIR3DL1",
        alias = "3DL1",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir3dl1: Option<bool>,
    #[serde(
        alias = "KIR3DL2",
        alias = "3DL2",
        deserialize_with = "crate::io::ser::optional_bool_deserialize",
        default
    )]
    pub kir3dl2: Option<bool>,
}

macro_rules! field_to_kirs {
    ($struct_name: ident, $($fname:ident),+) => {{
        let mut v = Vec::<Kir>::new();


        $(
        if let Some(field)  = $struct_name.$fname {
            if field {
                let name = stringify!($fname).to_uppercase();
                 v.push(name.parse::<Kir>().unwrap())
            }
        }
        )*

        v

    }
    };
}

impl From<CohortTemp> for Individual {
    fn from(cohort: CohortTemp) -> Self {
        let id = cohort.id;

        let kir_genotype = field_to_kirs!(
            cohort, kir2dl1, kir2dl2, kir2dl3, kir2dl4, kir2dl5, kir2ds1, kir2ds2, kir2ds3,
            kir2ds4, kir2ds5, kir3ds1, kir3dl1, kir3dl2
        );

        let hla_genotype = vec![
            cohort.a1,
            cohort.a2,
            cohort.b1.clone(),
            cohort.b2.clone(),
            cohort.c1.clone(),
            cohort.c2.clone(),
        ];

        Self {
            id,
            hla_genotype,
            kir_genotype,
        }
    }
}

#[derive(Debug)]
pub struct Individual {
    pub id: String,
    pub hla_genotype: Vec<ClassI>,
    pub kir_genotype: Vec<Kir>,
}

impl Individual {
    pub fn kir_bound_motifs<'a>(
        &self,
        motif_binding_map: &'a HashMap<Kir, Vec<LigandMotif>>,
    ) -> Vec<&'a LigandMotif> {
        self.kir_genotype
            .iter()
            .fold(Vec::<&LigandMotif>::new(), |mut bound_motifs, kir| {
                if let Some(motifs) = motif_binding_map.get(&kir) {
                    motifs.iter().for_each(|motif| bound_motifs.push(motif));
                }
                bound_motifs
            })
    }

    pub fn get_hla_motifs<'a>(&self, hla_ligand_map: &'a KirLigandMap) -> Vec<&'a LigandMotif> {
        self.hla_genotype
            .iter()
            .fold(Vec::<&LigandMotif>::new(), |mut hla_motifs, hla| {
                let mut info = hla_ligand_map.get_allele_info(&hla);
                info.sort();

                if let Some(info) = info.iter().next() {
                    hla_motifs.push(info.motif())
                };

                hla_motifs
            })
    }
}

#[cfg(test)]
mod tests {
    use crate::cohort::Individual;
    use crate::io::reader::{read_kir_motif_binding, read_temp_cohort};

    #[test]
    fn test_read_cohort() {
        let cohort = read_temp_cohort("tests/example_cohort.csv").unwrap();
    }

    #[test]
    fn test_create_individual() {
        let cohort = read_temp_cohort("tests/example_cohort.csv").unwrap();
        let individuals = cohort
            .into_iter()
            .map(|temp| Individual::from(temp))
            .collect::<Vec<Individual>>();
    }

    #[test]
    fn test_get_bound_motifs() {
        let cohort = read_temp_cohort("tests/example_cohort.csv").unwrap();
        let individuals = cohort
            .into_iter()
            .map(|temp| Individual::from(temp))
            .collect::<Vec<Individual>>();
        let motif_binding_map = read_kir_motif_binding();
    }
}
