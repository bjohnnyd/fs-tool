use immunoprot::mhc::hla::ClassI;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct CohortTemp {
    #[serde(
        alias = "ID",
        alias = "sample",
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
}

#[cfg(test)]
mod tests {
    use crate::io::reader::read_temp_cohort;

    #[test]
    fn test_read_cohort() {
        let cohort = read_temp_cohort("tests/example_cohort.csv").unwrap();
        dbg!(&cohort);
    }
}
