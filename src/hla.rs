use nom::{bytes::complete::tag, IResult};

pub mod mhc {
    use super::*;
    use crate::hla::mhc_meta::Locus::Unknown;

    pub(crate) fn sanitize_hla_name(hla_name: &str) -> String {
        hla_name
            .trim_start_matches("HLA")
            .trim_start_matches("-")
            .replace("*", "")
    }

    pub(crate) fn hla_name_to_vector<'a>(hla_name: &str) -> Vec<&str> {
            hla_name.split(":").collect()
    }

    pub(crate) fn hla_name_to_locus(hla_name: &str) -> mhc_meta::Locus {
        let mut locus: mhc_meta::Locus = mhc_meta::Locus::Unknown;
        hla_name.chars().fold(mhc_meta::Locus::Unknown, |mut locus, char | {
            match char {
                'A' => locus = mhc_meta::Locus::A,
                'B' => locus = mhc_meta::Locus::B,
                'C' => locus = mhc_meta::Locus::C,
                'P' => locus = mhc_meta::Locus::DP,
                'M' => locus = mhc_meta::Locus::DM,
                'O' => locus = mhc_meta::Locus::DO,
                'Q' => locus = mhc_meta::Locus::DQ,
                'R' => locus = mhc_meta::Locus::DR,
                _ => (),
            }
            locus
        })
    }

    pub struct HLA {
        gene: mhc_meta::Locus,
        allele_group: u8,
        hla_protein: u8,
        cds_synonymous_sub: u8,
        non_coding_difference: u8,
        expression_change: mhc_meta::ExpressionChange,
        ligand_group: mhc_meta::LigandGroup,
        mhc_class: mhc_meta::MHC,
    }

//        impl HLA {
//            pub fn new(name: &str) -> HLA {
//                let mut hla_name = sanitize_hla_name(name);
//                let locus: mhc_meta::Locus = hla_name_to_locus(&hla_name[..2]);
//                let hla_name_vec: Vec<&str> = if hla_name.contains(":") {
//                    hla_name.split(":").collect()
//                } else {
//                    vec!["BLA"]
//                }
//
//                return vec!["BLA"]
//
//            }
//        }
}

pub mod mhc_meta {

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) enum Locus {
        A,
        B,
        C,
        DP,
        DM,
        DO,
        DQ,
        DR,
        Unknown
    }

    pub(crate) enum MHC {
        I,
        II,
    }

    pub(crate) enum ExpressionChange {
        N,
        L,
        S,
        C,
        A,
        Q,
    }

    pub(crate) enum LigandGroup {
        I,
        II,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hla::mhc::{hla_name_to_vector, hla_name_to_locus};

    #[test]
    fn test_parse_hla_name() {
        let hla_name = mhc::sanitize_hla_name("HLA-A*03:01");
        assert_eq!(hla_name, "A03:01");
        let hla_name_split: Vec<&str> = mhc::hla_name_to_vector(hla_name.as_ref());
        assert_eq!(hla_name_split, vec!["A03", "01"] )
    }

    #[test]
    fn test_hla_to_locus() {
        println!("{:?}", hla_name_to_locus("A0301"));
        println!("{:?}", hla_name_to_locus("B0301"));
        println!("{:?}", hla_name_to_locus("C0301"));
        println!("{:?}", hla_name_to_locus("DP0301"));
        println!("{:?}", hla_name_to_locus("DR12"));
        println!("{:?}", hla_name_to_locus("DQ11"));
        println!("{:?}", hla_name_to_locus("DQ20"));
        println!("{:?}", hla_name_to_locus("01301"));
        assert_eq!(hla_name_to_locus("A0301"), mhc_meta::Locus::A);
        assert_eq!(hla_name_to_locus("0301"), mhc_meta::Locus::Unknown);
    }
}
