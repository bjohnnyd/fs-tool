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

    pub(crate) fn hla_name_to_locus(hla_name: &str) -> mhc_meta::Locus {
        let mut locus: mhc_meta::Locus = mhc_meta::Locus::Unknown;
        hla_name[..2].chars().fold(mhc_meta::Locus::Unknown, |mut locus, char| {
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
    pub(crate) fn hla_name_to_expression_change(change_char: char) -> mhc_meta::ExpressionChange {
            let mut e_change: mhc_meta::ExpressionChange = mhc_meta::ExpressionChange::Unknown;
            match change_char{
                    'N' => mhc_meta::ExpressionChange::N,
                    'L' => mhc_meta::ExpressionChange::L,
                    'S' => mhc_meta::ExpressionChange::S,
                    'C' => mhc_meta::ExpressionChange::C,
                    'A' => mhc_meta::ExpressionChange::A,
                    'Q' => mhc_meta::ExpressionChange::Q,
                    _ => mhc_meta::ExpressionChange::Unknown,
                }
    }

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) struct HLA {
        pub(crate) gene: mhc_meta::Locus,
        pub(crate) allele_group: String,
        pub(crate) hla_protein: Option<String>,
        pub(crate) cds_synonymous_sub: Option<String>,
        pub(crate) non_coding_difference: Option<String>,
        pub(crate) expression_change: mhc_meta::ExpressionChange,
        pub(crate) ligand_group: mhc_meta::LigandGroup,
        pub(crate) mhc_class: mhc_meta::MHC,
    }

    impl HLA {
        pub fn new(name: &str) -> HLA {
            let mut hla_name = sanitize_hla_name(name);

            let gene: mhc_meta::Locus = hla_name_to_locus(&hla_name[..2]);
            let expression_change: mhc_meta::ExpressionChange = hla_name_to_expression_change(hla_name.chars().last().unwrap());

            let hla_name: String = hla_name.chars().filter(|char| char.is_numeric()).collect();
            let hla_name_size = hla_name.len();
            let second_field_size = match hla_name_size % 2 {
                0 => 2,
                _ => 3,
            };

            let allele_group: String = hla_name.chars().take(2).collect();
            let hla_protein: Option<String> = if (2 + second_field_size) <= hla_name_size {
                Some(hla_name.chars().take(2 + second_field_size).skip(2).collect())
            } else {
                None
            };

            let cds_synonymous_sub: Option<String> = if (4 + second_field_size) <= hla_name_size {
                Some(hla_name.chars().take(4 + second_field_size).skip(2 + second_field_size).collect())
            } else {
                None
            };

            let non_coding_difference: Option<String> = if (6 + second_field_size) <= hla_name_size {
                Some(hla_name.chars().take(6 + second_field_size).skip(4 + second_field_size).collect())
            } else {
                None
            };

            HLA {
                gene,
                allele_group,
                hla_protein,
                cds_synonymous_sub,
                non_coding_difference,
                expression_change,
                ligand_group: mhc_meta::LigandGroup::I,
                mhc_class: mhc_meta::MHC::I
            }
        }
    }
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
        Unknown,
    }

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) enum MHC {
        I,
        II,
    }

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) enum ExpressionChange {
        N,
        L,
        S,
        C,
        A,
        Q,
        Unknown,
    }

    #[derive(Debug, Eq, PartialEq)]
    pub(crate) enum LigandGroup {
        I,
        II,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hla::mhc::{hla_name_to_locus, HLA};

    #[test]
    fn test_parse_hla_name() {
        let hla_name = mhc::sanitize_hla_name("HLA-A*03:01");
        assert_eq!(hla_name, "A03:01");
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

    #[test]
    fn test_create_hla() {
        let hla = HLA {
            gene: mhc_meta::Locus::A,
            allele_group: String::from("01"),
            hla_protein: Some("101".to_string()),
            cds_synonymous_sub: None,
            non_coding_difference: None,
            expression_change: mhc_meta::ExpressionChange::Unknown,
            ligand_group: mhc_meta::LigandGroup::I,
            mhc_class: mhc_meta::MHC::I,
        };
        assert_eq!(mhc::HLA::new("HLA-A*01:101"), hla);
    }
}
