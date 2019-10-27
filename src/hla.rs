pub mod mhc {
    use super::*;

    fn parse_hla_digits(
        name: String,
        second_field_size: usize,
    ) -> (String, Option<String>, Option<String>, Option<String>) {
        let mut result = Vec::with_capacity(3);
        let allele_group = name.chars().take(2).collect();
        let name_size = name.len();

        for idx in 0..=2 {
            let (start, end) = match idx {
                0 => (2, 2 + second_field_size),
                _ => (
                    (2 + second_field_size) + (idx * 2),
                    (2 + second_field_size) + (idx * 2) + 2,
                ),
            };
            if end <= name_size {
                result.push(Some(name.chars().take(end).skip(start).collect()))
            } else {
                result.push(None)
            };
        }

        (
            allele_group,
            result[0].clone(),
            result[1].clone(),
            result[2].clone(),
        )
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
            let mut hla_name = name.trim_start_matches("HLA-").replace("*", "");

            let gene: mhc_meta::Locus = mhc_meta::hla_name_to_locus(&hla_name[..2]);
            let expression_change: mhc_meta::ExpressionChange =
                mhc_meta::hla_name_to_expression_change(hla_name.chars().last().unwrap());

            hla_name = hla_name.chars().filter(|char| char.is_numeric()).collect();
            let hla_name_size = hla_name.len();
            let second_field_size = match hla_name_size % 2 {
                0 => 2,
                _ => 3,
            };

            let (allele_group, hla_protein, cds_synonymous_sub, non_coding_difference) =
                parse_hla_digits(hla_name, second_field_size);

            HLA {
                gene,
                allele_group,
                hla_protein,
                cds_synonymous_sub,
                non_coding_difference,
                expression_change,
                ligand_group: mhc_meta::LigandGroup::I,
                mhc_class: mhc_meta::MHC::I,
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

    pub(crate) fn sanitize_hla_name(hla_name: &str) -> String {
        hla_name
            .trim_start_matches("HLA")
            .trim_start_matches("-")
            .replace("*", "")
    }

    pub(crate) fn hla_name_to_locus(hla_name: &str) -> Locus {
        hla_name[..2]
            .chars()
            .fold(Locus::Unknown, |mut locus, char| {
                match char {
                    'A' => locus = Locus::A,
                    'B' => locus = Locus::B,
                    'C' => locus = Locus::C,
                    'P' => locus = Locus::DP,
                    'M' => locus = Locus::DM,
                    'O' => locus = Locus::DO,
                    'Q' => locus = Locus::DQ,
                    'R' => locus = Locus::DR,
                    _ => (),
                }
                locus
            })
    }
    pub(crate) fn hla_name_to_expression_change(change_char: char) -> ExpressionChange {
        match change_char {
            'N' => ExpressionChange::N,
            'L' => ExpressionChange::L,
            'S' => ExpressionChange::S,
            'C' => ExpressionChange::C,
            'A' => ExpressionChange::A,
            'Q' => ExpressionChange::Q,
            _ => ExpressionChange::Unknown,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hla::mhc::HLA;
    use crate::hla::mhc_meta::hla_name_to_locus;

    #[test]
    fn test_hla_to_locus() {
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
