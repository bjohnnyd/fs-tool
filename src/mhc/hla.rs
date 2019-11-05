use crate::mhc::error::HLAErr;

type Result<T> = std::result::Result<T, HLAErr>;

#[derive(Debug, Eq, PartialEq)]
pub(crate) struct HLA {
    pub gene: Gene,
    pub allele_group: String,
    pub hla_protein: Option<String>,
    pub cds_synonymous_sub: Option<String>,
    pub non_coding_diff: Option<String>,
    pub expression_change: ExpressionChange,
    // TODO: Need to implement ligand groups
    //    pub ligand_group: mhc_meta::LigandGroup,
}

impl HLA {
    pub fn new(name: &str) -> Result<HLA> {
        let mut hla_name = name.trim_start_matches("HLA-").replace("*", "");
        hla_name.parse::<HLA>()
    }
}

impl std::str::FromStr for HLA {
    type Err = HLAErr;

    fn from_str(s: &str) -> Result<Self> {
        let s = s.trim_start_matches("HLA-").replace("*", "");

        let mut first_two_chars = s.chars().take(2);
        let gene = match first_two_chars.next().ok_or(HLAErr::GeneNameTooShort)? {
            'A' => Gene::A,
            'B' => Gene::B,
            'C' => Gene::C,
            'D' => match first_two_chars.next().ok_or(HLAErr::GeneNameTooShort)? {
                'P' => Gene::DP,
                'M' => Gene::DM,
                'O' => Gene::DO,
                'Q' => Gene::DQ,
                'R' => Gene::DR,
                _ => Err(HLAErr::ParseError(s.clone()))?,
            },
            _ => Err(HLAErr::ParseError(s.clone()))?,
        };

        // TODO: Should it be Unknown at all?
        let expression_change = match s.chars().last().ok_or(HLAErr::GeneNameTooShort)? {
            'N' => ExpressionChange::N,
            'L' => ExpressionChange::L,
            'S' => ExpressionChange::S,
            'C' => ExpressionChange::C,
            'A' => ExpressionChange::A,
            'Q' => ExpressionChange::Q,
            _ => ExpressionChange::Unknown,
        };

        let s_digits: String = s.chars().filter(|c| c.is_numeric()).collect();
        let hla_protein_field_size = match s_digits.len() % 2 {
            0 => 2,
            _ => 3,
        };

        let allele_group = s_digits.chars().take(2).collect();

        let mut fields = [None, None, None];
        let mut current_pos = 2;

        for (i, slot) in fields.iter_mut().enumerate() {
            let start = current_pos;
            let end = current_pos + hla_protein_field_size + 2 * i;
            current_pos = end;

            *slot = if end <= s_digits.len() {
                Some(s_digits.chars().take(end).skip(start).collect())
            } else {
                None
            };
        }

        let [hla_protein, cds_synonymous_sub, non_coding_diff] = fields;

        Ok(HLA {
            gene,
            allele_group,
            hla_protein,
            cds_synonymous_sub,
            non_coding_diff,
            expression_change,
        })
    }
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum Gene {
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
pub(crate) enum ExpressionChange {
    N,
    L,
    S,
    C,
    A,
    Q,
    Unknown,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_hla_from_str() {
        let hla = HLA {
            gene: Gene::A,
            allele_group: String::from("01"),
            hla_protein: Some(String::from("101")),
            cds_synonymous_sub: None,
            non_coding_diff: None,
            expression_change: ExpressionChange::Unknown,
        };

        assert_eq!("HLA-A*01:101".parse::<HLA>(), Ok(hla));
        //        assert_eq!("A*01:101".parse::<HLA>(), Ok(hla));
        //        assert_eq!("A01101".parse::<HLA>(), Ok(hla));
    }
}
