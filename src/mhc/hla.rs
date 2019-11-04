use crate::mhc::error::HLAError;

type Result<T> = std::result::Result<T, HLAError>;

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
    pub fn new(name: &str) -> Result<HLA> {
        let mut hla_name = name.trim_start_matches("HLA-").replace("*", "");
        hla_name.to_hla()

        let gene: mhc_meta::Locus = mhc_meta::hla_name_to_locus(&hla_name[..2])?;
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

        Ok(HLA {
            gene,
            allele_group,
            hla_protein,
            cds_synonymous_sub,
            non_coding_difference,
            expression_change,
            ligand_group: mhc_meta::LigandGroup::I,
            mhc_class: mhc_meta::MHC::I,
        })
    }
}

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
pub(crate) enum ExpressionChange {
    N,
    L,
    S,
    C,
    A,
    Q,
    Unknown,
}

pub(crate) fn hla_name_to_locus(hla_name: &str) -> Result<Locus> {
    hla_name[..2]
        .chars()
        .fold(Ok(Locus::Unknown), |mut locus, c| {
            match c {
                'A' => locus = Ok(Locus::A),
                'B' => locus = Ok(Locus::B),
                'C' => locus = Ok(Locus::C),
                'P' => locus = Ok(Locus::DP),
                'M' => locus = Ok(Locus::DM),
                'O' => locus = Ok(Locus::DO),
                'Q' => locus = Ok(Locus::DQ),
                'R' => locus = Ok(Locus::DR),
                _ => {
                    if c.is_digit(10) & (locus == Ok(Locus::Unknown)) {
                        locus = Err(HLAError)
                    } else {
                        ()
                    }
                }
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

