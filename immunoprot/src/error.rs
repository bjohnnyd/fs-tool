pub(crate) static NOMENCLATURE_URL: &str = "http://hla.alleles.org/nomenclature/naming.html";
pub(crate) static MHCI_LIGAND_MOTIFS: &str = "A11, A3, Bw4-80T, Bw4-80I, Bw6, C1, C2, Unclassified";
pub(crate) static EXPRESSION_CHANGES: &str = "N, L, S, C, A, Q and '' (blank)";
pub(crate) static HLA_GENES: &str = "A, B, C, DP, DM, DO, DQ and DR";

/// Represent error types created by parsing to MHCI
#[derive(Debug)]
pub enum HLAError {
    UnknownLigandMotif(String),
    UnknownExpressionChangeTag(String),
    CouldNotParseMHCI(String),
    GeneUnknown(String),
    NoAlleleGroup(String),
}

impl std::fmt::Display for HLAError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let help_string = format!("Visit {} for more info.", NOMENCLATURE_URL);
        match self {
            HLAError::UnknownLigandMotif(lg_group) => {
                writeln!(
                    f, "Could not determine ligand group from '{}'\nAcceptable groups are: {}\n{}",
                    lg_group, MHCI_LIGAND_MOTIFS, help_string
                )
            },
            HLAError::UnknownExpressionChangeTag(exprs) => {
                writeln!(
                    f, "Unknown expression change tag '{}', acceptable values are {}\n{}",
                    exprs, EXPRESSION_CHANGES, help_string
                )
            },
            HLAError::CouldNotParseMHCI(mhc_I) => {
                writeln!(
                    f,
                    r#"Could not parse MHCI allele from '{}'\n Colons (:) are needed in MHCI names
                     while other separators and placeholders are optional.\n{}"#
                    ,mhc_I,  help_string
                )
            },
            HLAError::GeneUnknown(mhc_I) => {
                writeln!(
                    f, "Could not determine gene from '{}'. Recognized HLA genes are {}.\n{}",
                    mhc_I, HLA_GENES, help_string
                )
            },
            HLAError::NoAlleleGroup(mhc_I) => {
                writeln!(
                    f, "Could not parse allele group information{}\n{}",
                    mhc_I, help_string)
            },
        }
    }
}

impl std::error::Error for HLAError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            _ => None,
        }
    }
}
