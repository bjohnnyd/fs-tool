pub(crate) static NOMENCLATURE_URL: &str = "http://hla.alleles.org/nomenclature/naming.html";
pub(crate) static MHCI_LIGAND_GROUPS: &str = "A11, A3, Bw4-80T, Bw4-80I, Bw6, C1, C2, Unclassified";
pub(crate) static EXPRESSION_CHANGES: &str = "N, L, S, C, A, Q and '' (blank)";

/// Represent error types created by parsing to MHCI
#[derive(Debug)]
pub enum MHCIError {
    UnknownLigandGroup(String),
    UnknownExpressionChangeTag(String),
    CouldNotParseMHCI(String),
    GeneUnknown(String),
}

impl std::fmt::Display for MHCIError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let help_string = format!("Visit {} for more info.", NOMENCLATURE_URL);
        match self {
            MHCIError::UnknownLigandGroup(lg_group) => writeln!(f, "Could not determine ligand group from '{}'\nAcceptable groups are: {}\n{}", lg_group, MHCI_LIGAND_GROUPS, help_string),
            MHCIError::UnknownExpressionChangeTag(exprs) => writeln!(f, "Unknown expression change tag '{}', acceptable values are {}\n{}", exprs, EXPRESSION_CHANGES, help_string),
            MHCIError::CouldNotParseMHCI(mhc_I) => writeln!(f, "Could not parse MHCI allele from '{}'\n Colons (:) are needed in MHCI names while other separators and placeholders are optional.\n{}", mhc_I,  help_string),
            MHCIError::GeneUnknown(mhc_I) => writeln!(f, "Could not parse locus information from {}\n{}", mhc_I, help_string),
        }
    }
}

impl std::error::Error for MHCIError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            _ => None,
        }
    }
}
