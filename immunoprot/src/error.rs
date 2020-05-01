use thiserror::Error;
pub(crate) static NOMENCLATURE_URL: &str = "http://hla.alleles.org/nomenclature/naming.html";
pub(crate) static MHCI_LIGAND_MOTIFS: &str = "A11, A3, Bw4-80T, Bw4-80I, Bw6, C1, C2, Unclassified";
pub(crate) static EXPRESSION_CHANGES: &str = "N, L, S, C, A, Q and '' (blank)";
pub(crate) static HLA_GENES: &str = "A, B, C, DP, DM, DO, DQ and DR";

// Errors representing various NomenClature errors
#[derive(Debug, Error)]
pub enum NomenclatureError {
    #[error("Cound not determine the kir-ligand group from {0}.  Supported groups are {}", MHCI_LIGAND_MOTIFS)]
    UnknownLigandMotif(String),
    #[error("Unknown expression change tag '{0}'. ClassI nomenclature acceptable values are {}, for more details visit {}", EXPRESSION_CHANGES, NOMENCLATURE_URL)]
    UnknownExpressionChangeTag(String),
    #[error("Could not determine HLA Class I allele from '{0}'\n Colons (:) are necessary delimiters, other separators and placeholders are optional, visit {} for more details", NOMENCLATURE_URL)]
    CouldNotParseClassI(String),
    #[error("Could not determine gene from '{0}'. Recognized HLA genes are {}, for more details visit {}", HLA_GENES, NOMENCLATURE_URL)]
    GeneUnknown(String),
    #[error("Could not parse allele group information{0}, for naming guidelines see {}", NOMENCLATURE_URL)]
    NoAlleleGroup(String),
}


