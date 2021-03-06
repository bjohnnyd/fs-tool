use thiserror::Error;
pub(crate) static NOMENCLATURE_URL: &str = "http://hla.alleles.org/nomenclature/naming.html";
pub(crate) static MHCI_LIGAND_MOTIFS: &str = "A11, A3, Bw4-80T, Bw4-80I, Bw6, C1, C2, Unclassified";
pub(crate) static EXPRESSION_CHANGES: &str = "N, L, S, C, A, Q and '' (blank)";
pub(crate) static HLA_GENES: &str = "A, B, C, DP, DM, DO, DQ and DR";

use crate::ig_like::kir_ligand::IPD_KIR_URL;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum NomenclatureError {
    /* HLA Class I related */
    #[error(
        "Cound not determine the kir-ligand group from {0}.  Supported groups are {}",
        MHCI_LIGAND_MOTIFS
    )]
    #[doc(hidden)]
    UnknownLigandMotif(String),
    #[error("Unknown expression change tag '{0}'. ClassI nomenclature acceptable values are {}, for more details visit {}", EXPRESSION_CHANGES, NOMENCLATURE_URL)]
    #[doc(hidden)]
    UnknownExpressionChangeTag(String),
    #[error("Could not determine HLA Class I allele from '{0}'\n Colons (:) are necessary delimiters, other separators and placeholders are optional, visit {} for more details", NOMENCLATURE_URL)]
    #[doc(hidden)]
    CouldNotParseClassI(String),
    #[error("Could not determine gene from '{0}'. Recognized HLA genes are {}, for more details visit {}", HLA_GENES, NOMENCLATURE_URL)]
    #[doc(hidden)]
    GeneUnknown(String),
    #[error(
        "Could not parse allele group information{0}, for naming guidelines see {}",
        NOMENCLATURE_URL
    )]
    #[doc(hidden)]
    NoAlleleGroup(String),
    #[error("Empty string passed as an HLA allele please check your naming")]
    #[doc(hidden)]
    EmptyAlleleString,

    /* KIR related */
    #[error("The KIR type has an unknown tail. Tails can be either S, L or P but got '{0}'")]
    #[doc(hidden)]
    UnknownKirTail(String),
    #[error("The KIR type has an unknown domain. Domain can be either 2D or 3D but got '{0}'")]
    #[doc(hidden)]
    UnknownKirDomain(String),
    #[error("The KIR type protein is not known. Accepted KIR proteins are 1, 2, 3, 4, 5, 5A, 5B but got '{0}'")]
    #[doc(hidden)]
    UnknownKirProtein(String),
    #[error("Allele has insufficient information. A kir allele has to at least have a series specified (e.g. 003) but got '{0}'")]
    #[doc(hidden)]
    UnknownKirAllele(String),
}

/// Errors related to parsing IPD website
#[derive(Debug, Error)]
pub enum HtmlParseError {
    #[error("Could not create parser for HTML table")]
    #[doc(hidden)]
    CouldNotCreateParserForTable,
    #[error("Could not connect to {}", IPD_KIR_URL)]
    #[doc(hidden)]
    CouldNotConnect(#[from] attohttpc::Error),
    #[error("Could not read response from {}", IPD_KIR_URL)]
    #[doc(hidden)]
    CouldNotReadResponse(attohttpc::Error),
    #[error("Could not read HLA Class I allele from table row:\n{0}")]
    #[doc(hidden)]
    CouldNotReadClassI(String),
    #[error("Could not read motif from table row:\n{0}")]
    #[doc(hidden)]
    CouldNotReadMotif(String),
    #[error("Could not read frequency from table row:\n{0}")]
    #[doc(hidden)]
    CouldNotReadFreq(String),
    #[error("Expected 3 columns in HLA allele table but found {0}:\n")]
    #[doc(hidden)]
    IncorrectNumberOfColumns(usize, String),
}

/// Errors related to reading data from files
#[derive(Debug, Error)]
pub enum IoError {
    #[error("Could not read HLA allele at line {0}, column 1")]
    #[doc(hidden)]
    CouldNotReadAllele(usize),
    #[error("Could not read KIR Ligand Motif line {0}, column 2")]
    #[doc(hidden)]
    CouldNotReadMotif(usize),
    #[error("Could not read or open Kir Ligand Info file")]
    #[doc(hidden)]
    CouldNotReadOrOpenFile(#[from] csv::Error),
    #[error("Could not read line in Kir Ligand Info Table")]
    #[doc(hidden)]
    CouldNotParseLine,
}
