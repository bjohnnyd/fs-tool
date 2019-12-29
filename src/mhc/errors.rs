use std::fmt::{Error, Formatter};

#[derive(Debug)]
pub enum HLAError {
    IncorrectHLA(String),
    IncorrectGeneLocus(String),
    GeneNameTooShort,
    UnknownLigandGroup(String),
    IPDFrequencyIncorrect(String),
}

impl std::fmt::Display for HLAError {
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), Error> {
        match self {
            HLAError::IncorrectHLA(hla) => writeln!(f, "Incorrect HLA Nomenclature: {}", hla),
            HLAError::IncorrectGeneLocus(gene_locus) => writeln!(
                f,
                "Valid Gene Locus Names are\
                 A, B, C, DP, DM, DO, DQ, DR\
                 but the provided gene name was: {}",
                gene_locus
            ),
            HLAError::GeneNameTooShort => writeln!(f, "HLA gene name too short"),
            HLAError::UnknownLigandGroup(ligand_group) => {
                writeln!(f, "Unknown Ligand Group Specified: {}", ligand_group)
            }
            HLAError::IPDFrequencyIncorrect(freq) => {
                writeln!(f, "Unknown IPD Frequency Specified: {}", freq)
            }
        }
    }
}
