use thiserror::Error;

/// Errors related to immune protein nomenclature
#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not parse motif positions")]
    IncorrectMotifPositions(#[from] std::num::ParseIntError),
    #[error("Could not load any kir ligand information")]
    KirLigandMapError,
    #[error("Could not create output directory")]
    CouldNotCreateOutputDir,
    #[error("Could not create output file")]
    CouldNotCreateOutputFile,
    #[error("Could not write allele metadata to output.")]
    CouldNotWriteAlleleMeta,
    #[error("Could not write binding metadata to output.")]
    CouldNotWriteBindingMeta,
    #[error("Could not write fraction shared results to output.")]
    CouldNotWriteFsResult,
    #[error("Could not write cohort calculation results to output")]
    CouldNotWriteCohortResult,
    #[error("Could not open cohort file:\n{:?}", .0.kind().display())]
    CouldNotOpenCohortFile(#[from] csv::Error),
    #[error("No global config directory exists")]
    NoGlobalConfigDir,
    #[error("Unsorted vector of peptide indices found")]
    UnsortedIndices,
    #[error("Measures have to be of the format <name:idx,idx,idx...> got {0}")]
    IncorrectMeasure(String),
}

trait ErrorKindDisplay {
    fn display(&self) -> String;
}

impl ErrorKindDisplay for csv::ErrorKind {
    fn display(&self) -> String {
        use csv::DeserializeErrorKind::*;
        use csv::ErrorKind::*;
        match self {
            Deserialize { pos: Some(i), err } => {
                let n = i.line();
                let message = match err.kind() {
                    Message(s) | Unsupported(s) => s.to_string(),
                    _ => "".to_string(),
                };

                format!("{}, line {}", message, n)
            }
            _ => "".to_string(),
        }
    }
}
