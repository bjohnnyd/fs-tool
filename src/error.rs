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
    #[error("Could not open cohort file.")]
    CouldNotOpenCohortFile(#[from] csv::Error),
}
