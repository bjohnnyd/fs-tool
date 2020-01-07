pub use crate::prelude::snafu_error::*;

#[derive(Debug, Snafu)]
#[snafu(visibility = "pub")]
pub enum Error {
    #[snafu(display(r#"Incorrect HLA Nomenclature: {}"#, hla))]
    IncorrectHLA { hla: String },

    #[snafu(display(
        r#"Valid Gene Locus Names are A, B, C, DP, DM, DO, DQ, DR but the name provided was {}"#,
        locus
    ))]
    IncorrectGeneLocus { locus: String },

    #[snafu(display(r#"HLA name too short: {}"#, hla))]
    HLANameTooShort { hla: String },

    #[snafu(display(r#"The Ligand Group is unknown: {}"#, ligand_group))]
    UnknownLigandGroup { ligand_group: String },
    #[snafu(display(
        r#"IPD Frequency can either be "Common or Well Defined", "Rare", "Unknown" \
         but the provided frequency is {}"#,
        ipd_frequency
    ))]
    IPDFrequencyIncorrect { ipd_frequency: String },

    #[snafu(display(
        r#"Gene name has to contain at least one character found "{}""#,
        gene_name
    ))]
    GeneNameTooShort { gene_name: String },

    #[snafu(display(r#"Could Not Access IPD website at "{}": "{}""#, url, source))]
    WebsiteAccessError {
        source: attohttpc::Error,
        url: String,
    },

    #[snafu(display(r#"No Ligand Table found at {}"#, response))]
    NoLigandTableFound { response: String },

    #[snafu(display(r#"Could not read respocnse from "{}" : "{}""#, url, source))]
    CouldNotReadResponse {
        source: attohttpc::Error,
        url: String,
    },

    #[snafu(display(
        r#"Could not parse the content at Line: {} and Column {}"#,
        line,
        column
    ))]
    WebsiteParsingError { line: u32, column: u32 },

    #[snafu(display(r#"Could not access or create local data"#))]
    NoLocalData {},

    #[snafu(display(r#"No ligand information found at website {}"#, url))]
    NoLigandInformation { url: String },

    #[snafu(display(r#"Could not write ligand information "{}\t{}\t{}"\
         to file at path {}"#, info_1, info_2, info_3, p.display()))]
    WriteLigandInformationToFile {
        p: std::path::PathBuf,
        info_1: String,
        info_2: String,
        info_3: String,
        source: std::io::Error,
    },

    #[snafu(display(r#"Could not read ligand information from file "{}" : "{}""#, ligand_db_file.display(), source))]
    CouldNotReadLigandFile {
        ligand_db_file: std::path::PathBuf,
        source: std::io::Error,
    },

    #[snafu(display(r#"Could not open/create file "{}" : "{}""#, f_path.display(), source))]
    CouldNotOpenFile {
        f_path: std::path::PathBuf,
        source: std::io::Error,
    },

    #[snafu(display(
        r#"Could not create all the directories in path "{}" : "{}""#,
        out_dir.display(),
        source
        ))]
    CouldNotCreateDirectory {
        out_dir: std::path::PathBuf,
        source: std::io::Error,
    },

    #[snafu(display(r#"Could not start logging : "{}""#, source))]
    CouldNotStartLogging { source: log::SetLoggerError },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_could_not_open_file() {
        let p = std::path::PathBuf::new().join("blah");
        let result = std::fs::File::open(p.clone()).context(CouldNotOpenFile { f_path: p });
        match result {
            Ok(f) => (),
            Err(e) => {
                eprintln!("{}", e);
                ()
            }
        }
    }
}
