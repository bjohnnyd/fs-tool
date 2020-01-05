pub mod mhc {
    use crate::prelude::snafu_error::*;

    #[derive(Debug, Snafu)]
    pub enum Error {
        #[snafu(display("Incorrect HLA Nomenclature: {}", hla))]
        IncorrectHLA { hla: String },

        #[snafu(display(
            "Valid Gene Locus Names are A, B, C, DP, DM, DO, DQ, DR but the name provided was {}",
            locus
        ))]
        IncorrectGeneLocus { locus: String },

        #[snafu(display("HLA name too short: {}", hla))]
        HLANameTooShort { hla: String },

        #[snafu(display("The Ligand Group is unknown: {}", ligand_group))]
        UnknownLigandGroup { ligand_group: String },

        #[snafu(display(
            r#"IPD Frequency can either be "Common or Well Defined", "Rare", "Unknown" \
         but the provided frequency is {}"#,
            ipd_frequency
        ))]
        IPDFrequencyIncorrect { ipd_frequency: String },
    }
}

pub mod retrieve_ligands {
    use crate::prelude::snafu_error::*;

    #[derive(Debug, Snafu)]
    #[snafu(visibility = "pub(crate)")]
    pub enum Error {
        #[snafu(display(
            r#"\
        Please ensure that HLA alleles are of right format (e.g. "C*01:102" or "C*01:02")\
        for querying the IPD website, the passed HLA was {}"#,
            hla
        ))]
        IncorrectHLAInURL { hla: String },

        #[snafu(display(r#"Could Not Access IPD website at "{}": "{}""#, url, source))]
        WebsiteAccessError {
            source: attohttpc::Error,
            url: String,
        },

        #[snafu(display(r#"No Ligand Table found at {}"#, url))]
        NoLigandTableFound { url: String },

        #[snafu(display(r#"Could not read respocnse from "{}" : "{}""#, url, source))]
        CouldNotReadResponse {
            source: attohttpc::Error,
            url: String,
        },
    }
}
