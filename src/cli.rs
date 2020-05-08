use std::path::PathBuf;
use structopt::StructOpt;

use crate::calc::Measure;

pub const LIGAND_TABLE: &str = include_str!("resources/2019-12-29_lg.tsv");
pub const KIR_DEF: &str = "KIR:2,7,8,9";
pub const TCR_DEF: &str = "TCR:3,4,5,6,8,9";
pub const LOGGING_MODULES: [&str; 3] = ["immunoprot", "netmhcpan", "fstool"];

#[derive(Debug, StructOpt)]
#[structopt(
    name = "fstool",
    about = "Calculates fraction of shared bound motifs between HLA alleles while incorporating KIR ligand and LILRB binding information.",
    rename_all = "kebab-case"
)]
pub struct Opt {
    /// Determines verbosity of the processing, can be specified multiple times -vvv
    #[structopt(short, long, parse(from_occurrences))]
    verbose: u8,
    #[structopt(short, long)]
    /// Disables any information being printed to terminal (except errors)
    quiet: bool,
    #[structopt(short, long)]
    /// Updates the current kir ligand group data
    pub update: bool,
    #[structopt(short, long)]
    /// Returns the directory and files where the kir ligand data is stored
    location: bool,
    #[structopt(subcommand)]
    pub cmd: Command,
}

#[derive(Debug, StructOpt)]
pub enum Command {
    /// Allows calculation of fraction shared between Class I alleles from affinity predictions
    Allele {
        #[structopt(short, long, parse(from_os_str))]
        /// Path to file containing predicted Class I affinity data (NetMHCpan results)
        binding_predictions: Vec<PathBuf>,
        #[structopt(short, long, parse(from_os_str))]
        /// Location where to store allele vs allele fraction shared calculations
        output: Option<PathBuf>,
        #[structopt(long)]
        /// Whether to drop the default measures that are calculated by default
        /// on TCR and KIR motifs.
        drop_default: bool,
        #[structopt(short, long, value_delimiter = " ", default_value = "9")]
        /// Which length of input peptide sequence to consider
        peptide_length: Vec<usize>,
        /// Custom motif positions to use for calculations (format `Name:index,index..` e.g. KIR:2,7,8,9)
        #[structopt(short, long)]
        measure: Option<Vec<Measure>>,
    },
}

impl Opt {
    pub fn set_logging(&self) {
        use log::LevelFilter::{self, *};

        let log_level = if self.quiet {
            LevelFilter::Error
        } else {
            match self.verbose {
                level if level == 1 => Info,
                level if level == 2 => Debug,
                level if level > 2 => Trace,
                _ => Warn,
            }
        };

        env_logger::builder()
            .format_module_path(false)
            .filter_module(LOGGING_MODULES[0], log_level)
            .filter_module(LOGGING_MODULES[1], log_level)
            .filter_module(LOGGING_MODULES[2], log_level)
            .init()
    }
}
