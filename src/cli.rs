const PEPTIDE_LENGTHS: [&str; 4] = ["8", "9", "10", "11"];

use std::path::PathBuf;
use structopt::StructOpt;

use crate::calc::Measure;
use crate::error::Error;
use crate::io::reader::*;
use crate::io::writer::*;
use crate::LOGGING_MODULES;

use crate::meta::{AlleleMeta, BindingMeta};
use immunoprot::ig_like::kir_ligand::KirLigandMap;
use log::{info, warn};

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
    #[structopt(short, long)]
    /// Directory to store outputs
    output: PathBuf,
    /// Prefix to assign to all outputs
    #[structopt(long)]
    prefix: Option<String>,
    #[structopt(subcommand)]
    pub cmd: Command,
}

#[derive(Debug, StructOpt)]
pub enum Command {
    /// Allows calculation of fraction shared between Class I alleles from affinity predictions
    Allele {
        #[structopt(short, long, parse(from_os_str))]
        /// Path to file containing predicted Class I affinity data (NetMHCpan results)
        binding_predictions: PathBuf,
        #[structopt(long)]
        /// Whether to drop the default measures that are calculated by default
        /// on TCR and KIR motifs.
        drop_default: bool,
        #[structopt(short, long, possible_values = &PEPTIDE_LENGTHS, default_value = "9")]
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

    pub fn setup_kir_ligand_info(&self) -> Result<KirLigandMap, crate::error::Error> {
        if self.update {
            info!("Updating kir ligand information");
            let update = KirLigandMap::updated().ok();
            match update {
                Some(updated_info) => {
                    write_project_ligand_info(&updated_info);
                } ,
                _ => warn!("Failed to obtain updated kir ligand information from IPD/EBI.\n Trying to use the most recent downloaded version...")
            }
        }
        let kir_ligand_map = read_project_ligand_info();

        match kir_ligand_map {
            Some(ligand_map) => Ok(ligand_map),
            _ => Err(crate::error::Error::KirLigandMapError),
        }
    }

    pub fn output_writers(&self) -> std::result::Result<OutputWriters, Error> {
        let prefix = match &self.prefix {
            Some(prefix) => format!("{}_", prefix),
            _ => String::new(),
        };

        let output_dir = &self.output;

        if !output_dir.exists() {
            std::fs::create_dir_all(output_dir).or_else(|_| Err(Error::CouldNotCreateOutputDir))?
        }

        let allele_path = output_dir.join(format!("{}allele_metadata.csv", prefix));
        let binding_path = output_dir.join(format!("{}allele_binding_summary.csv", prefix));

        let allele_meta = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(crate::DEFAULT_DELIM)
            .from_path(allele_path)
            .or_else(|_| Err(Error::CouldNotCreateOutputFile))?;

        let binding_meta = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(crate::DEFAULT_DELIM)
            .from_path(binding_path)
            .or_else(|_| Err(Error::CouldNotCreateOutputFile))?;

        Ok(OutputWriters {
            allele_meta,
            binding_meta,
        })
    }
}

pub struct OutputWriters {
    pub allele_meta: csv::Writer<std::fs::File>,
    pub binding_meta: csv::Writer<std::fs::File>,
}

impl OutputWriters {
    pub fn write_allele_meta(
        &mut self,
        metadata: &[AlleleMeta],
    ) -> std::result::Result<Vec<()>, Error> {
        let write_result = metadata
            .iter()
            .map(|meta| self.allele_meta.serialize(meta))
            .collect::<Result<Vec<_>, _>>();

        Ok(write_result.or_else(|_| Err(Error::CouldNotWriteAlleleMeta))?)
    }

    pub fn write_binding_meta(
        &mut self,
        metadata: &[BindingMeta],
    ) -> std::result::Result<Vec<()>, Error> {
        let write_result = metadata
            .iter()
            .map(|meta| self.binding_meta.serialize(meta))
            .collect::<Result<Vec<_>, _>>();

        Ok(write_result.or_else(|_| Err(Error::CouldNotWriteBindingMeta))?)
    }
}
