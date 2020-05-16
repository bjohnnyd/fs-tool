const PEPTIDE_LENGTHS: [&str; 4] = ["8", "9", "10", "11"];

use std::path::PathBuf;
use structopt::StructOpt;

use crate::calc::Measure;
use crate::error::Error;
use crate::io::reader::*;
use crate::io::writer::*;
use crate::{KIR_DEF, LOGGING_MODULES, PROJECT_LIGAND_TABLE, TCR_DEF};

use immunoprot::ig_like::kir_ligand::KirLigandMap;
use immunoprot::mhc::hla::ClassI;
use log::{info, warn};

#[derive(Debug, StructOpt)]
#[structopt(
    name = "fstool",
    about = "Calculates fraction of shared bound motifs between HLA alleles while incorporating KIR ligand and LILRB binding information.",
    rename_all = "kebab-case"
)]
pub struct Opt {
    /// Number of threads
    #[structopt(short, long, default_value = "4")]
    pub threads: usize,
    /// Determines verbosity of the processing, can be specified multiple times -vvv
    #[structopt(short, long, parse(from_occurrences))]
    verbose: u8,
    #[structopt(short, long)]
    /// Disables any information being printed to terminal (except errors)
    quiet: bool,
    #[structopt(long)]
    /// Updates the current kir ligand group data
    pub update: bool,
    #[structopt(long)]
    /// Lists default measure names and motif positions as well as the default location
    /// updated kir ligand will be stored
    pub settings: bool,
    #[structopt(short, long, required_unless = "settings")]
    /// Directory to store outputs
    output: PathBuf,
    /// Prefix to assign to all outputs
    #[structopt(long)]
    prefix: Option<String>,
    #[structopt(short, long, parse(from_os_str), required_unless = "settings")]
    /// Path to file containing predicted Class I affinity data (NetMHCpan results)
    pub binding_predictions: Vec<PathBuf>,
    #[structopt(long)]
    /// Drop default measures based on TCR and KIR motifs.
    pub drop_default: bool,
    #[structopt(short, long, possible_values = &PEPTIDE_LENGTHS, default_value = "9")]
    /// Which length of input peptide sequence to consider
    pub peptide_length: Vec<usize>,
    /// Custom motif positions to use for calculations (format `Name:index,index..` e.g. KIR:2,7,8,9)
    #[structopt(short, long)]
    pub measure: Option<Vec<Measure>>,
    /// Whether only unique peptide/motif sequences should be considered in the calculations
    #[structopt(short, long)]
    pub unique: bool,
    /// Index allele used for cohort calculations only, all individuals will be compared to these alleles
    #[structopt(short, long, requires = "cohort")]
    pub index: Option<Vec<ClassI>>,
    /// Cohort of individuals for which all measures will be calculated
    #[structopt(short, long, requires = "index")]
    pub cohort: Option<PathBuf>,
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
        let allele_fs_path = output_dir.join(format!("{}allele_fs_result.csv", prefix));
        let cohort_result_path = output_dir.join(format!("{}cohort_result.csv", prefix));

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

        let allele_fs_result = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(crate::DEFAULT_DELIM)
            .from_path(allele_fs_path)
            .or_else(|_| Err(Error::CouldNotCreateOutputFile))?;

        let cohort_result = match self.cohort {
            Some(_) => Some(
                csv::WriterBuilder::new()
                    .has_headers(true)
                    .delimiter(crate::DEFAULT_DELIM)
                    .from_path(cohort_result_path)
                    .or_else(|_| Err(Error::CouldNotCreateOutputFile))?,
            ),
            _ => None,
        };

        Ok(OutputWriters {
            allele_meta,
            binding_meta,
            allele_fs_result,
            cohort_result,
        })
    }
}

pub fn print_defaults() -> Result<(), Box<dyn std::error::Error>> {
    println!("Default Measures:\n  - {}\n  - {}\n", TCR_DEF, KIR_DEF);

    if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
        let data_file = project_dir.data_dir().join(PROJECT_LIGAND_TABLE);
        let s = if data_file.exists() {
            "Data is"
        } else {
            "Updated data will be"
        };
        println!("{} stored at:  '{}'", s, data_file.display());
        Ok(())
    } else {
        Err("Could not create/or read the global data directory".into())
    }
}
pub fn get_measures(measures: Option<Vec<Measure>>, drop: bool) -> Vec<Measure> {
    let mut measures = measures.unwrap_or_default();

    if !drop {
        measures.push(TCR_DEF.parse().unwrap());
        measures.push(KIR_DEF.parse().unwrap());
    }

    measures
}
