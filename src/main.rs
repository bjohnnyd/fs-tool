// TODO: If binding predicitons file does not exist the error is not good;
// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
mod calc;
mod cli;
mod cohort;
mod error;
mod io;
mod meta;

pub static KIR_MOTIF_BINDING: &str = include_str!("resources/kir_motif_binding.tsv");
pub static LILRB_SIMSCORES: &str = include_str!("resources/lilrb_simscores.tsv");
pub static KIR_SIMSCORES: &str = include_str!("resources/kir_ligand_simscores.tsv");
pub const LIGAND_TABLE: &str = include_str!("resources/allele_motifs.tsv");
pub const PROJECT_LIGAND_TABLE: &str = "allele_motifs.tsv";
pub const KIR_DEF: &str = "KIR:2,7,8,9";
pub const TCR_DEF: &str = "TCR:2,3,4,5,6,9";
pub const LOGGING_MODULES: [&str; 3] = ["immunoprot", "netmhcpan", "fstool"];
pub const DEFAULT_DELIM: u8 = b',';

use crate::calc::{calculate_fs, create_calc_combs, create_index_fs_map};
use crate::cli::{get_measures, Opt};
use crate::cohort::Individual;
use crate::io::reader::read_temp_cohort;
use crate::meta::{create_allele_metadata, create_binding_metadata};

use netmhcpan::reader::read_raw_netmhcpan;
use structopt::StructOpt;

fn main() -> std::result::Result<(), ()> {
    match main_try() {
        Ok(()) => Ok(()),
        Err(e) => {
            eprintln!("{}", e);
            Err(())
        }
    }
}

fn main_try() -> Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();

    // if opt.location {
    //     if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
    //         println!("Current allele motif information is stored at {}", project_dir.data_dir().join(PROJECT_LIGAND_TABLE).display())
    //     } else {
    //         error!("Could not determine the global config directory");
    //         return Err(Box::new(crate::error::Error::NoGlobalConfigDir))
    //     }
    //     return Ok(())
    // }

    opt.set_logging();

    let kir_ligand_map = opt.setup_kir_ligand_info()?;
    let mut output_writers = opt.output_writers()?;

    let binding_data = read_raw_netmhcpan(opt.binding_predictions)?;

    let allele_meta = create_allele_metadata(&binding_data, &kir_ligand_map);
    let binding_meta = create_binding_metadata(&binding_data);

    output_writers.write_allele_meta(&allele_meta)?;
    output_writers.write_binding_meta(&binding_meta)?;

    let allele_combs = create_calc_combs(&binding_data);
    let measures = get_measures(opt.measure, opt.drop_default);

    let fs_result = calculate_fs(
        &allele_combs,
        &measures,
        &kir_ligand_map,
        &opt.peptide_length,
        binding_data.weak_threshold(),
        opt.unique,
    );
    output_writers.write_fs_result(&fs_result)?;

    if let (Some(index_alleles), Some(cohort_path)) = (opt.index, opt.cohort) {
        let cohort = read_temp_cohort(cohort_path)?
            .into_iter()
            .map(Individual::from)
            .collect::<Vec<Individual>>();

        let index_map = create_index_fs_map(index_alleles, fs_result);
        dbg!(index_map);
    }
    Ok(())
}
