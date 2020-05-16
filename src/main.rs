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
pub const TOOL_NAME: &str = "fs-tool";
pub const KIR_DEF: &str = "KIR:2,7,8,9";
pub const TCR_DEF: &str = "TCR:2,3,4,5,6,9";
pub const LOGGING_MODULES: [&str; 3] = ["immunoprot", "netmhcpan", "fs-tool"];
pub const DEFAULT_DELIM: u8 = b',';

use crate::calc::{calculate_fs, calculate_index_cohort_fs, create_calc_combs, IndexCache};
use crate::cli::{get_measures, print_defaults, Opt};
use crate::cohort::Individual;
use crate::io::reader::{read_kir_motif_binding, read_lilrb_scores, read_temp_cohort};
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
    let matches = Opt::clap().get_matches();

    if matches.is_present("settings") {
        return print_defaults();
    }

    let opt = Opt::from_args();
    opt.set_logging();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()?;

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

        let kir_motif_interactions = read_kir_motif_binding();
        let index_fs_cache = IndexCache::new(
            index_alleles,
            fs_result,
            &kir_motif_interactions,
            &measures,
            &opt.peptide_length,
        );
        let lilrb_scores = read_lilrb_scores();
        let cohort_result = calculate_index_cohort_fs(
            index_fs_cache,
            &cohort,
            &kir_motif_interactions,
            &lilrb_scores,
        );
        output_writers.write_cohort_result(&cohort_result)?;
    }
    Ok(())
}
