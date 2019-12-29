mod cmd_opts;
mod data;
mod mhc;
mod netmhcpan;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::{get_ligand_table, parse_ligand_table, LigandInfo};
use crate::mhc::hla::HLA;
use directories::ProjectDirs;
use rayon::prelude::*;
use std::convert::TryFrom;
use std::fs::File;
use std::path::PathBuf;
use structopt::StructOpt;

use std::io::Read;

pub const LIGAND_TABLE: &'static str = include_str!("resources/2019-12-29_lg.tsv");

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .unwrap();

    let mut ligand_data = Vec::<LigandInfo>::new();

    if let Ok(ligand_table) = get_ligand_table(opt.update_ligand_groups) {
        ligand_data = parse_ligand_table(ligand_table);
    } else {
        ligand_data = parse_ligand_table(LIGAND_TABLE)
    }

    let hla = ligand_data
        .into_iter()
        .filter_map(|lg| HLA::try_from(lg).ok())
        .collect::<Vec<HLA>>();

    Ok(())
}
