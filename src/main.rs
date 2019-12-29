mod cmd_opts;
mod data;
mod mhc;
mod netmhcpan;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::{
    get_ligand_db_file, retrieve_ligand_group, save_ligand_groups, LigandInfo,
};
use directories::ProjectDirs;
use rayon::prelude::*;
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

    // Need to make this a function and switch to other branch and set ligand data if fails
    if let Some(ligand_db_file) = get_ligand_db_file() {
        if opt.update_ligand_groups || !ligand_db_file.exists() {
            save_ligand_groups(&ligand_db_file)?;
        }

        let mut f = File::open(ligand_db_file)?;
        let mut ligand_table = String::new();
        f.read_to_string(&mut ligand_table);
        println!("{}", ligand_table);
    } else {
        println!("Using backup data");
        println!("{}", LIGAND_TABLE);
    }

    Ok(())
}
