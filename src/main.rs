mod cmd_opts;
mod data;
mod mhc;
mod netmhcpan;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::{obtain_hla_ligand_groups, retrieve_ligand_group, LigandInfo};
use rayon::prelude::*;
use structopt::StructOpt;

const LIGAND_TABLE: &'static str = include_str!("resources/2019-12-29_lg.tsv");

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .unwrap();
    if opt.update_ligand_groups {
        let hla = ["A", "B", "C"];
        //        obtain_hla_ligand_groups("resources/lg.tsv");

        let ligands = hla
            .par_iter()
            .filter_map(|gene| retrieve_ligand_group(gene).ok())
            .flatten()
            .collect::<Vec<LigandInfo>>();

        ligands.iter().for_each(|ligand| println!("{}", ligand))
    }
    Ok(())
}
