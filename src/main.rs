// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
mod calc;
mod cli;
mod error;
mod reader;
mod writer;

pub const LIGAND_TABLE: &str = include_str!("resources/2019-12-29_lg.tsv");
pub const PROJECT_LIGAND_TABLE: &str = "kir_ligand.tsv";
pub const KIR_DEF: &str = "KIR:2,7,8,9";
pub const TCR_DEF: &str = "TCR:3,4,5,6,8,9";
pub const LOGGING_MODULES: [&str; 3] = ["immunoprot", "netmhcpan", "fstool"];

use crate::cli::{Command::*, Opt};
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

    opt.set_logging();
    let kir_ligand_map = opt.setup_kir_ligand_info()?;

    match opt.cmd {
        Allele {
            binding_predictions,
            output,
            drop_default,
            peptide_length,
            measure,
        } => {
            let binding_data = binding_predictions
                .iter()
                .map(read_raw_netmhcpan)
                .collect::<Result<Vec<_>, _>>()?;
        }
    }

    Ok(())
}
