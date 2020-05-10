// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
mod calc;
mod cli;
mod error;
mod io;
mod meta;

pub const LIGAND_TABLE: &str = include_str!("resources/2019-12-29_lg.tsv");
pub const PROJECT_LIGAND_TABLE: &str = "kir_ligand.tsv";
pub const KIR_DEF: &str = "KIR:2,7,8,9";
pub const TCR_DEF: &str = "TCR:3,4,5,6,8,9";
pub const LOGGING_MODULES: [&str; 3] = ["immunoprot", "netmhcpan", "fstool"];
pub const DEFAULT_DELIM: u8 = b',';

use crate::cli::{Command::*, Opt};
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

    opt.set_logging();
    let kir_ligand_map = opt.setup_kir_ligand_info()?;
    let mut output_writers = opt.output_writers()?;

    match opt.cmd {
        Allele {
            binding_predictions,
            drop_default,
            peptide_length,
            measure,
        } => {
            let binding_data = read_raw_netmhcpan(binding_predictions)?;

            let allele_meta = create_allele_metadata(&binding_data, &kir_ligand_map);
            let binding_meta = create_binding_metadata(&binding_data);

            output_writers.write_allele_meta(&allele_meta)?;
            output_writers.write_binding_meta(&binding_meta)?;
        }
    }

    Ok(())
}
