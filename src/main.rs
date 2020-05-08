// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]

mod calc;
mod cli;
mod error;

use crate::cli::{Command::*, Opt};
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
    let update = opt.update;

    match opt.cmd {
        Allele {
            binding_predictions,
            output,
            drop_default,
            peptide_length,
        } => {}
    }

    Ok(())
}
