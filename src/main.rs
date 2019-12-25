mod cmd_opts;
mod data;

use structopt::StructOpt;
use tokio::fs::File;
use tokio::prelude::*;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::{obtain_hla_ligand_groups, retrieve_ligand_group};

#[tokio::main]
async fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    if opt.update_ligand_groups {
        obtain_hla_ligand_groups("resources/lg.tsv").await?;
    }
    Ok(())
}
