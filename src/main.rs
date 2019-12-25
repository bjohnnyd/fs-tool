mod cmd_opts;
mod data;

use structopt::StructOpt;
/// https://docs.rs/async-std/1.4.0/async_std/index.html#features
/// (effect on performance alternative below):
/// https://docs.rs/async-std/1.4.0/async_std/task/fn.block_on.html
use async_std::task;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::{obtain_hla_ligand_groups, retrieve_ligand_group};

async fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    if opt.update_ligand_groups {
        obtain_hla_ligand_groups("resources/lg.tsv").await?;
    }
    Ok(())
}
