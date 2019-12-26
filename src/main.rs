mod cmd_opts;
mod data;

use structopt::StructOpt;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::obtain_hla_ligand_groups;

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    if opt.update_ligand_groups {
        obtain_hla_ligand_groups("resources/lg.tsv");
    }
    Ok(())
}
