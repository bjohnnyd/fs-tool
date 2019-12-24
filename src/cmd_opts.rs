use std::path::PathBuf;
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "default")]
pub struct Opt {
    #[structopt(long)]
    pub update_ligand_groups: bool,
}
