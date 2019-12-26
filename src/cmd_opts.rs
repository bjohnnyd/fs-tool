use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "fs-tool")]
pub struct Opt {
    #[structopt(long)]
    pub update_ligand_groups: bool,
    #[structopt(short, long, default_value = "4")]
    pub threads: usize,
}
