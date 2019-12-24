mod cmd_opts;
mod data;

use structopt::StructOpt;
use tokio::fs::File;
use tokio::prelude::*;

use crate::cmd_opts::Opt;
use crate::data::retrieve_ligands::retrieve_ligand_group;

#[tokio::main]
async fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    if opt.update_ligand_groups {
        let mut lg = retrieve_ligand_group(&"A").await?;
        let mut b_lg = retrieve_ligand_group(&"B").await?;
        let mut c_lg = retrieve_ligand_group(&"C").await?;

        lg.append(&mut b_lg);
        lg.append(&mut c_lg);

        {
            let mut f = File::create("resources/lg.tsv").await?;
            for ligand_info in lg {
                f.write(
                    format!("{}\t{}\t{}\n", ligand_info.0, ligand_info.1, ligand_info.2).as_ref(),
                )
                .await?;
            }
            f.flush().await?;
        }
    }
    Ok(())
}
