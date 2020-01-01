use ::fs_tool::prelude::external::{Opt, StructOpt};
use ::fs_tool::prelude::fs_tool::{
    get_ligand_table, parse_ligand_table, read_netmhcpan, LigandInfo, NetMHCpanSummary, HLA,
};
use ::fs_tool::prelude::io::*;
use ::fs_tool::prelude::traits::TryFrom;
use std::io::BufRead;

pub const LIGAND_TABLE: &str = include_str!("resources/2019-12-29_lg.tsv");

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .unwrap();

    let mut ligand_data = Vec::<LigandInfo>::new();

    if let Ok(ligand_table) = get_ligand_table(opt.update_ligand_groups) {
        ligand_data = parse_ligand_table(ligand_table);
    } else {
        ligand_data = parse_ligand_table(LIGAND_TABLE)
    }

    let hla = ligand_data
        .into_iter()
        .filter_map(|lg| HLA::try_from(lg).ok())
        .collect::<Vec<HLA>>();

    if let Some(netmhcout) = opt.netmhcpan {
        let f = File::open(netmhcout)?;
        let netmhcpan_summary = read_netmhcpan(f);
        dbg!(netmhcpan_summary);
    }

    Ok(())
}
