use ::fs_tool::prelude::external::{Opt, StructOpt};
use ::fs_tool::prelude::fs_tool::{
    get_ligand_table, parse_ligand_table, read_netmhcpan, LigandInfo, Measure, NetMHCpanSummary,
    HLA,
};
use ::fs_tool::prelude::io::*;
use ::fs_tool::prelude::traits::TryFrom;
use fs_tool::prelude::fs_tool::Calculator;
use std::io::BufRead;

pub const LIGAND_TABLE: &str = include_str!("resources/2019-12-29_lg.tsv");

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();
    let mut ligand_data = Vec::<LigandInfo>::new();
    let mut measures = vec![
        Measure {
            name: "CD8".to_string(),
            pos: vec![2, 3, 4, 5, 6, 9],
        },
        Measure {
            name: "NK".to_string(),
            pos: vec![2, 7, 8, 9],
        },
    ];

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .unwrap();

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
        let netmhcpan_summary = read_netmhcpan(f)?;

        if opt.drop_default_measures {
            measures = Vec::<Measure>::new();
        }

        if let Some(custom_measures) = opt.measures {
            custom_measures
                .iter()
                .filter_map(|measure| measure.parse::<Measure>().ok())
                .for_each(|measure| measures.push(measure))
        }

        let mut calculations = Calculator::new(&netmhcpan_summary, measures);
        calculations.process_measures();

        if let Some(output) = opt.output {
            let mut output = File::create(output)?;
            calculations.write_calculations(&mut output);
        } else {
            let mut output = io::stdout();
            calculations.write_calculations(&mut output);
        }
    }

    Ok(())
}
