use ::fs_tool::error::CouldNotOpenFile;
use ::fs_tool::prelude::external::{Opt, StructOpt};
use ::fs_tool::prelude::fs_tool::{
    get_ligand_table, parse_ligand_table, read_netmhcpan, LigandInfo, Measure, NetMHCpanSummary,
    HLA,
};
use ::fs_tool::prelude::io::*;
use ::fs_tool::prelude::snafu_error::ResultExt;
use ::fs_tool::prelude::traits::TryFrom;
use fs_tool::prelude::fs_tool::Calculator;

fn main() -> std::result::Result<(), ()> {
    match main_try() {
        Ok(()) => Ok(()),
        Err(e) => {
            eprintln!("Error: {}", e);
            Err(())
        }
    }
}

fn main_try() -> Result<(), Box<dyn std::error::Error>> {
    let mut opt = Opt::from_args();

    rayon::ThreadPoolBuilder::new()
        .num_threads(opt.threads)
        .build_global()
        .unwrap();

    let ligand_hla = opt.get_ligand_data();
    let mut output = opt.get_output()?;

    opt.set_measures();

    if let Some(netmhcout) = opt.netmhcpan {
        let f = File::open(netmhcout.clone()).context(CouldNotOpenFile { f_path: netmhcout })?;
        let netmhcpan_summary = read_netmhcpan(f)?;

        if let Some(measures) = opt.measures {
            let mut calculations =
                Calculator::new(&netmhcpan_summary, measures, opt.peptide_length);
            calculations.process_measures();
            calculations.write_calculations(&mut output)?;
        }
    }

    Ok(())
}
