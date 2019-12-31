use crate::prelude::fs_tool::*;
use std::env::var;
use structopt::StructOpt;

const RANK_TAG: &str = "# Rank";
const NN_TAG: &str = "HLA-";

#[derive(StructOpt, Debug)]
#[structopt(name = "fs-tool")]
pub struct Opt {
    #[structopt(long)]
    pub update_ligand_groups: bool,
    #[structopt(short, long, default_value = "4")]
    pub threads: usize,
}

/* Need to deal nom error */
fn read_netmhcpan<T: AsRef<str>>(input: T) -> Result<(), Box<dyn std::error::Error>> {
    let iter = input
        .as_ref()
        .lines()
        .filter(|line| !line.is_empty())
        .map(str::trim);
    let mut netmhcpan_summary = NetMHCpanSummary::new();

    iter.for_each(|line| {
        if line.starts_with(NN_TAG) {
            let (_, nn) = get_nn(line).unwrap();
            netmhcpan_summary.add_hla(nn);
        }

        if !netmhcpan_summary.is_threshold_set() && line.starts_with(RANK_TAG) {
            let (_, (mut rank_type, threshold)) = get_rank_threshold(line).unwrap();
            netmhcpan_summary.add_threshold(rank_type, threshold);
        }

        /* bytes version */
        if line.as_bytes()[0].is_ascii_digit() {
            if let Ok((_, (pep_info, variant_info, binding_info))) = process_netmhcpan_record(line)
            {
                netmhcpan_summary.add_sequence(
                    pep_info.4.clone(),
                    pep_info.0.clone(),
                    pep_info.1.clone(),
                );
                let peptide = Peptide::from((pep_info, variant_info));
                netmhcpan_summary.add_peptide(peptide);
            }
        }

        /* ascii compatible */
        //if line.chars().take(1).all(false, |mut a,c| {a = c.is_ascii_digit(); a}) {}

        /* Could be done on bytes and this one works only on acii */
        //            if let Some(_) = line.chars().next().filter(|c| c.is_dec_digit()) {}
    });

    dbg!(netmhcpan_summary);
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::reader::read_netmhcpan;

    const netmhcpan: &str =
            "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\
            \n\
            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000\n\
            -----------------------------------------------------------------------------------\n\
              Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel\n\
            -----------------------------------------------------------------------------------\n\
                1  HLA-A*03:01     TPQDLNTMLNT  TPQDTMLNT  0  4  2  0  0  TPQDLNTMLNT     Gag_180_209 0.0000160 83.3333\n\
                2  HLA-A*03:01     PQDLNTMLNTV  PQDLNTMLV  0  8  2  0  0  PQDLNTMLNTV     Gag_180_209 0.0000120 87.0000\n\
                3  HLA-A*03:01     QDLNTMLNTVG  QDLNTNTVG  0  5  2  0  0  QDLNTMLNTVG     Gag_180_209 0.0000040 96.0000\
            ";

    #[test]
    fn nn_line() {
        read_netmhcpan(netmhcpan).unwrap();
    }
}
