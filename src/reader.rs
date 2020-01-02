use crate::prelude::fs_tool::*;
use crate::prelude::io::{BufReader, Cursor, File, PathBuf, Read};
use std::io::{self, BufRead, Write};
use structopt::StructOpt;

const RANK_TAG: &str = "# Rank";
const NN_TAG: &str = "HLA-";

#[derive(StructOpt, Debug)]
#[structopt(
    name = "fs-tool",
    about = "Calculates fraction of shared peptides between HLA alleles based on NetMHCpan predictions"
)]
pub struct Opt {
    #[structopt(long)]
    pub update_ligand_groups: bool,

    #[structopt(short, long, default_value = "4")]
    pub threads: usize,

    #[structopt(short, long, parse(from_os_str))]
    pub netmhcpan: Option<PathBuf>,

    #[structopt(short, long, parse(from_os_str))]
    pub output: Option<PathBuf>,

    #[structopt(short, long)]
    pub measures: Option<Vec<String>>,

    #[structopt(long)]
    pub drop_default_measures: bool,
}

impl Opt {
    pub fn get_output(&self) -> Result<Box<dyn Write>, std::io::Error> {
        match &self.output {
            Some(path) => File::create(path).map(|f| Box::new(f) as Box<dyn Write>),
            None => Ok(Box::new(io::stdout())),
        }
    }
}

pub trait ToRead {
    type ToRead: Read;
    fn to_read(self) -> Self::ToRead;
}

impl<'a, T> ToRead for &'a T
where
    T: AsRef<str>,
{
    type ToRead = Cursor<&'a [u8]>;

    fn to_read(self) -> Self::ToRead {
        Cursor::new(self.as_ref().as_bytes())
    }
}

impl ToRead for File {
    type ToRead = File;

    fn to_read(self) -> File {
        self
    }
}

/* Need to deal nom error */
pub fn read_netmhcpan<T: ToRead>(input: T) -> Result<NetMHCpanSummary, Box<dyn std::error::Error>> {
    let reader = BufReader::new(input.to_read());
    let iter = reader
        .lines()
        .filter_map(|line| line.ok())
        .filter(|line| !line.is_empty());

    let mut netmhcpan_summary = NetMHCpanSummary::new();

    iter.for_each(|mut line| {
        line = line.trim().to_string();
        if line.starts_with(NN_TAG) {
            let (_, nn) = get_nn(&line).unwrap();
            netmhcpan_summary.add_hla(nn);
        }

        if !netmhcpan_summary.is_threshold_set() && line.starts_with(RANK_TAG) {
            let (_, (mut rank_type, threshold)) = get_rank_threshold(&line).unwrap();
            netmhcpan_summary.add_threshold(rank_type, threshold);
        }

        /* bytes version */
        if line.as_bytes()[0].is_ascii_digit() {
            if let Ok((_, (pep_info, variant_info, binding_info))) = process_netmhcpan_record(&line)
            {
                if netmhcpan_summary.alleles.len() <= 1 {
                    netmhcpan_summary.add_sequence(pep_info.4, pep_info.0, pep_info.1);
                    let peptide = Peptide::from((pep_info, variant_info));
                    netmhcpan_summary.add_peptide(peptide);
                }

                /* Need to deal with this error */
                let peptide_identity = PeptideIdentity::from(&pep_info);

                netmhcpan_summary.insert_hla_record(peptide_identity, binding_info);
            }
        }
    });

    Ok(netmhcpan_summary)
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
                3  HLA-A*03:01     QDLNTMLNTVG  QDLNTNTVG  0  5  2  0  0  QDLNTMLNTVG     Gag_180_209 0.0000040 96.0000\n\
                4  HLA-A*03:01     DLNTMLNTVGG  DLNTNTVGG  0  4  2  0  0  DLNTMLNTVGG     Gag_180_209 0.0000040 96.0000\n\
                5  HLA-A*03:01     LNTMLNTVGGH  LMLNTVGGH  0  1  2  0  0  LNTMLNTVGGH     Gag_180_209 0.0001090 48.8333\n\
                6  HLA-A*03:01     NTMLNTVGGHQ  NTMTVGGHQ  0  3  2  0  0  NTMLNTVGGHQ     Gag_180_209 0.0001260 46.1429\n\
                7  HLA-A*03:01     TMLNTVGGHQA  TMLNGGHQA  0  4  2  0  0  TMLNTVGGHQA     Gag_180_209 0.0001260 46.1429\n\
                8  HLA-A*03:01     MLNTVGGHQAA  MLNGGHQAA  0  3  2  0  0  MLNTVGGHQAA     Gag_180_209 0.0002300 36.4000\n\
                9  HLA-A*03:01     LNTVGGHQAAM  NTVGGHQAM  1  7  1  0  0   NTVGGHQAAM     Gag_180_209 0.0000530 62.5000\n\
               10  HLA-A*03:01     NTVGGHQAAMQ  NTVGGAAMQ  0  5  2  0  0  NTVGGHQAAMQ     Gag_180_209 0.0001420 44.1250\n\
               11  HLA-A*03:01     TVGGHQAAMQM  TVHQAAMQM  0  2  2  0  0  TVGGHQAAMQM     Gag_180_209 0.0004120 28.5000\n\
               12  HLA-A*03:01     VGGHQAAMQML  VGGAAMQML  0  3  2  0  0  VGGHQAAMQML     Gag_180_209 0.0000120 87.0000\n\
               13  HLA-A*03:01     GGHQAAMQMLK  GQAAMQMLK  0  1  2  0  0  GGHQAAMQMLK     Gag_180_209 0.0313010  3.8215\n\
            ";

    #[test]
    fn nn_line() {
        read_netmhcpan(&netmhcpan).unwrap();
    }
}
