use crate::prelude::collections::{HashMap, HashSet};
/// # TODO
/// 1. Implement conversion from Vec<&str>
/// 2. Add addition of Ligand Info using some sort of matching
/// 3. Implement comparison between/Retreiving motifs as Vec of strings
use crate::prelude::fs_tool::{BindLevel, Peptide, PeptideIdentity, Proteome, HLA};
use crate::prelude::traits::FromStr;

#[derive(Debug)]
pub struct NetMHCpanRecord<'a> {
    hla: &'a HLA,
    peptide: &'a Peptide<'a>,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

#[derive(Debug)]
pub struct NetMHCpanSummary<'a> {
    alleles: HashSet<HLA>,
    proteome: Proteome,
    peptides: HashMap<PeptideIdentity<'a>, Peptide<'a>>,
}
impl<'a> NetMHCpanRecord<'a> {
    fn new(
        hla: &'a HLA,
        peptide: &'a Peptide<'a>,
        bind_data: (f32, Option<f32>, f32, BindLevel),
    ) -> Self {
        let (score, aff, rank, bind_level) = bind_data;

        Self {
            hla,
            peptide,
            score,
            aff,
            rank,
            bind_level,
        }
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_create_netmhcpan_record() {}
}
