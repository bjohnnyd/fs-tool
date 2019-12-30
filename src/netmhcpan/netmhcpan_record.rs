/// # TODO
/// 1. Implement conversion from Vec<&str>
/// 2. Add addition of Ligand Info using some sort of matching
/// 3. Implement comparison between/Retreiving motifs as Vec of strings
use crate::prelude::fs_tool::{BindLevel, Peptide, HLA};
use crate::prelude::traits::FromStr;

#[derive(Debug)]
pub struct NetMHCpanRecord<'a> {
    hla: &'a HLA,
    peptide: Peptide<'a>,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

impl<'a> NetMHCpanRecord<'a> {
    fn new(
        hla: &'a HLA,
        peptide: Peptide<'a>,
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
