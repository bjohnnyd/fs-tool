/// # TODO
/// 1. Implement conversion from Vec<&str>
/// 2. Add addition of Ligand Info using some sort of matching
/// 3. Implement comparison between/Retreiving motifs as Vec of strings
use crate::prelude::fs_tool::{BindLevel, Peptide, HLA};

#[derive(Debug)]
pub struct NetMHCpanRecord<'a> {
    hla: &'a HLA,
    peptide: Peptide<'a>,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}
