use crate::prelude::collections::{HashMap, HashSet};
/// # TODO
/// 1. Implement conversion from Vec<&str>
/// 2. Add addition of Ligand Info using some sort of matching
/// 3. Implement comparison between/Retreiving motifs as Vec of strings
use crate::prelude::fs_tool::{BindLevel, Peptide, PeptideIdentity, Proteome, HLA};
use crate::prelude::traits::FromStr;

pub struct StrongThreshold(f32);
pub struct WeakThreshold(f32);
pub struct Distance(f32);

#[derive(Debug)]
pub struct NearestNeighbour {
    index: HLA,
    distance: f32,
    nearest_neighbour: HLA,
}

impl PartialEq for NearestNeighbour {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.nearest_neighbour == other.nearest_neighbour
    }
}

impl std::hash::Hash for NearestNeighbour {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.nearest_neighbour.hash(state);
    }
}

impl Eq for NearestNeighbour {}

impl From<(HLA, f32, HLA)> for NearestNeighbour {
    fn from(t: (HLA, f32, HLA)) -> Self {
        Self {
            index: t.0,
            distance: t.1,
            nearest_neighbour: t.2,
        }
    }
}
#[derive(Debug)]
pub struct NetMHCpanRecord<'a> {
    hla: &'a HLA,
    peptide: &'a Peptide<'a>,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
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

#[derive(Debug)]
pub struct NetMHCpanSummary<'a> {
    alleles: HashSet<NearestNeighbour>,
    records: HashMap<&'a HLA, Vec<&'a Peptide<'a>>>,
    proteome: Proteome,
    peptides: HashMap<PeptideIdentity<'a>, Peptide<'a>>,
    pub weak_threshold: Option<f32>,
    pub strong_threshold: Option<f32>,
}

impl<'a> NetMHCpanSummary<'a> {
    pub(crate) fn new() -> Self {
        Self {
            alleles: HashSet::<NearestNeighbour>::new(),
            records: HashMap::<&HLA, Vec<&Peptide>>::new(),
            proteome: Proteome::new(),
            peptides: HashMap::<PeptideIdentity, Peptide>::new(),
            weak_threshold: None,
            strong_threshold: None,
        }
    }

    pub fn add_hla(&mut self, nn: NearestNeighbour) -> bool {
        self.alleles.insert(nn)
    }

    pub fn add_peptide(&mut self, id: PeptideIdentity<'a>, pep: Peptide<'a>) -> Option<Peptide> {
        self.peptides.insert(id, pep)
    }

    pub fn add_sequence(&mut self, id: &'a str, pos: usize, seq: &'a str) {
        self.proteome.add_peptide(id, pos, seq)
    }

    pub fn is_threshold_set(&self) -> bool {
        self.weak_threshold.is_some() && self.strong_threshold.is_some()
    }

    pub fn add_threshold(&mut self, rank_type: &str, threshold: f32) {
        if rank_type.trim() == "Strong" {
            self.strong_threshold = Some(threshold);
        } else if rank_type.trim() == "Weak" {
            self.weak_threshold = Some(threshold)
        }
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_create_netmhcpan_record() {}
}
