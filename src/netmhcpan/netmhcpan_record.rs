use crate::netmhcpan::parser::{BindingInfo, PepInfo};
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
pub struct NetMHCpanRecord {
    peptide_identity: PeptideIdentity,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

impl NetMHCpanRecord {
    fn new(
        peptide_identity: PeptideIdentity,
        bind_data: (f32, Option<f32>, f32, BindLevel),
    ) -> Self {
        let (score, aff, rank, bind_level) = bind_data;

        Self {
            peptide_identity,
            score,
            aff,
            rank,
            bind_level,
        }
    }
}

#[derive(Debug)]
pub struct NetMHCpanSummary {
    pub alleles: HashSet<NearestNeighbour>,
    records: HashMap<HLA, Vec<NetMHCpanRecord>>,
    proteome: Proteome,
    peptides: HashMap<PeptideIdentity, Peptide>,
    pub weak_threshold: Option<f32>,
    pub strong_threshold: Option<f32>,
}

impl NetMHCpanSummary {
    pub(crate) fn new() -> Self {
        Self {
            alleles: HashSet::<NearestNeighbour>::new(),
            records: HashMap::<HLA, Vec<NetMHCpanRecord>>::new(),
            proteome: Proteome::new(),
            peptides: HashMap::<PeptideIdentity, Peptide>::new(),
            weak_threshold: None,
            strong_threshold: None,
        }
    }

    pub fn add_hla(&mut self, nn: NearestNeighbour) -> bool {
        self.alleles.insert(nn)
    }

    pub fn add_peptide(&mut self, pep: Peptide) -> Option<Peptide> {
        let peptide_identity = PeptideIdentity::from(&pep);
        self.peptides.insert(peptide_identity, pep)
    }

    pub fn add_sequence<T: AsRef<str>>(&mut self, id: T, pos: usize, seq: T) {
        self.proteome
            .add_peptide(id.as_ref().to_string(), pos, seq.as_ref().to_string())
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

    /* Need to deal with error */
    pub fn insert_hla_record(
        &mut self,
        peptide_identity: PeptideIdentity,
        BindingInfo(hla_id, score, aff, rank, bind_level): BindingInfo,
    ) {
        let hla = hla_id.parse::<HLA>().expect("could not parse hla");

        let netmhcpan_record =
            NetMHCpanRecord::new(peptide_identity, (score, aff, rank, bind_level));
        let records = self.records.entry(hla).or_insert(Vec::new());
        records.push(netmhcpan_record);
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn test_create_netmhcpan_record() {}
}
