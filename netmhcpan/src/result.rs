// output info at http://www.cbs.dtu.dk/services/NetMHC/output.php

use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::ops::Range;

use crate::error::Error;

use immunoprot::mhc::hla::ClassI;
use log::warn;

pub const WEAK_TRESHOLD: f32 = 2.0;
pub const STRONG_THRESHOLD: f32 = 0.5;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum RankThreshold {
    Strong(f32),
    Weak(f32),
}

#[derive(Debug, Clone, PartialOrd)]
pub struct NearestNeighbour {
    pub(crate) index: ClassI,
    pub(crate) distance: f32,
    pub(crate) nn: ClassI,
}

impl NearestNeighbour {
    pub fn new(index: ClassI, distance: f32, nn: ClassI) -> Self {
        Self {
            index,
            distance,
            nn,
        }
    }

    pub fn info(&self) -> (&ClassI, f32, &ClassI) {
        (&self.index, self.distance, &self.nn)
    }
}

//TODO: FIX need to use `cmp`
impl PartialEq for NearestNeighbour {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.nn == other.nn
    }
}

impl Eq for NearestNeighbour {}

impl Hash for NearestNeighbour {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.nn.hash(state);
    }
}

// The core is a 9mer always used for alignment and identification
#[derive(Debug, Clone, Eq)]
pub struct Peptide {
    // 0-based in documents but never represented as such in the results so need to `-1` shift before
    pub(crate) pos: usize,
    pub(crate) seq: String,
    pub(crate) identity: String,
    // Represents the epitope. This might be possible to deduce??
    pub(crate) icore: String,
    pub(crate) offset: usize,
    pub(crate) gap: Range<usize>,
    pub(crate) ins: Range<usize>,
}

impl Peptide {
    pub fn new(
        pos: usize,
        seq: String,
        identity: String,
        icore: String,
        alignment_mods: &[usize],
    ) -> Self {
        let offset = alignment_mods[0];
        let gap = Range {
            start: alignment_mods[1],
            end: alignment_mods[1..3].iter().sum(),
        };
        let ins = Range {
            start: alignment_mods[3],
            end: alignment_mods[3..5].iter().sum(),
        };

        Self {
            pos,
            seq,
            identity,
            icore,
            offset,
            gap,
            ins,
        }
    }

    pub fn seq(&self) -> &str {
        &self.seq
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn gap_region(&self) -> &Range<usize> {
        &self.gap
    }

    pub fn ins_region(&self) -> &Range<usize> {
        &self.ins
    }

    pub fn sequence(&self) -> &str {
        &self.seq
    }

    pub fn protein(&self) -> &str {
        &self.identity
    }

    // TODO: Need to see what happens or if there are cases with gap + ins
    /// Converts a peptide sequence to core representation.  Undefined behaviour with cases where there
    /// are both gaps and insertions in the alignment
    pub fn core(&self) -> String {
        if self.gap_region().len() != 0 {
            self.seq()
                .chars()
                .enumerate()
                .filter(|(i, _)| !self.gap.contains(i))
                .map(|(_, c)| c)
                .collect()
        } else if self.ins.len() != 0 {
            self.seq().chars().enumerate().map(|(_, c)| c).collect()
        } else {
            self.seq().to_string()
        }
    }

    pub fn icore(&self) -> &str {
        self.icore.as_ref()
    }

    pub fn sequence_motif(&self, aa_pos: &[usize]) -> String {
        self.seq()
            .chars()
            .enumerate()
            .filter(|(i, _)| aa_pos.contains(i))
            .map(|(_, c)| c)
            .collect()
    }

    /// Returns the differences from original sequence in the format
    /// `(offset, gap range, insertion range)`
    pub fn aa_diff(&self) -> (usize, &Range<usize>, &Range<usize>) {
        (self.offset, &self.gap, &self.ins)
    }
}

impl PartialEq for Peptide {
    fn eq(&self, other: &Self) -> bool {
        self.pos == other.pos
            && self.len() == other.len()
            && self.identity == other.identity
            && self.seq() == other.seq()
    }
}

// NOTE: Might not be needed anymore due to PartialEq
impl Hash for Peptide {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.pos.hash(state);
        self.len().hash(state);
        self.seq().hash(state);
        self.identity.hash(state);
    }
}

#[derive(Debug, Clone, Eq, Ord, PartialOrd)]
pub struct Protein {
    identity: String,
    sequence: String,
}

impl Protein {
    pub fn new<T>(identity: T) -> Self
    where
        T: AsRef<str>,
    {
        let sequence = String::from("");
        let identity = String::from(identity.as_ref());

        Self { sequence, identity }
    }

    pub fn identity(&self) -> &str {
        &self.identity
    }

    pub fn seq(&self) -> &str {
        &self.sequence
    }

    /// Adds a peptide sequence to the current protein at position specified.
    /// 0-based so if passing NetMHCpan 4.0 `Pos` column `-1` needs to be subtracted.
    /// The position + sequence length can only be larger than 1 for it to be added.
    /// If the sequence + position are within current sequence length the sequence is overwritten
    pub fn add_sequence_at_pos<T>(&mut self, pos: usize, sequence: T) -> Result<(), Error<()>>
    where
        T: AsRef<str>,
    {
        let sequence = sequence.as_ref();
        let region = Range {
            start: pos,
            end: pos + sequence.len(),
        };

        if region.start > self.sequence.len() {
            return Err(Error::ProteinTooShort(pos, self.sequence.len()));
        } else if region.end > self.sequence.len() {
            self.sequence.replace_range(region.start.., sequence)
        } else {
            self.sequence.replace_range(region, sequence)
        }

        Ok(())
    }

    /// Returns protein sequence at specified positions (0 based and right open-ended)
    ///
    pub fn sequence(&self, start: usize, end: usize) -> Result<&str, Error<()>> {
        if end > self.sequence.len() {
            Err(Error::ProteinTooShort(end, self.sequence.len()))
        } else {
            Ok(&self.sequence[start..end])
        }
    }
}

impl PartialEq for Protein {
    fn eq(&self, other: &Self) -> bool {
        self.identity == other.identity && self.sequence == other.sequence
    }
}
impl Hash for Protein {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.identity.hash(state);
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BindingInfo {
    pub(crate) peptide: Peptide,
    pub(crate) score: f32,
    pub(crate) affinity: Option<f32>,
    pub(crate) rank: f32,
}

impl BindingInfo {
    pub fn motif(&self, aa_pos: &[usize]) -> String {
        self.peptide.sequence_motif(aa_pos)
    }

    pub fn peptide(&self) -> &Peptide {
        &self.peptide
    }

    pub fn seq(&self) -> &str {
        self.peptide.seq()
    }
    pub fn len(&self) -> usize {
        self.peptide.len()
    }

    pub fn rank(&self) -> f32 {
        self.rank
    }
}

#[derive(Debug, Clone)]
pub struct BindingData {
    pub(crate) alleles: HashSet<NearestNeighbour>,
    pub(crate) allele_binding: HashMap<ClassI, Vec<BindingInfo>>,
    pub(crate) proteome: HashMap<String, Protein>,
    pub(crate) peptides: HashSet<Peptide>,
    pub(crate) weak_threshold: Option<f32>,
    pub(crate) strong_threshold: Option<f32>,
}

impl Default for BindingData {
    fn default() -> Self {
        Self {
            alleles: HashSet::<NearestNeighbour>::new(),
            allele_binding: HashMap::<ClassI, Vec<BindingInfo>>::new(),
            proteome: HashMap::<String, Protein>::new(),
            peptides: HashSet::<Peptide>::new(),
            weak_threshold: None,
            strong_threshold: None,
        }
    }
}

impl BindingData {
    pub fn new() -> Self {
        BindingData::default()
    }

    pub fn get_binding_info(&self, allele: &ClassI) -> Option<&Vec<BindingInfo>> {
        let binding_info = self.allele_binding.get(allele);
        if binding_info.is_none() {
            warn!("{} has no associated binding data", allele);
        };
        binding_info
    }

    pub fn list_alleles(&self) -> Vec<&ClassI> {
        self.alleles.iter().map(|nn| &nn.index).collect()
    }

    pub fn list_nn(&self) -> &HashSet<NearestNeighbour> {
        &self.alleles
    }

    pub fn proteins(&self) -> Vec<String> {
        self.proteome
            .keys()
            .map(|identity| identity.to_string())
            .collect::<Vec<String>>()
    }

    pub fn pep_lengths(&self) -> Vec<usize> {
        let mut pep_lengths = self
            .peptides
            .iter()
            .map(|pep| pep.len())
            .collect::<Vec<usize>>();
        pep_lengths.sort();
        pep_lengths.dedup();
        pep_lengths
    }

    pub fn strong_threshold(&self) -> f32 {
        if let Some(threshold) = self.strong_threshold {
            threshold
        } else {
            warn!(
                "No binding threshold found in NetMHCpan result using default of {}",
                STRONG_THRESHOLD
            );
            STRONG_THRESHOLD
        }
    }

    pub fn weak_threshold(&self) -> f32 {
        if let Some(threshold) = self.weak_threshold {
            threshold
        } else {
            warn!(
                "No binding threshold found in NetMHCpan result using default of {}",
                WEAK_TRESHOLD
            );
            WEAK_TRESHOLD
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::result::{Peptide, Protein};

    #[test]
    fn test_peptide_core() {
        let identity = String::from("Gag_180_209");
        let seq1 = String::from("GHQAAMQMLK");
        let seq2 = String::from("HQAAMQMLK");
        let mut protein = Protein::new(identity.clone());
        protein.add_sequence_at_pos(0, seq1.clone()).unwrap();
        protein.add_sequence_at_pos(1, seq2.clone()).unwrap();

        let alignments_mode_none = [0; 5];
        let alignments_mode_one = [0, 1, 1, 0, 0];
        let pep_identical = Peptide::new(
            1,
            seq2.clone(),
            identity.clone(),
            seq1.clone(),
            &alignments_mode_none,
        );
        let pep_diff = Peptide::new(
            0,
            seq1.clone(),
            identity.clone(),
            seq2.clone(),
            &alignments_mode_one,
        );

        let expected_core_diff = String::from("GQAAMQMLK");

        assert_eq!(pep_identical.core(), pep_identical.sequence().to_string());
        assert_ne!(pep_diff.core(), pep_diff.sequence().to_string());
        assert_eq!(pep_diff.core(), expected_core_diff);
    }
    #[test]
    fn test_add_seq_to_protein() {
        let mut protein = Protein::new("Gag");
        let expected = Protein {
            identity: "Gag".to_string(),
            sequence: "ABCDEF".to_string(),
        };
        protein.add_sequence_at_pos(0, "ABCDEF").unwrap();
        assert_eq!(protein, expected);

        protein.add_sequence_at_pos(2, "DEFG").unwrap();
        let expected = Protein {
            identity: "Gag".to_string(),
            sequence: "ABDEFG".to_string(),
        };
        assert_eq!(protein, expected);
    }
}
