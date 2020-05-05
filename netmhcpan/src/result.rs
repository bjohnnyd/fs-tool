use std::collections::{HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::ops::Range;

use crate::error::Error;
use immunoprot::mhc::hla::ClassI;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum RankThreshold {
    Strong(f32),
    Weak(f32),
}

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub struct NearestNeighbour {
    index: ClassI,
    distance: f32,
    nn: ClassI,
}

impl NearestNeighbour {
    pub fn new(index: ClassI, distance: f32, nn: ClassI) -> Self {
        Self {
            index,
            distance,
            nn,
        }
    }
}

// The core is a 9mer always used for alignment and identification
#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct Peptide<'a> {
    // 0-based in documents but never represented as such in the results so need to `-1` shift before
    pos: usize,
    len: usize,
    protein: &'a Protein,
    // Represents the epitope. This might be possible to deduce??
    icore: String,
    // predicts an n-terminal protrusion if > 0
    offset: usize,
    gap_start: usize,
    gap_len: usize,
    ins_start: usize,
    ins_len: usize,
}

impl<'a> Peptide<'a> {
    fn new(
        pos: usize,
        len: usize,
        protein: &'a Protein,
        icore: String,
        offset: usize,
        gap_start: usize,
        gap_len: usize,
        ins_start: usize,
        ins_len: usize,
    ) -> Self {
        Self {
            pos,
            len,
            protein,
            icore,
            offset,
            gap_start,
            gap_len,
            ins_start,
            ins_len,
        }
    }

    pub fn sequence(&self) -> &str {
        &self.protein.sequence[self.pos..(self.pos + self.len)]
    }

    pub fn protein(&self) -> &Protein {
        &*self.protein
    }

    // TODO: Are there any cases where there is a gap and insertion?
    pub fn core(&self) -> String {
        let mut seq = self.sequence().to_string();

        if self.gap_len != 0 {
            let gap = Range {
                start: self.gap_start,
                end: self.gap_len + 1,
            };

            seq = seq
                .chars()
                .enumerate()
                .filter(|(i, c)| !gap.contains(i))
                .map(|(_, c)| c)
                .collect();
        }

        if self.ins_len != 0 {
            let ins = Range {
                start: self.ins_start,
                end: self.ins_len + 1,
            };

            seq = seq
                .chars()
                .enumerate()
                .map(|(i, c)| if ins.contains(&i) { '-' } else { c })
                .collect();
        }
        seq
    }

    pub fn icore(&self) -> &str {
        self.icore.as_ref()
    }

    pub fn sequence_motif(&self, aa_pos: &[usize]) -> String {
        self.sequence()
            .chars()
            .enumerate()
            .filter(|(i, c)| aa_pos.contains(i))
            .map(|(_, c)| c)
            .collect()
    }

    /// Returns the differences from original sequence in the format
    /// `(offset, gap start, gap length, insertion start, insertion length)`
    pub fn aa_diff(&self) -> (usize, usize, usize, usize, usize) {
        (
            self.offset,
            self.gap_start,
            self.gap_len,
            self.ins_start,
            self.ins_len,
        )
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct Protein {
    identity: String,
    sequence: String,
}

impl Protein {
    fn new<T>(identity: T) -> Self
    where
        T: AsRef<str>,
    {
        let sequence = String::from("");
        let identity = String::from(identity.as_ref());

        Self { sequence, identity }
    }

    /// Adds a peptide sequence to the current protein at position specified.
    /// 0-based so if passing NetMHCpan 4.0 `Pos` column `-1` needs to be subtracted.
    /// The position + sequence length can only be larger than 1 for it to be added.
    /// If the sequence + position are within current sequence length the sequence is overwritten
    fn add_sequence_at_pos<T>(&mut self, pos: usize, sequence: T) -> Result<(), Error>
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
        } else {
            if region.end > self.sequence.len() {
                self.sequence.replace_range(region.start.., sequence)
            } else {
                self.sequence.replace_range(region, sequence)
            }
        }
        Ok(())
    }
}

impl Hash for Protein {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.identity.hash(state);
    }
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub struct Proteome<'a> {
    pub(crate) proteins: HashSet<Protein>,
    pub(crate) peptides: HashMap<Protein, Vec<Peptide<'a>>>,
}
#[derive(Debug, Clone, PartialEq)]
pub struct BindingData<'a> {
    pub(crate) alleles: HashSet<ClassI>,
    pub(crate) proteome: Proteome<'a>,
    pub(crate) weak_threshold: Option<f32>,
    pub(crate) strong_threshold: Option<f32>,
}

#[cfg(test)]
mod tests {
    use crate::result::{Peptide, Protein};

    #[test]
    fn test_peptide_core() {
        let mut protein = Protein::new("Gag_180_209");
        protein.add_sequence_at_pos(0, "GHQAAMQMLK");
        protein.add_sequence_at_pos(1, "HQAAMQMLK");

        let pep_identical = Peptide::new(1, 9, &protein, String::from("HQAAMQMLK"), 0, 0, 0, 0, 0);
        let pep_diff = Peptide::new(0, 10, &protein, String::from("GHQAAMQMLK"), 0, 1, 1, 0, 0);

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
