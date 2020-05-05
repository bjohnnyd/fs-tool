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
    // 0-based in documents but never represented as such in the results
    pos: usize,
    protein: &'a Protein,
    sequence: &'a [u8],
    // Represents the epitope. This might be possible to deduce??
    interaction_core: String,
    // predicts an n-terminal protrusion if > 0
    offset: usize,
    del_pos: usize,
    del_len: usize,
    ins_pos: usize,
    ins_len: usize,
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
    /// The position + sequence length can only be larger than 1 for it to be added.
    /// If the sequence + position are within current sequence length the sequence is overwritten
    fn add_sequence_at_pos<T>(&mut self, pos: usize, sequence: T) -> Result<(), Error>
    where
        T: AsRef<str>,
    {
        let sequence = sequence.as_ref();
        let region = Range {
            start: pos - 1,
            end: pos - 1 + sequence.len(),
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
    pub(crate) peptides: HashMap<Protein, Peptide<'a>>,
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
    use crate::result::Protein;

    #[test]
    fn test_add_seq_to_protein() {
        let mut protein = Protein::new("Gag");
        let expected = Protein {
            identity: "Gag".to_string(),
            sequence: "ABCDEF".to_string(),
        };
        protein.add_sequence_at_pos(1, "ABCDEF").unwrap();
        assert_eq!(protein, expected);

        protein.add_sequence_at_pos(3, "DEFG").unwrap();
        let expected = Protein {
            identity: "Gag".to_string(),
            sequence: "ABDEFG".to_string(),
        };
        assert_eq!(protein, expected);
    }
}
