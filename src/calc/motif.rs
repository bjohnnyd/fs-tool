use std::collections::HashSet;

use netmhcpan::result::Peptide;

// Multiple ways to exract and compare motifs but there are two defining parameters
// position and length. I want to go over peptides if unique than occurence can be set to 1 and
// otherwise not

/// Represents a sequence motif bound by
pub struct SequenceMotif<'a> {
    peptides: HashSet<&'a Peptide>,
    /// A boolean mapping to peptide lengths 8,9,10 and 11 so peptide length can be got as 8 + idx
    /// == TRUE
    peptide_lengths: [usize; 4],
    sequence: &'a str,
    unique: bool,
}

impl SequenceMotif<'_> {
    pub fn new(peptides: &Vec<Peptide>, sequence: &str) -> Self {
        todo!()
    }
    pub fn occurence_in_lengths(&self) -> Vec<&usize> {
        self.peptide_lengths
            .iter()
            .filter(|length| **length != 0usize)
            .collect()
    }
}
