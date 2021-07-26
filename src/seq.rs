/// Represents a protein sequence
struct SeqProtein {
    id: String,
    seq: Vec<u8>,
}

/// Peptide sequence that points to a protein and has length
struct SeqPeptide<'a> {
    protein: &'a SeqProtein,
    start: usize,
    end: usize,
}

impl SeqProtein {
    pub fn new(id: String, seq: Vec<u8>) -> Self {
        Self { id, seq }
    }

    pub fn id(&self) -> &str {
        &self.id
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
}

impl<'a> SeqPeptide<'a> {
    fn new(protein: &'a SeqProtein, start: usize, end: usize) -> Self {
        Self {
            protein,
            start,
            end,
        }
    }

    fn seq(&self) -> &[u8] {
        &self.protein.seq()[self.start..self.end]
    }

    fn motif(&self, pos: &[usize]) -> Vec<&u8> {
        todo!()
    }
}
