#[derive(Debug)]
pub struct Deletion(usize, usize);
#[derive(Debug)]
pub struct Insertion(usize, usize);

#[derive(Debug)]
pub enum BindLevel {
    SB,
    WB,
    NB,
}

#[derive(Debug)]
pub struct Peptide<'a> {
    pos: usize,
    length: usize,
    core: &'a [u8; 9],
    interaction_core: &'a [u8],
    core_start_offset: usize,
    deletion: Deletion,
    insertion: Insertion,
    protein: &'a Protein,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

#[derive(Debug)]
pub struct Protein {
    identity: String,
    sequence: String,
}

#[cfg(test)]
mod tests {}
