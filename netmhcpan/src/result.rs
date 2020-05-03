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
