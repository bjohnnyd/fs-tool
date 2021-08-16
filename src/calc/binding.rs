use std::collections::HashMap;

use crate::calc::motif::Measure;
use fxhash::FxBuildHasher;
use immunoprot::mhc::hla::ClassI;

pub type MotifsBound = HashMap<String, usize, FxBuildHasher>;

#[derive(Debug)]
pub struct BindingMhc<'a> {
    hla: &'a ClassI,
    binding_data: HashMap<&'a Measure, MotifsBound>,
}

impl<'a> BindingMhc<'a> {
    pub fn new(hla: &'a ClassI) -> Self {
        Self {
            hla,
            binding_data: HashMap::new(),
        }
    }

    pub fn add_measure(&mut self, measure: &'a Measure, motifs: MotifsBound) {
        self.binding_data.insert(measure, motifs);
    }
}

#[cfg(test)]
mod test {

    use super::*;

    const cd8_indices: &str = "2,3,4,5,6,7,9";
    const kir_indices: &str = "2,7,8,9";

    const orig_pep1: &str = "TPQDLNTML";
    const orig_pep2: &str = "PQDLNTMLN";

    #[test]
    fn create_binding_hla() {
        let hla = ClassI::default();
        let binding_mhc = BindingMhc::new(&hla);
        assert_eq!(binding_mhc.hla, &hla);
    }
}
