use std::collections::HashMap;

use crate::calc::motif::Measure;
use fxhash::{FxBuildHasher, FxHashMap};
use immunoprot::mhc::hla::ClassI;

#[derive(Debug)]
pub struct BindingMhc<'a> {
    hla: &'a ClassI,
    binding_data: HashMap<&'a Measure, HashMap<String, usize, FxBuildHasher>>,
}

impl<'a> BindingMhc<'a> {
    pub fn new(hla: &'a ClassI) -> Self {
        Self {
            hla,
            binding_data: HashMap::new(),
        }
    }

    pub fn add_measure(
        &mut self,
        measure: &'a Measure,
        motifs: HashMap<String, usize, FxBuildHasher>,
    ) {
        self.binding_data.insert(measure, motifs);
    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn create_binding_hla() {
        let hla = ClassI::default();
        let binding_mhc = BindingMhc::new(&hla);
        assert_eq!(binding_mhc.hla, &hla);
    }
}
