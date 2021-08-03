mod factory;

use crate::io::reader::read_kir_motif_binding;
use immunoprot::ig_like::kir_ligand::KirLigandMap;
use once_cell::sync::Lazy;
use rand::prelude::*;
use std::collections::HashMap;

pub use crate::calc::get_bound_kirs;
pub use crate::calc::motif::Measure;
pub use immunoprot::ig_like::kir::Kir;
pub use immunoprot::ig_like::kir_ligand::LigandMotif;
pub use immunoprot::mhc::hla::ClassI;
pub use netmhcpan::result::{BindingInfo, Peptide};

/// Generic type for erros in tests where we don't care about the underlying error kind
pub type TestResult<T, E = Box<dyn std::error::Error>> = Result<T, E>;

/// Ensures that all tests are using the same kir ligand info
pub static KIR_MAP: Lazy<KirLigandMap> = Lazy::new(|| KirLigandMap::init().unwrap());

/// Set of measures to be used in tests with the last measure being a subset of the other 2
/// ensuring that it is always at least as large both of the other measures
pub static TEST_MEASURES: Lazy<[Measure; 3]> = Lazy::new(|| {
    [
        Measure::new(String::from("TCR"), "2,3,4,5,6,9").unwrap(),
        Measure::new(String::from("KIR"), "2,7,8,9").unwrap(),
        Measure::new(String::from("ANCHOR"), "2,9").unwrap(),
    ]
});
