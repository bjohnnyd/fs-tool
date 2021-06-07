mod factory;

pub use immunoprot::ig_like::kir_ligand::KirLigandMap;
use once_cell::sync::Lazy;

use crate::calc::Measure;

pub static KIR_MAP: Lazy<KirLigandMap> = Lazy::new(|| KirLigandMap::init().unwrap());
pub static TEST_MEASURES: Lazy<[Measure; 2]> = Lazy::new(|| {
    [
        Measure {
            name: "TCR".to_string(),
            motif_pos: vec![2, 3, 4, 5, 6, 9],
        },
        Measure {
            name: "KIR".to_string(),
            motif_pos: vec![2, 7, 8, 9],
        },
    ]
});
