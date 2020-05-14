// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
// #![allow(dead_code, unused_variables)]

///  Contains errors associated with HLA/KIR nomenclature parsing, errors relating to parsing
/// the IPD/EBI website HTML and errors related to I/O
pub mod error;

/// Path to default kir ligand table for HLA ClassI allleles
pub const LIGAND_MAP_DEF: &str = include_str!("resources/allele_motifs.tsv");

/// Contains representations associated with MHC/HLA representations and nomenclature
pub mod mhc {
    /// MHC ClassI and ClassII related
    pub mod hla;
}

/// Contains representations associated with KIR/LILRB representations and nomenclature
pub mod ig_like {
    /// KIR related
    pub mod kir;
    /// Obtaining and representing KIR ligand information
    pub mod kir_ligand;
}
