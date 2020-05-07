// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
#![allow(dead_code, unused_variables)]
pub mod error;

pub mod mhc {
    pub mod hla;
}

pub mod ig_like {
    pub mod kir;
    pub mod kir_ligand;
    pub mod lilrb;
}

#[cfg(test)]
mod tests {}
