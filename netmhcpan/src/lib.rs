// #![warn(missing_debug_implementations, rust_2018_idioms, missing_docs)]
//#![allow(dead_code, unused_variables, clippy::len_without_is_empty )]
#![allow(clippy::len_without_is_empty)]
pub mod error;
pub mod parser;
pub mod reader;
pub mod result;

pub const WEAK_TRESHOLD: f32 = 2.0;
pub const STRONG_THRESHOLD: f32 = 0.5;
