extern crate nom;
extern crate reqwest;
extern crate select;

mod hla;
mod ligands_groups;
mod parsers;

use hla::mhc::HLA;

fn main() -> std::result::Result<(), std::boxed::Box<dyn std::error::Error>> {
    //    hla = HLA::new(String::from("HLA-A*01:02"));
    Ok(())
}
