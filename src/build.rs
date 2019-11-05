// build.rs

use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

fn main() {
    let data_dir = PathBuf::from("generated");
    fs::create_dir(data_dir);

    let netmhcpan_alleles = File::open("data/MHC_allele_names.txt");

    let mut f = File::create(data_dir.join("data.rs")).unwrap();


}
