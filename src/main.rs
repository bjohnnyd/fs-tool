mod hla;

use hla::mhc::HLA;

fn main() -> std::result::Result<(), std::boxed::Box<dyn std::error::Error>> {
    hla = HLA::new(String::from("HLA-A*01:02"));
}
