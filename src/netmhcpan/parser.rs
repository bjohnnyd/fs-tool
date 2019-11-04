use nom::{
    branch, bytes, character, combinator, error::ErrorKind, multi, sequence, Err::Error, IResult,
};

use super::*;
use crate::hla::mhc::HLA;
use nom::Err;

pub fn create_hla(hla_name: &str) -> () {
    if hla_name.contains(":") {
        let mut hla_information = hla_name.split(":");
        for val in hla_information {
            println!("{:?}", val);
        }
        ()
        //            let (
        //                gene_name,
        //                allele_group,
        //                hla_protein,
        //                cds_synonymous_sub,
        //                non_coding_difference,
        //                expression_change,
        //                ligand_group,
        //                mhc_class,
        //            ) = hla_information.unwrap();
        //            gene_name
    }
}

impl std::fmt::Debug for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <ParseError as std::fmt::Display>::fmt(self, f)
    }
}

impl std::error::Error for ParseError {}

//#[cfg(test)]
//mod test {
//    use super::*;
//    use crate::parsers::parsers::sanitize_hla_name;
//
//    #[test]
//    fn test_string_split() {
//        let with_delim_split: Vec<&str> = "A*01:101".split(":").collect();
//        let without_delim_split: Vec<&str> = "A*01101".split(":").collect();
//        assert_eq!(with_delim_split, vec!["A*01", "101"]);
//        assert_eq!(without_delim_split, vec!["A*01101"]);
//    }
//
//    #[test]
//    fn test_create_hla() {
//        assert_eq!(parsers::create_hla("HLA-A*01:101"), ());
//    }
//
//    #[test]
//    fn test_sanitize() {
//        assert_eq!(parsers::sanitize_hla_name("HLA-A*01:101"), "A01:101");
//    }
//}
