#[allow(unused)]
pub(in crate::tests) mod test_helpers;

// use assert_cmd::prelude::*;
// use predicates::str::contains;
// use std::process::Command;

// #[test]
// No arguments should mention help and required field
// fn cli_no_args_help() {
// Command::cargo_bin("fs-tool")
// .unwrap()
// .assert()
// .stderr(contains("help"))
// .stderr(contains("--binding_predictions"))
// .stderr(contains("--output"));
// }

// #[test]
// Should throw an error for incorrect binding prediction file
// fn cli_wrong_binding_predictions() {
// Command::cargo_bin("fs-tool")
// .unwrap()
// .args(&[
// "-b",
// "tests/input/binding_predictions/large_cohort_batch1.txt",
// "-o",
// "tests/output/alleles",
// ])
// .assert()
// .stderr(contains("No such file"));
// }
