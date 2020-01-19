mod calc;
mod cohort;
mod data;
pub mod error;
mod mhc;
mod netmhcpan;
mod reader;

pub mod prelude {

    pub mod io {
        pub use std::fs::{self, File};
        pub use std::io::{self, BufReader, Cursor, Read, Write};
        pub use std::path::PathBuf;
    }

    pub mod error {
        pub use crate::data::errors as data_errors;
        pub use crate::mhc::errors as mhc_errors;
        pub use crate::netmhcpan::errors as netmhcpan_errors;
        pub use snafu::{ensure, Backtrace, ErrorCompat, ResultExt, Snafu};

        pub use data_errors::RetrieveLigandError;
        pub use mhc_errors::HLAError;
        pub use netmhcpan_errors::*;
    }

    pub mod snafu_error {
        pub use snafu::{ensure, Backtrace, ErrorCompat, ResultExt, Snafu};
    }

    pub mod fs_tool {
        pub use crate::calc::*;
        pub use crate::data::retrieve_ligands::*;
        pub use crate::mhc::hla::*;
        pub use crate::netmhcpan::netmhcpan_record::*;
        pub use crate::netmhcpan::parser::*;
        pub use crate::netmhcpan::proteome::*;
        pub use crate::reader::*;
    }

    pub mod fs_trait {
        pub use crate::netmhcpan::netmhcpan_record::FindSimilarHLA;
    }

    pub mod external {
        pub use crate::reader::Opt;
        pub use scraper::{Html, Selector};
        pub use structopt::StructOpt;
    }

    pub mod traits {
        pub use crate::mhc::hla::ToDisplay;
        pub use std::convert::TryFrom;
        pub use std::iter::FromIterator;
        pub use std::str::FromStr;
    }

    pub mod nom_tools {
        pub use nom::branch::alt;
        pub use nom::bytes::complete::{tag, take_till, take_until, take_while};
        pub use nom::character::complete::space0;
        pub use nom::character::is_space;
        pub use nom::combinator::{map_res, opt};
        pub use nom::sequence::tuple;
        pub use nom::{AsChar, IResult};
    }

    pub mod collections {
        pub use std::collections::{HashMap, HashSet};
    }

    pub mod logging {
        pub use log::{debug, error, info, trace, warn};
    }
}
