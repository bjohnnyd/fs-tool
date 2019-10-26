use nom::{
    IResult,
    sequence,
    bytes::complete,
    multi::separated_list,
    combinator::map_parser,
    Err::Error,
    error::ErrorKind
};

pub(crate) mod parsers {
    use super::*;
    use nom::Err;


    pub(crate) fn not_hla_prefix(hla_name: &str) -> IResult<&str, &str> {
       match sequence::preceded(complete::tag("HLA"),complete::tag("-"))(hla_name) {
           Ok((nomenclature, hla_prefix)) => Ok((nomenclature,"")),
           Err(Error((hla_name, ErrorKind::Tag))) => Ok((hla_name,"")),
           Err(e) => Err(e),
       }
    }

    pub(crate) fn split_hla_fields(hla_name: &str) -> IResult<&str, Vec<&str>> {
        separated_list(
            complete::is_not("-"),complete::tag("A")
        )(hla_name)
    }

}

mod parser_errors {
    #[derive(Default)]
    pub struct ParseError;
    impl std::fmt::Display for ParseError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "A parsing error occurred.")
        }
    }
    impl std::fmt::Debug for ParseError {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            <ParseError as std::fmt::Display>::fmt(self, f)
        }
    }
    impl std::error::Error for ParseError { }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_remove_hla_prefix() {
        assert_eq!(parsers::not_hla_prefix("HLA-A*01:101"), Ok(("A*01:101","")));
        assert_eq!(parsers::not_hla_prefix("A*01:101"), Ok(("A*01:101","")));
//        assert_eq!(parsers::split_hla_fields("HLA-A*01:101"), Ok(("HLA-A*01:101",vec!["HLA-"])));
    }
}
