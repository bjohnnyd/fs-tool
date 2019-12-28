use nom::{
    branch::alt,
    bytes::complete::{is_not, tag, take_until, take_while},
    character::{
        complete::{char, digit1},
        is_digit, is_space,
    },
    combinator::{map_res, opt},
    IResult,
};

fn hla_start(input: &str) -> IResult<&str, &str> {
    tag("HLA-")(input)
}

fn optional_asterisk(input: &str) -> IResult<&str, Option<char>> {
    opt(char('*'))(input)
}

fn colon(input: &str) -> IResult<&str, char> {
    char(':')(input)
}
fn is_gene_locus(input: &str) -> IResult<&str, char> {
    alt((char('A'), char('B'), char('C')))(input)
}

fn take_digits(input: &str) -> IResult<&str, &str> {
    digit1(input)
}

fn take_until_space(input: &str) -> IResult<&str, &str> {
    is_not(" ")(input)
}

fn take_until_digit(input: &str) -> IResult<&str, &str> {
    take_while(|c| !char::is_digit(c, 10))(input)
}

fn take_nn(input: &str) -> IResult<&str, f32> {
    map_res(take_while(|c| char::is_digit(c, 10) || c == '.'), as_float)(input)
}

fn as_float(input: &str) -> Result<f32, std::num::ParseFloatError> {
    input.parse::<f32>()
}

//fn (input: &str) -> IResult<&str, char> {
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hla() {
        let hla =
            "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)";

        let (i, _) = hla_start(hla).unwrap();
        let (i, hla) = take_until_space(i).unwrap();
        let (i, pre_digit) = take_until_digit(i).unwrap();
        let (i, nn) = take_nn(i).unwrap();
        dbg!(i);
        dbg!(nn);
    }
}
