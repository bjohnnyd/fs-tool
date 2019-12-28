use nom::sequence::tuple;
use nom::{
    branch::alt,
    bytes::complete::{is_not, tag, take_until, take_while},
    character::{
        complete::{char, digit1, space0},
        is_digit, is_space,
    },
    combinator::{map_res, opt},
    IResult,
};

fn hla_start(input: &str) -> IResult<&str, &str> {
    tag("HLA-")(input)
}

fn take_to_hla_start(input: &str) -> IResult<&str, &str> {
    take_until("HLA-")(input)
}
fn take_until_space(input: &str) -> IResult<&str, &str> {
    is_not(" ")(input)
}

fn take_until_closing_bracket(input: &str) -> IResult<&str, &str> {
    is_not(")")(input)
}
fn take_until_digit(input: &str) -> IResult<&str, &str> {
    take_while(|c| !char::is_digit(c, 10))(input)
}

fn take_float(input: &str) -> IResult<&str, f32> {
    map_res(take_while(|c| char::is_digit(c, 10) || c == '.'), as_float)(input)
}

fn as_float(input: &str) -> Result<f32, std::num::ParseFloatError> {
    input.parse::<f32>()
}

fn is_commented(input: &str) -> IResult<&str, char> {
    char('#')(input)
}

fn space(input: &str) -> IResult<&str, &str> {
    space0(input)
}

fn rank(input: &str) -> IResult<&str, &str> {
    tag("Rank")(input)
}

fn get_nearest_neighbour<T: AsRef<str>>(line: T) -> (String, f32, String) {
    if let Ok((_, (index, _, distance, _, nn))) = tuple((
        take_until_space,
        take_until_digit,
        take_float,
        take_to_hla_start,
        take_until_closing_bracket,
    ))(line.as_ref())
    {
        (index.to_string(), distance, nn.to_string())
    } else {
        ("".to_string(), 0.00, "".to_string())
    }
}
//fn (input: &str) -> IResult<&str, char> {
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hla() {
        let hla =
            "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\

            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000";

        let (index, distance, nn) = get_nearest_neighbour(hla);
        dbg!((index, distance, nn));
        //        let (i, _) = hla_start(hla).unwrap();
        //        let (i, hla) = take_until_space(i).unwrap();
        //        let (i, pre_digit) = take_until_digit(i).unwrap();
        //        let (i, nn) = take_float(i).unwrap();
        //        let (i, spaces) = space(i).unwrap();
        //        dbg!(i);
        //        dbg!(nn);
    }
}
