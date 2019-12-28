use nom::sequence::tuple;
use nom::AsChar;
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

mod helpers {
    pub fn to_float(input: &str) -> Result<f32, std::num::ParseFloatError> {
        input.parse::<f32>()
    }
}

fn get_hla(input: &str) -> IResult<&str, &str> {
    let find_hla = take_until("HLA-");
    let extract_hla = take_while(|c: char| c.is_alphanum() || c == '-' || c == '*' || c == ':');

    let (i, (_, hla)) = tuple((find_hla, extract_hla))(input)?;
    Ok((i, hla))
}

fn get_number(input: &str) -> IResult<&str, f32> {
    let (i, _) = nom::bytes::complete::take_till(&nom::AsChar::is_dec_digit)(input)?;
    map_res(digit1, self::helpers::to_float)(i)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_hla() {
        let hla =
            "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\

            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000";

        let (i, (index, distance, nn)) = tuple((get_hla, get_number, get_hla))(hla).unwrap();
        dbg!(i);
        dbg!(distance);
        dbg!(index);
        dbg!(nn);
    }
}
