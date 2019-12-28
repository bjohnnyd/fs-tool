use nom::{
    branch::alt,
    bytes::complete::{is_not, tag, take_till, take_until, take_while},
    character::complete::{char, digit1, space0},
    combinator::{map_res, opt},
    AsChar, IResult,
};

mod helpers {
    pub fn to_float(input: &str) -> Result<f32, std::num::ParseFloatError> {
        input.parse::<f32>()
    }
}

fn get_hla(input: &str) -> IResult<&str, &str> {
    let find_hla = take_until("HLA-");
    let extract_hla = take_while(|c: char| c.is_alphanum() || c == '-' || c == '*' || c == ':');

    let (i, (_, hla)) = nom::sequence::tuple((find_hla, extract_hla))(input)?;
    Ok((i, hla))
}

fn get_number(input: &str) -> IResult<&str, f32> {
    let (i, _) = nom::bytes::complete::take_till(&nom::AsChar::is_dec_digit)(input)?;
    let take_digits = take_while(|c: char| c.is_dec_digit() || c == '.');
    map_res(take_digits, self::helpers::to_float)(i)
}

fn get_rank(input: &str) -> IResult<&str, (&str, f32)> {
    let find_threshold_type = alt((take_until("Strong"), take_until("Weak")));
    let threshold_type = alt((tag("Strong"), tag("Weak")));
    let (i, (_, threshold_type, threshold)) =
        nom::sequence::tuple((find_threshold_type, threshold_type, get_number))(input)?;
    Ok((i, (threshold_type, threshold)))
}

fn get_nn(input: &str) -> IResult<&str, (&str, f32, &str)> {
    let (i, (index, distance, nn)) = nom::sequence::tuple((get_hla, get_number, get_hla))(input)?;
    Ok((i, (index, distance, nn)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use nom::InputTake;

    #[test]
    fn parse_hla() {
        let hla =
            "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\

            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000\n\
            -----------------------------------------------------------------------------------\n\
              Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel\n\
            -----------------------------------------------------------------------------------\n\
                1  HLA-A*03:01     TPQDLNTMLNT  TPQDTMLNT  0  4  2  0  0  TPQDLNTMLNT     Gag_180_209 0.0000160 83.3333\n\
                2  HLA-A*03:01     PQDLNTMLNTV  PQDLNTMLV  0  8  2  0  0  PQDLNTMLNTV     Gag_180_209 0.0000120 87.0000\n\
                3  HLA-A*03:01     QDLNTMLNTVG  QDLNTNTVG  0  5  2  0  0  QDLNTMLNTVG     Gag_180_209 0.0000040 96.0000\
            ";

        hla.lines().map(|line| line.trim()).for_each(|line| {
            if line.starts_with("HLA-") {
                let (_, (index, distance, nn)) = get_nn(line).unwrap();
                dbg!(index);
                dbg!(distance);
                dbg!(nn);
            }
            if line.starts_with("# Rank") {
                let (_, (threshold_type, threshold)) = get_rank(line).unwrap();
                dbg!(threshold_type);
                dbg!(threshold);
            }

            if let Some(first_char) = line.chars().take(1).next() {
                if first_char.is_dec_digit() {
                    let netmhcpan_line = line.split_ascii_whitespace().collect::<Vec<&str>>();
                    dbg!(netmhcpan_line.len());
                    println!("STARTS WITH DIGIT\n{}", line);
                }
            }
        })
    }
}
