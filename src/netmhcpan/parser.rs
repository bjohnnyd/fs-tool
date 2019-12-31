use crate::prelude::fs_tool::{BindLevel, Deletion, Insertion, NetMHCpanRecord, Peptide};
use crate::prelude::nom_tools::*;
use nom::sequence::tuple;

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

fn get_rank_threshold(input: &str) -> IResult<&str, (&str, f32)> {
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

fn get_element(input: &str) -> IResult<&str, &str> {
    let (i, (_, element)) =
        nom::sequence::tuple((space0, take_while(|c: char| !c.is_whitespace())))(input)?;
    Ok((i, element))
}

/// Need to deal with error
fn get_usize(input: &str) -> IResult<&str, usize> {
    let (i, digit) = nom::combinator::map(get_element, |gp| gp.parse::<usize>())(input)?;

    Ok((i, digit.unwrap()))
}

fn get_bind_level(input: &str) -> IResult<&str, Option<&str>> {
    let get_non_alpha = take_till(&AsChar::is_alpha);
    let get_bind_level = opt(alt((tag("WB"), tag("SB"))));
    let (i, (_, mut bind_level)) = tuple((get_non_alpha, get_bind_level))(input)?;
    if bind_level.is_none() {
        bind_level = Some("NB")
    }

    Ok((i, bind_level))
}

fn process_measures<'a>(input: &'a str) -> IResult<&str, (f32, Option<f32>, f32, &'a str)> {
    let (i, (score, measure1, measure2, mut bind_level)) = tuple((
        opt(get_number),
        opt(get_number),
        opt(get_number),
        get_bind_level,
    ))(input)?;

    if measure2.is_none() {
        Ok(((
            i,
            (
                score.unwrap(),
                measure2,
                measure1.unwrap(),
                bind_level.unwrap(),
            ),
        )))
    } else {
        Ok(((
            i,
            (
                score.unwrap(),
                measure1,
                measure2.unwrap(),
                bind_level.unwrap(),
            ),
        )))
    }
}

/* Need to add creation of Netmhcpan record */
pub fn process_netmhcpan_record<'a>(
    input: &'a str,
) -> IResult<&'a str, (usize, &'a str, usize, usize)> {
    let (i, (pos, hla, pep_seq, core_seq)) =
        tuple((get_usize, get_element, get_element, get_element))(input)?;
    let (i, (offset, del_gp, del_gl, ins_gp, ins_gl)) =
        tuple((get_usize, get_usize, get_usize, get_usize, get_usize))(i)?;
    let (i, (icore, identity)) = tuple((get_element, get_element))(i)?;
    let (i, (score, aff, rank, bind_level)) = process_measures(i)?;

    let pep_info = (pos, pep_seq, core_seq, icore, identity);
    let variation_info = (offset, del_gp, del_gl, ins_gp, ins_gl);
    let binding_info = (hla, score, aff, rank, bind_level);

    Ok((i, (pos, hla, del_gp, del_gl)))
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
                let (_, (threshold_type, threshold)) = get_rank_threshold(line).unwrap();
                dbg!(threshold_type);
                dbg!(threshold);
            }

            if let Some(first_char) = line.chars().take(1).next() {
                if first_char.is_dec_digit() {
                    //                    let netmhcpan_line = line.split_ascii_whitespace().collect::<Vec<&str>>();
                    //                    dbg!(netmhcpan_line.len());
                    println!("STARTS WITH DIGIT\n{}", line);
                    let (i, result) = process_netmhcpan_record(line).unwrap();
                    let (i, _) = tuple((get_element, get_element))(i).unwrap();
                }
            }
        })
    }
}
