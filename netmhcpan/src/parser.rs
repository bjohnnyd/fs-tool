pub const NETMHCPAN_VERSION: &str = "4.0";

use nom::{
    bytes::complete::{tag, take_until, take_while, take_while1},
    character::is_alphanumeric,
    combinator::opt,
    sequence::tuple,
    IResult,
};

use crate::result::{NearestNeighbour, RankThreshold};
use immunoprot::mhc::hla::ClassI;

/* Basic Parsers */
fn take_first_numeric(i: &[u8]) -> IResult<&[u8], f32> {
    let take_until_digit = take_while(|c: u8| c.is_ascii_alphabetic() || c.is_ascii_whitespace());
    let take_digits = take_while1(|c: u8| c.is_ascii_digit() || c.is_ascii_punctuation());

    let (remainder, (_, numeric_word)) = tuple((take_until_digit, take_digits))(i)?;

    let numeric = String::from_utf8(numeric_word.to_vec())
        .unwrap()
        .parse::<f32>()
        .unwrap();

    Ok((remainder, numeric))
}

fn take_word(i: &[u8]) -> IResult<&[u8], &[u8]> {
    let word = take_while1(is_alphanumeric);
    let space = take_while1(|c| c == b' ');

    let (remainder, (_, word)) = tuple((space, word))(i)?;

    Ok((remainder, word))
}

// TODO: Need to deal with error
fn take_hla_allele(i: &[u8]) -> IResult<&[u8], ClassI> {
    let until_allele = take_until("HLA-");
    let (remainder, (_, hla)) = tuple((until_allele, take_word))(i)?;
    dbg!(&hla);

    let hla_allele = String::from_utf8(hla.to_vec())
        .unwrap()
        .parse::<ClassI>()
        .unwrap();

    Ok((remainder, hla_allele))
}

/* Line Identifiers */

fn is_nn_line(i: &[u8]) -> IResult<&[u8], Option<&[u8]>> {
    opt(tag(b"HLA-"))(i)
}

fn is_rank_line(i: &[u8]) -> IResult<&[u8], Option<&[u8]>> {
    opt(tag(b"# Rank Threshold"))(i)
}

/* Line Parsers */

// TODO: Error fixing needed, especially regarding match
fn get_rank_info(i: &[u8]) -> IResult<&[u8], RankThreshold> {
    use RankThreshold::*;

    let take_until_rank_threshold =
        take_while(|c: u8| c.is_ascii_alphabetic() || c.is_ascii_whitespace());
    let rank_threshold = take_while1(|c: u8| c.is_ascii_digit() || c.is_ascii_punctuation());
    let (remainder, (_, rank_type, _, rank_threshold)) = tuple((
        take_word,
        take_word,
        take_until_rank_threshold,
        rank_threshold,
    ))(i)
    .unwrap();

    let rank_threshold = String::from_utf8(rank_threshold.to_vec())
        .unwrap()
        .parse::<f32>()
        .unwrap();
    let rank = match rank_type {
        b"Strong" => Strong(rank_threshold),
        b"Weak" => Weak(rank_threshold),
        _ => Weak(0f32),
    };

    Ok((remainder, rank))
}

fn get_nn_info(i: &[u8]) -> IResult<&[u8], NearestNeighbour> {
    let (remainder, (index, distance, nn)) =
        tuple((take_hla_allele, take_first_numeric, take_hla_allele))(i).unwrap();
    let nn_info = NearestNeighbour::new(index, distance, nn);

    Ok((remainder, nn_info))
}

#[cfg(test)]
mod tests {
    use crate::parser::{
        get_nn_info, get_rank_info, is_nn_line, is_rank_line, take_first_numeric, take_word,
    };
    use crate::result::{NearestNeighbour, RankThreshold};

    #[test]
    fn test_parse_numeric() {
        assert_eq!(take_first_numeric(b"word5"), Ok((&b""[..], 5f32)));
        assert_eq!(take_first_numeric(b"5word"), Ok((&b"word"[..], 5f32)));
    }

    #[test]
    fn test_identify_nn_line() {
        let netmhcout = std::fs::read_to_string("tests/netmhcpan_woBA.txt").unwrap();

        let nn_neighbours =
            netmhcout
                .lines()
                .fold(Vec::<NearestNeighbour>::new(), |mut nn_neighbours, line| {
                    match is_nn_line(line.as_bytes()).unwrap() {
                        (nn_line, Some(_)) => {
                            let (_, mut nn) = get_nn_info(nn_line).unwrap();
                            nn_neighbours.push(nn)
                        }
                        _ => (),
                    }
                    nn_neighbours
                });

        dbg!(nn_neighbours);
    }

    #[test]
    fn test_rank() {
        use RankThreshold::*;
        let netmhcout = std::fs::read_to_string("tests/netmhcpan_woBA.txt").unwrap();
        let expected = vec![
            Strong(0.5_f32),
            Weak(2.0_f32),
            Strong(0.5_f32),
            Weak(2.0_f32),
        ];

        let mut thresholds =
            netmhcout
                .lines()
                .fold(Vec::<RankThreshold>::new(), |mut thresholds, line| {
                    match is_rank_line(line.as_bytes()).unwrap() {
                        (rank_line, Some(_)) => {
                            let (_, mut rank) = get_rank_info(rank_line).unwrap();
                            thresholds.push(rank)
                        }
                        _ => (),
                    }
                    thresholds
                });
        assert_eq!(thresholds, expected);
    }
}
