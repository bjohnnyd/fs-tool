pub const NETMHCPAN_VERSION: &str = "4.0";
pub const HLA_GENES: [u8; 3] = [b'A', b'B', b'C'];
pub const HLA_GENE_SEPARATORS: [u8; 3] = [b':', b'*', b'-'];

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
    let take_until_digit = take_while(|c: u8| !c.is_ascii_digit());
    let take_digits = take_while1(|c: u8| c.is_ascii_digit() || c.is_ascii_punctuation());

    let (remainder, (_, numeric_word)) = tuple((take_until_digit, take_digits))(i)?;

    let numeric = String::from_utf8(numeric_word.to_vec())
        .unwrap()
        .parse::<f32>()
        .unwrap();

    Ok((remainder, numeric))
}

fn take_word(i: &[u8]) -> IResult<&[u8], &[u8]> {
    let word = take_while(is_alphanumeric);
    let space = take_while(|c| c == b' ');

    let (remainder, (_, word)) = tuple((space, word))(i)?;

    Ok((remainder, word))
}

// TODO: Need to deal with error
fn take_hla_allele(i: &[u8]) -> IResult<&[u8], ClassI> {
    let allele_prefix = opt(tag("HLA-"));
    let take_allele = take_while(|c: u8| {
        c.is_ascii_digit() || HLA_GENES.contains(&c) || HLA_GENE_SEPARATORS.contains(&c)
    });

    let (remainder, (_, hla)) = tuple((allele_prefix, take_allele))(i)?;
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

fn is_peptide_line(i: &[u8]) -> IResult<&[u8], bool> {
    if let (non_space, Some(_)) = opt(take_while(|c: u8| c.is_ascii_whitespace()))(i)? {
        if !non_space.is_empty() {
            let is_digit = non_space[0].is_ascii_digit();
            return Ok((non_space, is_digit));
        }
    }

    Ok((i, false))
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
    let (remainder, (index, distance, _, nn)) = tuple((
        take_hla_allele,
        take_first_numeric,
        take_until("HLA-"),
        take_hla_allele,
    ))(i)
    .unwrap();
    let nn_info = NearestNeighbour::new(index, distance, nn);

    Ok((remainder, nn_info))
}

#[cfg(test)]
mod tests {
    use crate::parser::*;
    use crate::result::{NearestNeighbour, RankThreshold};

    #[test]
    fn test_parse_numeric() {
        assert_eq!(take_first_numeric(b"word5"), Ok((&b""[..], 5f32)));
        assert_eq!(take_first_numeric(b"5word"), Ok((&b"word"[..], 5f32)));
    }

    #[test]
    fn test_pep_line() {
        let netmhcout = std::fs::read_to_string("tests/netmhcpan_woBA.txt").unwrap();

        netmhcout
            .lines()
            .for_each(|line| match is_peptide_line(line.as_bytes()) {
                Ok((pep_line, true)) => {
                    println!("{}", String::from_utf8(pep_line.to_vec()).unwrap())
                }
                _ => (),
            });
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
                            let (_, nn) = get_nn_info(nn_line).unwrap();
                            nn_neighbours.push(nn)
                        }
                        _ => (),
                    }
                    nn_neighbours
                });
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

        let thresholds =
            netmhcout
                .lines()
                .fold(Vec::<RankThreshold>::new(), |mut thresholds, line| {
                    match is_rank_line(line.as_bytes()).unwrap() {
                        (rank_line, Some(_)) => {
                            let (_, rank) = get_rank_info(rank_line).unwrap();
                            thresholds.push(rank)
                        }
                        _ => (),
                    }
                    thresholds
                });
        assert_eq!(thresholds, expected);
    }
}
