// TODO: Deal with errors from nom to custom
pub const NETMHCPAN_VERSION: &str = "4.0";
pub const HLA_GENES: [u8; 3] = [b'A', b'B', b'C'];
pub const HLA_GENE_SEPARATORS: [u8; 3] = [b':', b'*', b'-'];

use nom::{
    bytes::complete::{tag, take_until, take_while, take_while1},
    combinator::opt,
    multi::many_m_n,
    sequence::tuple,
    IResult,
};

use crate::result::{BindingInfo, NearestNeighbour, Peptide, RankThreshold};

use immunoprot::mhc::hla::ClassI;

/* Basic Parsers */
pub fn take_first_numeric(i: &[u8]) -> IResult<&[u8], String> {
    let take_until_digit = take_while(|c: u8| !c.is_ascii_digit());
    let take_digits = take_while1(|c: u8| c.is_ascii_digit() || c.is_ascii_punctuation());

    let (remainder, (_, numeric_word)) = tuple((take_until_digit, take_digits))(i)?;

    let numeric = String::from_utf8(numeric_word.to_vec()).unwrap();

    Ok((remainder, numeric))
}

pub fn take_word(i: &[u8]) -> IResult<&[u8], String> {
    let word = take_while(|c: u8| !c.is_ascii_whitespace());
    let space = take_while(|c| c == b' ');

    let (remainder, (_, word)) = tuple((space, word))(i)?;
    let word = String::from_utf8(word.to_vec()).unwrap();

    Ok((remainder, word))
}

// TODO: Need to deal with error
pub fn take_hla_allele(i: &[u8]) -> IResult<&[u8], ClassI> {
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

pub fn is_nn_line(i: &[u8]) -> IResult<&[u8], Option<&[u8]>> {
    opt(tag(b"HLA-"))(i)
}

pub fn is_rank_line(i: &[u8]) -> IResult<&[u8], Option<&[u8]>> {
    opt(tag(b"# Rank Threshold"))(i)
}

pub fn is_peptide_line(i: &[u8]) -> IResult<&[u8], bool> {
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
pub fn get_rank_info(i: &[u8]) -> IResult<&[u8], RankThreshold> {
    use RankThreshold::*;

    let take_until_rank_threshold =
        take_while(|c: u8| c.is_ascii_alphabetic() || c.is_ascii_whitespace());
    let (remainder, (_, rank_type, _, rank_threshold)) = tuple((
        take_word,
        take_word,
        take_until_rank_threshold,
        take_first_numeric,
    ))(i)
    .unwrap();

    let rank_threshold = rank_threshold.parse::<f32>().unwrap();

    let rank = match rank_type.as_ref() {
        "Strong" => Strong(rank_threshold),
        "Weak" => Weak(rank_threshold),
        _ => Weak(0f32),
    };

    Ok((remainder, rank))
}

pub fn get_nn_info(i: &[u8]) -> IResult<&[u8], NearestNeighbour> {
    let (remainder, (index, distance, _, nn)) = tuple((
        take_hla_allele,
        take_first_numeric,
        take_until("HLA-"),
        take_hla_allele,
    ))(i)
    .unwrap();

    let distance = distance.parse::<f32>().unwrap();
    let nn_info = NearestNeighbour::new(index, distance, nn);

    Ok((remainder, nn_info))
}

pub fn get_netmhc_entry_info(i: &[u8]) -> IResult<&[u8], (usize, ClassI, String)> {
    let (remainder, (pos, _, allele, pep_seq)) = tuple((
        take_first_numeric,
        take_until("HLA-"),
        take_hla_allele,
        take_word,
    ))(i)?;

    // Shifts position by `-1` due to NetMHCpan entry representation being 1-based
    let pos = pos.parse::<usize>().unwrap() - 1;

    Ok((remainder, (pos, allele, pep_seq)))
}

pub fn get_netmhc_align_info(i: &[u8]) -> IResult<&[u8], (Vec<usize>, String, String)> {
    let (i, alignment_mods) = many_m_n(5, 5, take_first_numeric)(i)?;
    let (remainder, (icore, identity)) = tuple((take_word, take_word))(i)?;

    let alignment_mods = alignment_mods
        .iter()
        .map(|num| num.parse::<usize>().unwrap())
        .collect::<Vec<usize>>();

    Ok((remainder, (alignment_mods, icore, identity)))
}

pub fn get_netmhc_binding_info<'a, 'b>(
    i: &'b [u8],
    peptide: Peptide,
) -> IResult<&'b [u8], BindingInfo> {
    let (remainder, binding_info) = many_m_n(2, 3, take_first_numeric)(i)?;

    let binding_info = binding_info
        .iter()
        .map(|num| num.parse::<f32>().unwrap())
        .collect::<Vec<f32>>();

    let score = binding_info[0];

    let (affinity, rank) = if binding_info.len() == 2 {
        (None, binding_info[1])
    } else {
        (Some(binding_info[1]), binding_info[2])
    };

    Ok((
        remainder,
        BindingInfo {
            peptide,
            score,
            affinity,
            rank,
        },
    ))
}

#[cfg(test)]
mod tests {
    use crate::parser::*;
    use crate::result::{NearestNeighbour, Peptide, Protein, RankThreshold};

    static TEST_ENTRY: &[u8] = b"  1  HLA-A*03:01     TPQDLNTMLNT  TPLNTMLNT  0  2  2  0  0  TPQDLNTMLNT     Gag_180_209 0.0190370 40692.6 77.6355";

    #[test]
    fn test_parse_numeric() {
        assert_eq!(
            take_first_numeric(b"word5"),
            Ok((&b""[..], String::from("5")))
        );
        assert_eq!(
            take_first_numeric(b"5word"),
            Ok((&b"word"[..], String::from("5")))
        );
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

    #[test]
    fn test_parse_entry() {
        let (i, entry_info) = get_netmhc_entry_info(TEST_ENTRY).unwrap();
        let (i, alignment_info) = get_netmhc_align_info(i).unwrap();

        let mut protein = Protein::new(alignment_info.2.clone());
        protein.add_sequence_at_pos(entry_info.0, &entry_info.2);
        let peptide = Peptide::new(
            entry_info.0,
            entry_info.2.len(),
            entry_info.2,
            alignment_info.2,
            alignment_info.1,
            &alignment_info.0,
        );

        let (i, binding_info) = get_netmhc_binding_info(i, peptide).unwrap();
        dbg!(&binding_info);
    }
}
