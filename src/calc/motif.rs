use crate::error::Error;
use std::{convert::TryFrom, fmt::Display, ops::Range, str::Utf8Error};

/// Represents the motif positions to be used for calculating fraction of shared peptides.
/// Might be extended by a field representing whether the calculations should take KIR genotypes into
/// consideration.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Measure {
    pub name: String,
    pub ranges: Vec<Range<usize>>,
}

impl Measure {
    pub fn new(name: String, s: &str) -> Result<Self, Error> {
        let mut indices = Self::parse_indices(s)?;
        indices.sort();

        let ranges = Self::indices_to_ranges(0, 1, &indices, Vec::new())?;

        Ok(Self { name, ranges })
    }

    pub fn get_motif<'a>(&self, peptide: &'a str) -> Motif<'a> {
        let pep_len = peptide.len();
        Motif(
            self.ranges
                .iter()
                .flat_map(|range| {
                    if range.end < pep_len {
                        Some(&peptide.as_bytes()[*range])
                    } else {
                        None
                    }
                })
                .collect(),
        )
    }

    fn parse_indices(s: &str) -> Result<Vec<usize>, Error> {
        Ok(s.split(',')
            .map(|digit| digit.parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?)
    }

    fn indices_to_ranges(
        start_idx: usize,
        end_idx: usize,
        values: &[usize],
        mut result: Vec<Range<usize>>,
    ) -> Result<Vec<Range<usize>>, Error> {
        if end_idx == values.len() {
            result.push(Range {
                start: values[start_idx],
                end: values[end_idx - 1],
            });
            Ok(result)
        } else if end_idx == 0 || end_idx == start_idx {
            Self::indices_to_ranges(start_idx, end_idx + 1, values, result)
        } else {
            let diff = values[end_idx] - values[end_idx - 1];
            if diff == 1 {
                Self::indices_to_ranges(start_idx, end_idx + 1, values, result)
            } else if diff > 1 {
                result.push(Range {
                    start: values[start_idx],
                    end: values[end_idx - 1],
                });
                Self::indices_to_ranges(end_idx, end_idx + 1, values, result)
            } else {
                return Err(Error::UnsortedIndices);
            }
        }
    }
}

impl std::str::FromStr for Measure {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let name_indices: Vec<&str> = s.split(':').collect();

        if name_indices.len() != 2 {
            Err(Error::IncorrectMeasure(s.to_string()))
        } else {
            let name = &name_indices[0];
            let indices = &name_indices[1];

            Measure::new(name.to_string(), indices)
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Motif<'a>(Vec<&'a [u8]>);

impl TryFrom<Motif<'_>> for String {
    type Error = std::string::FromUtf8Error;

    fn try_from(value: Motif<'_>) -> Result<Self, Self::Error> {
        String::from_utf8(
            value
                .0
                .iter()
                .map(|seq| seq.to_vec())
                .flatten()
                .collect::<Vec<u8>>(),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_motif() {
        let input_measure = "CD8:2,3,4,5,6,9";
        let measure = input_measure.parse::<Measure>().unwrap();

        let input_peptide = "ABCDEFGHIJ";
        let expected = String::from("BCDEFI");

        let motif = measure.get_motif(&input_peptide);

        assert_eq!(String::try_from(motif).unwrap(), expected);
    }

    #[test]
    fn test_get_motif_peptide_too_short() {
        let input_measure = "CD8:2,3,4,5,6,9";
        let measure = input_measure.parse::<Measure>().unwrap();

        let input_peptide = "ABCDEFGH";
        let expected = String::from("BCDEF");

        let motif = measure.get_motif(&input_peptide);

        assert_eq!(String::try_from(motif).unwrap(), expected);
    }
}
