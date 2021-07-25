use crate::error::Error;

/// Represents the motif positions to be used for calculating fraction of shared peptides.
/// Might be extended by a field representing whether the calculations should take KIR genotypes into
/// consideration.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Measure {
    pub name: String,
    pub motif_pos: Vec<usize>,
}

impl Measure {
    fn parse_indices(s: &str) -> Result<Vec<usize>, Error> {
        Ok(s.split(',')
            .map(|digit| digit.parse::<usize>())
            .collect::<Result<Vec<_>, _>>()?)
    }

    fn get_motif(&self, peptide: &str) -> Vec<u8> {
        let pep_len = peptide.len();
        self.motif_pos
            .iter()
            .flat_map(|idx| {
                if *idx < pep_len {
                    Some(peptide.as_bytes()[*idx])
                } else {
                    None
                }
            })
            .collect()
    }
}

impl std::str::FromStr for Measure {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut name = String::new();
        let mut motif_pos = Vec::<usize>::new();

        let mut name_pos = s.split(':');

        if let Some(field) = name_pos.next() {
            match name_pos.next() {
                Some(measure_pos) => {
                    name = field.to_string();
                    motif_pos = Measure::parse_indices(measure_pos)?
                }
                _ => {
                    motif_pos = Measure::parse_indices(field)?;
                    name = motif_pos
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<String>>()
                        .join("_");
                }
            }
        }
        Ok(Self { name, motif_pos })
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
        let expected = input_peptide
            .bytes()
            .enumerate()
            .filter(|(i, _)| measure.motif_pos.contains(i))
            .map(|(_, aa)| aa)
            .collect::<Vec<u8>>();

        let motif = measure.get_motif(&input_peptide);

        assert_eq!(motif, expected);
    }

    #[test]
    fn test_get_motif_peptide_too_short() {
        let input_measure = "CD8:2,3,4,5,6,9";
        let measure = input_measure.parse::<Measure>().unwrap();

        let input_peptide = "ABCDEFGH";
        let expected = input_peptide
            .bytes()
            .enumerate()
            .filter(|(i, _)| measure.motif_pos.contains(i))
            .map(|(_, aa)| aa)
            .collect::<Vec<u8>>();

        let motif = measure.get_motif(&input_peptide);

        assert_eq!(motif, expected);
    }
}
