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
