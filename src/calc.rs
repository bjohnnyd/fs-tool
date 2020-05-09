// TODO: Need a function to create all combinations
use crate::error::Error;

use immunoprot::mhc::hla::ClassI;
use netmhcpan::result::{BindingInfo, Peptide};

/// Represents the motif positions to be used for calculating fraction of shared peptides.
/// Might be extended by a field representing whether the calculations should take KIR genotypes into
/// consideration.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Measure {
    pub name: String,
    pub motif_pos: Vec<usize>,
}

#[derive(Debug)]
pub struct CalculatorComb<'a> {
    pub alleles: (&'a ClassI, &'a ClassI),
    pub binding_data: (&'a [BindingInfo], &'a [BindingInfo]),
}

impl<'a> CalculatorComb<'a> {
    fn new(
        first_allele: &'a ClassI,
        second_allele: &'a ClassI,
        first_bd: &'a Vec<BindingInfo>,
        second_bd: &'a Vec<BindingInfo>,
    ) -> Self {
        Self {
            alleles: (first_allele, second_allele),
            binding_data: (first_bd, second_bd),
        }
    }

    /// Calculates shared peptides motifs based on a threshold for determining bound and whether unique motifs should only be considered,
    /// the peptides considered can also optionally be filtered based on length
    fn calculate_shared_motifs(
        &self,
        motif: &[usize],
        threshold: f32,
        unique: bool,
        length: Option<usize>,
    ) -> (f32, f32) {
        let get_bound_motifs = |item: &BindingInfo| {
            if item.rank() < threshold {
                match length {
                    Some(length) => {
                        let pep = item.peptide();
                        if pep.len() == length {
                            Some(pep.sequence_motif(motif))
                        } else {
                            None
                        }
                    }
                    None => Some(item.peptide().sequence_motif(motif)),
                }
            } else {
                None
            }
        };

        let mut first_motifs = self
            .binding_data
            .0
            .iter()
            .filter_map(get_bound_motifs)
            .collect::<Vec<String>>();
        first_motifs.sort();

        let mut second_motifs = self
            .binding_data
            .0
            .iter()
            .filter_map(get_bound_motifs)
            .collect::<Vec<String>>();
        second_motifs.sort();

        if unique {
            first_motifs.dedup();
            second_motifs.dedup();
        }

        let first_shared = first_motifs
            .iter()
            .filter(|motif| second_motifs.contains(motif))
            .count() as f32;
        let second_shared = second_motifs
            .iter()
            .filter(|motif| first_motifs.contains(motif))
            .count() as f32;

        let total_bound = (first_motifs.len() + second_motifs.len()) as f32;

        (first_shared / total_bound, second_shared / total_bound)
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    const NETMHCPAN_RESULT: &str = r#"HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\
            \n\
            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000\n\
            -----------------------------------------------------------------------------------\n\
              Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel\n\
            -----------------------------------------------------------------------------------\n\
                1  HLA-A*03:01     TPQDLNTMLNT  TPQDTMLNT  0  4  2  0  0  TPQDLNTMLNT     Gag_180_209 0.0000160 83.3333\n\
                2  HLA-A*03:01     PQDLNTMLNTV  PQDLNTMLV  0  8  2  0  0  PQDLNTMLNTV     Gag_180_209 0.0000120 87.0000\n\
                3  HLA-A*03:01     QDLNTMLNTVG  QDLNTNTVG  0  5  2  0  0  QDLNTMLNTVG     Gag_180_209 0.0000040 96.0000\n\
                4  HLA-A*03:01     DLNTMLNTVGG  DLNTNTVGG  0  4  2  0  0  DLNTMLNTVGG     Gag_180_209 0.0000040 96.0000\n\
                5  HLA-A*03:01     LNTMLNTVGGH  LMLNTVGGH  0  1  2  0  0  LNTMLNTVGGH     Gag_180_209 0.0001090 48.8333\n\
                6  HLA-A*03:01     NTMLNTVGGHQ  NTMTVGGHQ  0  3  2  0  0  NTMLNTVGGHQ     Gag_180_209 0.0001260 46.1429\n\
                7  HLA-A*03:01     TMLNTVGGHQA  TMLNGGHQA  0  4  2  0  0  TMLNTVGGHQA     Gag_180_209 0.0001260 46.1429\n\
                8  HLA-A*03:01     MLNTVGGHQAA  MLNGGHQAA  0  3  2  0  0  MLNTVGGHQAA     Gag_180_209 0.0002300 36.4000\n\
                9  HLA-A*03:01     LNTVGGHQAAM  NTVGGHQAM  1  7  1  0  0   NTVGGHQAAM     Gag_180_209 0.0000530 62.5000\n\
               10  HLA-A*03:01     NTVGGHQAAMQ  NTVGGAAMQ  0  5  2  0  0  NTVGGHQAAMQ     Gag_180_209 0.0001420 44.1250\n\
               11  HLA-A*03:01     TVGGHQAAMQM  TVHQAAMQM  0  2  2  0  0  TVGGHQAAMQM     Gag_180_209 0.0004120 28.5000\n\
               12  HLA-A*03:01     VGGHQAAMQML  VGGAAMQML  0  3  2  0  0  VGGHQAAMQML     Gag_180_209 0.0000120 87.0000\n\
               13  HLA-A*03:01     GGHQAAMQMLK  GQAAMQMLK  0  1  2  0  0  GGHQAAMQMLK     Gag_180_209 0.0313010  3.8215\n\
               14  HLA-A*03:01     GHQAAMQMLKE  GQAAMQMLK  0  1  1  0  0   GHQAAMQMLK     Gag_180_209 0.0004080 28.6176\n\
               15  HLA-A*03:01     HQAAMQMLKET  HQAAMQMLK  0  0  0  0  0    HQAAMQMLK     Gag_180_209 0.0003110 32.0000\n\
               16  HLA-A*03:01     QAAMQMLKETI  QAAMQMLTI  0  7  2  0  0  QAAMQMLKETI     Gag_180_209 0.0000140 85.0000\n\
               17  HLA-A*03:01     AAMQMLKETIN  AAMQKETIN  0  4  2  0  0  AAMQMLKETIN     Gag_180_209 0.0000060 93.7500\n\
               18  HLA-A*03:01     AMQMLKETINE  AMLKETINE  0  1  2  0  0  AMQMLKETINE     Gag_180_209 0.0001620 41.9000\n\
               19  HLA-A*03:01     MQMLKETINEE  MLKETINEE  0  1  2  0  0  MQMLKETINEE     Gag_180_209 0.0000850 53.5417\n\
               20  HLA-A*03:01     QMLKETINEEA  QMLKETINA  0  8  2  0  0  QMLKETINEEA     Gag_180_209 0.0000410 67.2727\n\
               21  HLA-A*03:01       MLKETINEE  MLKETINEE  0  0  0  0  0    MLKETINEE     Gag_180_209 0.0079270  7.4157\n\
               22  HLA-A*03:01       LKETINEEA  LKETINEEA  0  0  0  0  0    LKETINEEA     Gag_180_209 0.0000450 65.4545\n\
                "#;

    #[test]
    fn test_create_measure() {
        let input_measure = "CD8:2,3,4,5,6,9";
        let measure = input_measure.parse::<Measure>().unwrap();

        assert_eq!(
            measure,
            Measure {
                name: "CD8".to_string(),
                motif_pos: vec![2, 3, 4, 5, 6, 9]
            }
        )
    }
}
