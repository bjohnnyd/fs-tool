// TODO: Need a function to create all combinations
use crate::error::Error;

use immunoprot::ig_like::kir_ligand::{KirLigandMap, LigandMotif};
use immunoprot::mhc::hla::ClassI;
use netmhcpan::result::{BindingData, BindingInfo};

use serde::{Deserialize, Serialize};

/// Represents the motif positions to be used for calculating fraction of shared peptides.
/// Might be extended by a field representing whether the calculations should take KIR genotypes into
/// consideration.
#[derive(Debug, PartialEq, Eq, Hash)]
pub struct Measure {
    pub name: String,
    pub motif_pos: Vec<usize>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct CalcFsResult {
    pub measure: String,
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub index: ClassI,
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub non_index: ClassI,
    #[serde(
        serialize_with = "crate::io::ser::optional_motif_serialize",
        deserialize_with = "crate::io::ser::optional_motif_deserialize"
    )]
    pub index_ligand_motif: Option<LigandMotif>,
    #[serde(
        serialize_with = "crate::io::ser::optional_motif_serialize",
        deserialize_with = "crate::io::ser::optional_motif_deserialize"
    )]
    pub non_index_ligand_motif: Option<LigandMotif>,
    pub fraction_shared: f32,
    pub peptide_length: usize,
    pub index_bound: usize,
    pub non_index_bound: usize,
}

#[derive(Debug)]
pub struct CalculatorComb<'a> {
    pub alleles: (&'a ClassI, &'a ClassI),
    pub binding_data: (Vec<&'a BindingInfo>, Vec<&'a BindingInfo>),
}

impl<'a> CalculatorComb<'a> {
    pub fn new(
        index_allele: &'a ClassI,
        non_index_allele: &'a ClassI,
        index_bd: Vec<&'a BindingInfo>,
        non_index_bd: Vec<&'a BindingInfo>,
    ) -> Self {
        Self {
            alleles: (index_allele, non_index_allele),
            binding_data: (index_bd, non_index_bd),
        }
    }

    pub fn get_motifs(
        &self,
        threshold: f32,
        length: usize,
        aa_pos: &[usize],
    ) -> (Vec<String>, Vec<String>) {
        let bound_motifs = |item: &&BindingInfo| {
            if item.rank() < threshold && item.len() == length {
                Some(item.motif(aa_pos))
            } else {
                None
            }
        };

        let index_motifs = self
            .binding_data
            .0
            .iter()
            .filter_map(bound_motifs)
            .collect::<Vec<String>>();
        let non_index_motifs = self
            .binding_data
            .0
            .iter()
            .filter_map(bound_motifs)
            .collect::<Vec<String>>();

        (index_motifs, non_index_motifs)
    }

    /// Counts number of unique bound motifs (if positions provided) or peptides
    pub fn count_bound(
        &self,
        threshold: f32,
        unique: bool,
        length: usize,
        aa_pos: Option<&[usize]>,
    ) -> (usize, usize) {
        let is_bound = |item: &&&BindingInfo| item.rank() < threshold && item.len() == length;

        let (mut index_motifs, mut non_index_motifs) = match aa_pos {
            Some(aa_pos) => self.get_motifs(threshold, length, aa_pos),
            _ => {
                let index_bound = self
                    .binding_data
                    .0
                    .iter()
                    .filter(is_bound)
                    .map(|info| info.seq().to_string())
                    .collect();
                let non_index_bound = self
                    .binding_data
                    .0
                    .iter()
                    .filter(is_bound)
                    .map(|info| info.seq().to_string())
                    .collect();

                (index_bound, non_index_bound)
            }
        };

        if unique {
            index_motifs.sort();
            non_index_motifs.sort();

            index_motifs.dedup();
            non_index_motifs.dedup();
        }

        (index_motifs.len(), non_index_motifs.len())
    }

    /// Calculates shared peptides motifs based on a threshold for determining bound and whether unique motifs should only be considered,
    /// the peptides considered can also optionally be filtered based on length
    pub fn calculate_shared_motifs(
        &self,
        motif: &[usize],
        threshold: f32,
        unique: bool,
        length: usize,
    ) -> (f32, f32) {
        let (mut index_motifs, mut non_index_motifs) = self.get_motifs(threshold, length, motif);

        if unique {
            index_motifs.sort();
            non_index_motifs.sort();

            index_motifs.dedup();
            non_index_motifs.dedup();
        }

        let index_shared = index_motifs
            .iter()
            .filter(|motif| non_index_motifs.contains(motif))
            .count() as f32;

        let non_index_shared = non_index_motifs
            .iter()
            .filter(|motif| index_motifs.contains(motif))
            .count() as f32;

        let total_bound = (index_motifs.len() + non_index_motifs.len()) as f32;

        (
            index_shared / index_motifs.len() as f32,
            non_index_shared / non_index_motifs.len() as f32,
        )
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

/// Creates all possible allele combinations from NetMHCpan predictions
pub fn create_calc_combs(binding_data: &BindingData) -> Vec<CalculatorComb> {
    binding_data.list_alleles().iter().fold(
        Vec::<CalculatorComb>::new(),
        |mut comb, index_allele| {
            for non_index_allele in binding_data.list_alleles() {
                if *index_allele != non_index_allele {
                    let calc_comb = CalculatorComb {
                        alleles: (index_allele, non_index_allele),
                        binding_data: (
                            binding_data.get_bound_info(index_allele),
                            binding_data.get_bound_info(non_index_allele),
                        ),
                    };

                    comb.push(calc_comb)
                }
            }

            comb
        },
    )
}

/// Make calculations for specific measures and peptide lengths
pub fn calculate_fs(
    combinations: &[CalculatorComb],
    measures: &[Measure],
    ligand_map: &KirLigandMap,
    pep_lengths: &[usize],
    threshold: f32,
    unique: bool,
) -> Vec<CalcFsResult> {
    combinations
        .iter()
        .fold(Vec::<CalcFsResult>::new(), |mut results, comb| {
            let index = comb.alleles.0.clone();
            let non_index = comb.alleles.1.clone();
            let mut index_motifs = ligand_map.get_allele_info(&index);
            let mut non_index_motifs = ligand_map.get_allele_info(&non_index);

            index_motifs.sort();
            non_index_motifs.sort();

            let index_ligand_motif = match index_motifs.iter().next() {
                Some(info) => Some(info.motif().clone()),
                _ => None,
            };

            let non_index_ligand_motif = match non_index_motifs.iter().next() {
                Some(info) => Some(info.motif().clone()),
                _ => None,
            };

            measures.iter().for_each(|measure_group| {
                pep_lengths.iter().for_each(|pep_length| {
                    let measure = measure_group.name.to_string();
                    let motif = &measure_group.motif_pos;
                    let (index_bound, non_index_bound) =
                        comb.count_bound(threshold, unique, *pep_length, Some(motif));
                    let (fraction_shared, _) =
                        comb.calculate_shared_motifs(motif, threshold, unique, *pep_length);

                    let result = CalcFsResult {
                        measure,
                        index: index.clone(),
                        non_index: non_index.clone(),
                        index_ligand_motif: index_ligand_motif.clone(),
                        non_index_ligand_motif: non_index_ligand_motif.clone(),
                        fraction_shared,
                        peptide_length: *pep_length,
                        index_bound,
                        non_index_bound,
                    };

                    results.push(result);
                });
            });

            results
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use netmhcpan::reader::read_raw_netmhcpan;

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

    #[test]
    fn test_calculate_fs() {
        let binding_data = read_raw_netmhcpan("tests/netmhcpan/netmhcpan_wBA.txt").unwrap();
        let comb = create_calc_combs(&binding_data);

        dbg!(&comb);
        dbg!(&comb.len());
    }
}
