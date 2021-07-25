pub mod hla;
pub mod motif;

use std::collections::{HashMap, HashSet};

use crate::cohort::Individual;
use immunoprot::ig_like::kir::Kir;
use immunoprot::ig_like::kir_ligand::{KirLigandMap, LigandMotif};
use immunoprot::mhc::hla::ClassI;
use motif::Measure;
use netmhcpan::result::{BindingData, BindingInfo};

use log::{debug, warn};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/* FS */

#[derive(Debug, Serialize, Deserialize, Clone)]
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
    pub binding_data: (&'a [BindingInfo], &'a [BindingInfo]),
}

impl<'a> CalculatorComb<'a> {
    pub fn new(
        index_allele: &'a ClassI,
        non_index_allele: &'a ClassI,
        index_bd: &'a [BindingInfo],
        non_index_bd: &'a [BindingInfo],
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
        let bound_motifs = |item: &BindingInfo| {
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
            .1
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
        let is_bound = |item: &&BindingInfo| item.rank() < threshold && item.len() == length;

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
                    .1
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

        (
            index_shared / index_motifs.len() as f32,
            non_index_shared / non_index_motifs.len() as f32,
        )
    }
}

/// Creates all possible allele combinations from NetMHCpan predictions
/// TODO: IMPORTANT this is memory intensive
pub fn create_calc_combs(binding_data: &BindingData) -> Vec<CalculatorComb> {
    binding_data.list_alleles().par_iter().fold(||
        Vec::<CalculatorComb>::new(),
        |mut comb, index_allele| {
            for non_index_allele in binding_data.list_alleles() {
                if *index_allele != non_index_allele {
                    debug!("Storing binding data info for {}, {}", &index_allele, &non_index_allele);
                    match (binding_data.get_binding_info(index_allele), binding_data.get_binding_info(non_index_allele)) {
                        (Some(index_data), Some(non_index_data)) => {
                            let calc_comb = CalculatorComb {
                                alleles: (index_allele, non_index_allele),
                                binding_data: (
                                    index_data,
                                    non_index_data,
                                ),
                            };
                            comb.push(calc_comb);
                        },
                        _ => warn!("Since binding data is not present for both, no fs will be made between {} and {}", &index_allele, &non_index_allele)
                    }

                }
            }

            comb
        },
    )
        .reduce(|| Vec::new(), |mut a , b| {a.extend(b); a})
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
        .par_iter()
        .fold(
            || Vec::<CalcFsResult>::new(),
            |mut results, comb| {
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
                        debug!(
                            "Calculating FS for index {}, non index {}, measure {}  and length {}",
                            &index, &non_index, &measure_group.name, &pep_length
                        );
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
            },
        )
        .reduce(
            || Vec::new(),
            |mut a, b| {
                a.extend(b);
                a
            },
        )
}

/* Cohort */

// TODO: Need to unit test
pub fn create_index_fs_map(
    index_alleles: Vec<ClassI>,
    fs_results: Vec<CalcFsResult>,
) -> HashMap<ClassI, Vec<CalcFsResult>> {
    index_alleles
        .into_par_iter()
        .fold(|| HashMap::new(), |mut index_calc_map, index_allele| {
            let index_results = fs_results.iter().filter(|result| result.index == index_allele).cloned().collect::<Vec<CalcFsResult>>();
            if index_results.is_empty() {
                warn!("Index allele '{}' was not present in the NetMHCpan results and no calculations for this allele will be performed.", &index_allele);
            } else {
                index_calc_map.insert(index_allele, index_results);
            }
            index_calc_map
        })
        .reduce(|| HashMap::new(), | mut a , b| {a.extend(b); a})
}

#[derive(Debug)]
pub struct IndexCache {
    pub indexes: HashSet<ClassI>,
    pub index_motifs: HashMap<ClassI, LigandMotif>,
    pub index_act_kirs: HashMap<ClassI, Vec<Kir>>,
    pub index_inh_kirs: HashMap<ClassI, Vec<Kir>>,
    pub fs_cache: HashMap<(String, usize), HashMap<(ClassI, ClassI), CalcFsResult>>,
}

impl IndexCache {
    pub fn new(
        index_alleles: Vec<ClassI>,
        fs_result: Vec<CalcFsResult>,
        kir_interactions: &HashMap<Kir, Vec<LigandMotif>>,
        measures: &[Measure],
        pep_lengths: &[usize],
    ) -> Self {
        let mut indexes = HashSet::new();
        let mut index_motifs = HashMap::new();
        let mut index_act_kirs = HashMap::<ClassI, Vec<Kir>>::new();
        let mut index_inh_kirs = HashMap::<ClassI, Vec<Kir>>::new();
        let mut fs_cache =
            HashMap::<(String, usize), HashMap<(ClassI, ClassI), CalcFsResult>>::new();

        let measure_combs = measures.iter().fold(Vec::new(), |mut cache_keys, measure| {
            pep_lengths
                .iter()
                .for_each(|length| cache_keys.push((measure.name.to_string(), *length)));
            cache_keys
        });

        fs_result.into_iter().for_each(|result| {
            if index_alleles.contains(&result.index) {
                let index_allele = result.index.clone();
                if !index_motifs.contains_key(&index_allele) {
                    indexes.insert(index_allele.clone());
                    if let Some(index_motif) = result.index_ligand_motif.clone() {
                        kir_interactions.iter().for_each(|(kir, kir_motifs)| {
                            if kir_motifs.contains(&index_motif) {
                                if kir.is_activating() {
                                    let act_kirs =
                                        index_act_kirs.entry(index_allele.clone()).or_default();
                                    act_kirs.push(kir.clone())
                                } else if kir.is_inhibitory() {
                                    let inh_kirs =
                                        index_inh_kirs.entry(index_allele.clone()).or_default();
                                    inh_kirs.push(kir.clone())
                                }
                            }
                        });

                        index_motifs.insert(index_allele.clone(), index_motif);
                    }
                }

                let index_cache = fs_cache
                    .entry((result.measure.to_string(), result.peptide_length))
                    .or_default();
                index_cache.insert((index_allele, result.non_index.clone()), result);
            }
        });

        index_alleles
            .iter()
            .for_each(|allele| {
                if !indexes.contains(allele) {
                    warn!("Index allele '{}' has no associated NetMHCpan data and will not have cohort results produced", allele);
                }
            });

        Self {
            indexes,
            index_motifs,
            index_act_kirs,
            index_inh_kirs,
            fs_cache,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CohortResult {
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub index: ClassI,
    pub id: String,
    pub measure: String,
    #[serde(serialize_with = "crate::io::ser::optional_float_serialize")]
    pub fs: Option<f32>,
    #[serde(serialize_with = "crate::io::ser::optional_float_serialize")]
    pub ikir_fs: Option<f32>,
    #[serde(serialize_with = "crate::io::ser::optional_float_serialize")]
    pub akir_fs: Option<f32>,
    // TODO: Need to move LILRB outside of this output as it is independent of measure or peptide length
    // and is inefficiently calculated
    #[serde(serialize_with = "crate::io::ser::optional_float_serialize")]
    pub lilrb1: Option<f32>,
    #[serde(serialize_with = "crate::io::ser::optional_float_serialize")]
    pub lilrb2: Option<f32>,
    pub peptide_length: usize,
    pub alleles_considered: usize,
}

/// Temporary function to get kirs bound by allele
pub fn get_bound_kirs(
    motif_interactions: &HashMap<Kir, Vec<LigandMotif>>,
    motif: &LigandMotif,
) -> Vec<Kir> {
    motif_interactions
        .iter()
        .fold(Vec::new(), |mut result, (kir, kir_motifs)| {
            if kir_motifs.contains(motif) {
                result.push(kir.clone())
            };
            result
        })
}

// TODO: Deal with possible errors and also some are never going to return an error
// TODO: Issue when no FS result still LILRB should be produced
pub fn calculate_index_cohort_fs(
    index_cache: IndexCache,
    cohort: &[Individual],
    kir_motif_interactions: &HashMap<Kir, Vec<LigandMotif>>,
    lilrb_scores: &[LilrbScore],
) -> Vec<CohortResult> {
    use LigandMotif::*;

    cohort.par_iter().fold(|| Vec::new(), |mut results, individual| {
        let genotype = &individual.hla_genotype;
        let kir_bound = individual.kir_bound_motifs(&kir_motif_interactions);
        debug!("Started processing individual {}", &individual.id);


        index_cache
            .fs_cache
            .iter()
            .for_each(|((measure, length), calc_result)| {
                index_cache.indexes.iter().for_each(|index| {
                    debug!("Started processing {} for peptide lengths {} with index allele {} and individual {}", &measure, &length, &index, &individual.id);
                    let index_motif = index_cache.index_motifs.get(index);
                    let index_lilrb_scores = lilrb_scores.iter().get_matching(&index);
                    if index_lilrb_scores.is_empty() {
                        warn!("Index allele '{}', has no associated LILRB binding similarity scores", &index);
                    }

                    debug!("Got lilrb scores for index allele");
                    let index_act_kirs = match index_cache.index_act_kirs.get(index) {
                        Some(act_kirs) => act_kirs.clone(),
                        _ => Vec::new(),
                    };

                    let index_inh_kirs = match index_cache.index_inh_kirs.get(index) {
                        Some(inh_kirs) => inh_kirs.clone(),
                        _ => Vec::new(),
                    };

                    debug!("The index HLA being compared {} binds the following activating KIRs: {:?}", &index, &index_act_kirs);
                    debug!("The index HLA being compared {} binds the following inhibitory KIRs: {:?}", &index, &index_act_kirs);

                    let (fs, ikir_fs, akir_fs, lilrb1, lilrb2) = genotype.iter().fold(
                        (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                        |(mut fs, mut ikir_fs, mut akir_fs, mut lilrb1, mut lilrb2), genotype_allele| {

                            if !index_lilrb_scores.is_empty() && index != genotype_allele {
                                let lilrb_scores = index_lilrb_scores.iter().copied().get_matching(genotype_allele);

                                match lilrb_scores.len() {
                                    0 => { warn!("No LILRB binding scores found for allele '{}' in individual {}", &genotype_allele, individual.id) },
                                    1 => {
                                        lilrb1.push(lilrb_scores[0].lilrb1_score);
                                        lilrb2.push(lilrb_scores[0].lilrb2_score);
                                    },
                                    n => {
                                        lilrb1.push(lilrb_scores.iter().map(|score| score.lilrb1_score ).sum::<f32>() / n as f32);
                                        lilrb2.push(lilrb_scores.iter().map(|score| score.lilrb2_score ).sum::<f32>() / n as f32);
                                    },
                                }
                            }

                            if let Some(fs_result) =
                                calc_result.get(&(index.clone(), genotype_allele.clone()))
                            {
                                let initial = fs_result.fraction_shared;
                                let mut akir = initial;
                                let mut ikir = initial;

                                match (index_motif, &fs_result.non_index_ligand_motif) {
                                    (Some(index_motif), Some(gene_motif)) if !index_motif.any_kirs_bound() && !gene_motif.any_kirs_bound()  => {
                                        akir = 1.0;
                                        ikir = 1.0
                                    }
                                    (Some(index_motif), Some(gene_motif)) if (index_motif.any_kirs_bound() && !gene_motif.any_kirs_bound()) || (!index_motif.any_kirs_bound() && gene_motif.any_kirs_bound())  => {
                                        akir = 0.0;
                                        ikir = 0.0
                                    },
                                    (Some(index_motif), Some(gene_motif)) => {
                                        let genotype_bound_kirs =
                                            get_bound_kirs(&kir_motif_interactions, &gene_motif);

                                        debug!("Comparing index motif {} with gene motif {} in individual {}", &index_motif, &gene_motif, &individual.id);
                                        debug!("Individual {} has genotype with HLA {} that binds the following KIRs {:?}.", &individual.id, genotype_allele, &genotype_bound_kirs);
                                        let act_bound = genotype_bound_kirs
                                            .iter()
                                            .filter(|kir| index_act_kirs.contains(kir))
                                            .collect::<Vec<&Kir>>();
                                        debug!("In individual {} index {} shares the following activating kirs with gene allele {}: {:?}", &individual.id, &index, genotype_allele, &act_bound);
                                        let inh_bound = genotype_bound_kirs
                                            .iter()
                                            .filter(|kir| index_inh_kirs.contains(kir))
                                            .collect::<Vec<&Kir>>();
                                        debug!("In individual {} index {} shares the following inhibitory kirs with gene allele {}: {:?}", &individual.id, &index, genotype_allele, &act_bound);

                                        let act_n = act_bound
                                            .iter()
                                            .filter(|kir| individual.kir_genotype.contains(kir))
                                            .count();

                                        debug!("In individual of the activating KIRs bound by the index and the genotype allele {:?}, individual {} has {} of them present.", &act_bound, &individual.id, &act_n);

                                        let inh_n = inh_bound
                                            .iter()
                                            .filter(|kir| individual.kir_genotype.contains(kir))
                                            .count();
                                        debug!("In individual of the inhibitory KIRs bound by the index and the genotype allele {:?}, individual {} has {} of them present.", &inh_bound, &individual.id, &inh_n);

                                        if act_n == 0
                                            && !(index_act_kirs.is_empty() && act_bound.is_empty())
                                        {
                                            akir = 0.0;
                                        }

                                        if inh_n == 0
                                            && !(index_inh_kirs.is_empty() && inh_bound.is_empty())
                                        {
                                            ikir = 0.0;
                                        }
                                    }
                                    // Should this be only (None, None) and throw or ignore otherwise
                                    (None, _) => {
                                        warn!("Index {} has no known ligand motif, setting kir bound calculations to 0 for both activating and inhibitory KIRs", &index);
                                        akir = 0.0;
                                        ikir = 0.0;
                                    }
                                    (_, None) => {
                                        warn!("Allele {} in individual {} has no ligand motif information, setting kir bound calculations to 0 for both activating and inhibitory KIRs", &genotype_allele, &individual.id);
                                        akir = 0.0;
                                        ikir = 0.0;
                                    }
                                }

                                fs.push(initial);
                                ikir_fs.push(ikir);
                                akir_fs.push(akir);
                            }

                            (fs, ikir_fs, akir_fs, lilrb1, lilrb2)
                        },
                    );

                    let alleles_considered = fs.len();

                    let fs = fs
                        .into_iter()
                        .max_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"));
                    let ikir_fs = ikir_fs
                        .into_iter()
                        .max_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"));
                    let akir_fs = akir_fs
                        .into_iter()
                        .max_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"));
                    let lilrb1 = lilrb1
                        .into_iter()
                        .max_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"));
                    let lilrb2 = lilrb2
                        .into_iter()
                        .max_by(|a, b| a.partial_cmp(b).expect("Tried to compare a NaN"));

                    let result = CohortResult {
                        index: index.clone(),
                        id: individual.id.to_string(),
                        measure: measure.to_string(),
                        fs,
                        ikir_fs,
                        akir_fs,
                        lilrb1,
                        lilrb2,
                        peptide_length: *length,
                        alleles_considered,
                    };
                    results.push(result);
                });
            });

        results
    })
        .reduce(|| Vec::new(), | mut a , b| {a.extend(b); a})
}

/* LILRB */

#[derive(Debug)]
pub struct LilrbScore {
    pub first_allele: ClassI,
    pub second_allele: ClassI,
    pub lilrb1_score: f32,
    pub lilrb2_score: f32,
}

trait ExhaustiveSearch<'a> {
    fn get_matching(&self, allele: &ClassI) -> Vec<&'a LilrbScore>;
}

impl<'a, I> ExhaustiveSearch<'a> for I
where
    I: Iterator<Item = &'a LilrbScore> + Clone,
{
    fn get_matching(&self, allele: &ClassI) -> Vec<&'a LilrbScore> {
        let (allele_group, exact) =
            self.clone()
                .fold((Vec::new(), Vec::new()), |(mut group, mut exact), score| {
                    if score.first_allele.allele_group() == allele.allele_group()
                        || score.second_allele.allele_group() == allele.allele_group()
                    {
                        group.push(score);

                        if score.first_allele == *allele || score.second_allele == *allele {
                            exact.push(score)
                        }
                    }

                    (group, exact)
                });

        if exact.is_empty() {
            allele_group
        } else {
            exact
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::reader::read_lilrb_scores;
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
    fn test_create_calculation_combinations() {
        let binding_data = read_raw_netmhcpan(vec![
            "src/tests/input/binding_predictions/netmhcpan_wBA.txt",
        ])
        .unwrap();
        let comb = create_calc_combs(&binding_data);
        dbg!(&comb);

        let index_allele = comb[0].alleles.0.to_nomenclature_string();
        let nonindex_allele = comb[0].alleles.1.to_nomenclature_string();

        assert_eq!(index_allele, "HLA-A*03:01");
        assert_eq!(nonindex_allele, "HLA-B*27:05");
    }

    #[test]
    fn test_lilrb_search() {
        let lilrb_scores = read_lilrb_scores();
        let test_allele = "A*03:02".parse::<ClassI>().unwrap();

        assert!(lilrb_scores
            .iter()
            .get_matching(&test_allele)
            .iter()
            .all(|lilr_info| {
                lilr_info.first_allele.allele_group() == test_allele.allele_group()
                    || lilr_info.second_allele.allele_group() == test_allele.allele_group()
            }));
    }
}
