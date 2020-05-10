use immunoprot::ig_like::kir_ligand::{AlleleFreq, KirLigandInfo, KirLigandMap, LigandMotif};
use immunoprot::mhc::hla::ClassI;
use netmhcpan::result::BindingData;

use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;

#[derive(Debug)]
pub struct LigandMeta {
    pub kir_ligand_allele: ClassI,
    pub kir_ligand_motif: LigandMotif,
    pub kir_ligand_allele_freq: AlleleFreq,
}

impl LigandMeta {
    pub(crate) fn new(ligand_info: &KirLigandInfo) -> Self {
        Self {
            kir_ligand_allele: ligand_info.allele().clone(),
            kir_ligand_motif: ligand_info.motif().clone(),
            kir_ligand_allele_freq: ligand_info.freq().clone(),
        }
    }
}

#[skip_serializing_none]
#[derive(Debug, Deserialize)]
pub struct AlleleMeta {
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub allele: ClassI,
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub netmhcpan_nn: ClassI,
    pub netmhcpan_nn_distance: f32,
    pub ligand_meta: Option<LigandMeta>,
}

#[skip_serializing_none]
#[derive(Debug, Serialize, Deserialize)]
pub struct BindingMeta {
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub allele: ClassI,
    #[serde(with = "serde_with::rust::display_fromstr")]
    pub protein: String,
    pub n_strong_bound: usize,
    pub n_weak_bound: usize,
    pub pep_length: usize,
}

pub fn create_allele_metadata(
    binding_data: &BindingData,
    kir_ligand_map: &KirLigandMap,
) -> Vec<AlleleMeta> {
    binding_data
        .list_nn()
        .iter()
        .map(|nn| {
            let (allele, netmhcpan_nn_distance, netmhcpan_nn) = nn.info();
            let mut ligand_info = kir_ligand_map.get_allele_info(allele);
            ligand_info.sort();

            let ligand_meta = match ligand_info.iter().next() {
                Some(ligand_info) => Some(LigandMeta::new(ligand_info)),
                None => None,
            };

            AlleleMeta {
                allele: allele.clone(),
                netmhcpan_nn: netmhcpan_nn.clone(),
                netmhcpan_nn_distance,
                ligand_meta,
            }
        })
        .collect::<Vec<AlleleMeta>>()
}

pub fn create_binding_metadata(binding_data: &BindingData) -> Vec<BindingMeta> {
    let mut binding_meta = Vec::<BindingMeta>::new();
    let proteins = binding_data.proteins();
    let pep_lengths = binding_data.pep_lengths();
    let alleles = binding_data.list_alleles();

    let strong_threshold = binding_data.strong_threshold();
    let weak_threshold = binding_data.weak_threshold();

    alleles.iter().for_each(|allele| {
        let binding_info = binding_data.get_bound_info(allele);

        proteins.iter().for_each(|protein| {
            pep_lengths.iter().for_each(|pep_length| {
                let mut n_strong_bound = 0;
                let mut n_weak_bound = 0;

                binding_info.iter().for_each(|binding_info| {
                    let peptide = binding_info.peptide();

                    if peptide.len() == *pep_length && peptide.protein() == protein {
                        let binding_rank = binding_info.rank();

                        if binding_rank < weak_threshold {
                            n_weak_bound += 1;
                        }

                        if binding_rank < strong_threshold {
                            n_strong_bound += 1;
                        }
                    }
                });

                binding_meta.push(BindingMeta {
                    allele: (*allele).clone(),
                    protein: protein.clone(),
                    n_strong_bound,
                    n_weak_bound,
                    pep_length: *pep_length,
                })
            })
        })
    });

    binding_meta
}
