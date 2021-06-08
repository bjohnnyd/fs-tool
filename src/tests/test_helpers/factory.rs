use immunoprot::ig_like::kir_ligand::KirLigandMap;
use netmhcpan::result::{BindingInfo, Peptide};
use rand::{prelude::SliceRandom, Rng};

pub struct AlleleFactory {
    kir_ligand_map: KirLigandMap,
}

enum RankThreshold {
    Strong,
    Weak,
}

pub struct BindingFactory {
    peptides: Vec<Peptide>,
}

impl BindingFactory {
    fn generate_binding_info(
        &self,
        n: usize,
        rank_class: Option<RankThreshold>,
    ) -> Vec<BindingInfo> {
        let mut rng = rand::thread_rng();
        let mut peptides = self.peptides.to_vec();
        let score = rng.gen::<f32>();
        let affinity = None;

        peptides.shuffle(&mut rng);
        peptides = peptides.into_iter().take(n).collect();

        peptides
            .into_iter()
            .map(|peptide| {
                let rank: f32 = match rank_class {
                    Some(RankThreshold::Strong) => rng.gen::<f32>() * 0.5f32,
                    Some(RankThreshold::Weak) => rng.gen::<f32>() * 2f32,
                    _ => rng.gen::<f32>() + 3f32,
                };

                BindingInfo {
                    peptide,
                    score,
                    affinity,
                    rank,
                }
            })
            .collect()
    }
}

pub fn generate_peptides(n: usize) -> Vec<Peptide> {
    todo!()
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_allele_factory_generation() {
        // let allele_factory = AlleleFactory::new();
        // allele_factory.get_allele_with_motif(SOME_MOTIF);
        // allele_factory.get_allele_that_binds_kir(SOME_KIR);
    }
}
