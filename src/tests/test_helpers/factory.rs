use immunoprot::ig_like::kir_ligand::KirLigandMap;
use netmhcpan::result::{BindingInfo, Peptide};
use rand::{distributions::Alphanumeric, prelude::SliceRandom, Rng};

const AA: [&str; 20] = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
    "V",
];

const IDENTITY: &str = "Protein";
const ALIGNMENT_MODS: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

pub struct AlleleFactory {
    kir_ligand_map: KirLigandMap,
}

pub enum RankThreshold {
    Strong,
    Weak,
}

pub struct BindingFactory {
    peptides: Vec<Peptide>,
}

impl BindingFactory {
    pub fn new(peptides: Vec<Peptide>) -> Self {
        Self { peptides }
    }
    pub fn generate_binding_info(
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

pub fn generate_peptides(n: usize, length: usize) -> Vec<Peptide> {
    let mut aa = AA.clone();
    let mut rng = rand::thread_rng();

    (0..n)
        .map(|_| {
            aa.shuffle(&mut rng);
            let seq = aa.iter().cloned().take(length).collect::<String>();
            Peptide::new(
                0,
                seq.clone(),
                IDENTITY.to_string(),
                seq.to_string(),
                &ALIGNMENT_MODS,
            )
        })
        .collect()
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_peptide_generation() {
        let peptides = generate_peptides(3, 9);
        assert_eq!(peptides.len(), 3);
        assert!(peptides.iter().all(|peptide| peptide.len() == 9));
    }

    #[test]
    fn test_binding_info_generation() {
        let peptides = generate_peptides(5, 9);
        let binding_factory = BindingFactory::new(peptides);

        let binding_infos = binding_factory.generate_binding_info(10, Some(RankThreshold::Strong));

        assert!(binding_infos
            .iter()
            .all(|binding_info| binding_info.rank < 0.5));
    }

    #[test]
    fn test_allele_factory_generation() {
        // let allele_factory = AlleleFactory::new();
        // allele_factory.get_allele_with_motif(SOME_MOTIF);
        // allele_factory.get_allele_that_binds_kir(SOME_KIR);
    }
}
