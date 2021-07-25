use crate::tests::test_helpers::*;
use rand::{distributions::Alphanumeric, prelude::SliceRandom, Rng};
use std::collections::hash_map::RandomState;
use std::collections::{HashMap, HashSet};
use std::iter::FromIterator;

const AA: [&str; 20] = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
    "V",
];

const IDENTITY: &str = "Protein";
const ALIGNMENT_MODS: [usize; 12] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];

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

    /// Generates n number of binding info data of a given rank if provided
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

/// Generate n peptides of specified length using 20 AA alphabet
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

pub struct AlleleFactory {
    kir_ligand_map: &'static KirLigandMap,
    kir_binding: HashMap<Kir, Vec<LigandMotif>>,
}

impl AlleleFactory {
    pub fn new() -> TestResult<Self> {
        let kir_ligand_map = &KIR_MAP;
        let kir_binding = read_kir_motif_binding();
        Ok(Self {
            kir_ligand_map,
            kir_binding,
        })
    }

    pub fn gen_allele_with_motif(&self, motif: &LigandMotif, n: usize) -> Vec<ClassI> {
        let mut rng = rand::thread_rng();
        let mut alleles_with_motif = self
            .kir_ligand_map
            .cache
            .iter()
            .filter(|(_, kir_ligand_info)| kir_ligand_info.motif() == motif)
            .map(|(k, v)| k)
            .collect::<Vec<&ClassI>>();

        alleles_with_motif.shuffle(&mut rng);
        alleles_with_motif.into_iter().take(n).cloned().collect()
    }

    pub fn gen_allele_with_kir_bound(&self, kir: &Kir, n: usize) -> TestResult<Vec<ClassI>> {
        let mut rng = rand::thread_rng();

        let target_motifs: HashSet<&LigandMotif, RandomState> = HashSet::from_iter(
            self.kir_binding
                .get(kir)
                .ok_or(format!("No motifs found for {}", kir))?,
        );

        let mut alleles_that_bind_kir = self
            .kir_ligand_map
            .cache
            .iter()
            .filter(|(_, kir_ligand_info)| target_motifs.contains(kir_ligand_info.motif()))
            .map(|(k, v)| k)
            .collect::<Vec<&ClassI>>();

        alleles_that_bind_kir.shuffle(&mut rng);
        Ok(alleles_that_bind_kir.into_iter().take(n).cloned().collect())
    }
}

#[cfg(test)]
mod tests {

    use std::str::FromStr;

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
    fn test_allele_motif_factory_generation() {
        use LigandMotif::*;
        let allele_factory = AlleleFactory::new().unwrap();
        let a11_alleles = allele_factory.gen_allele_with_motif(&A11, 5);

        assert_eq!(a11_alleles.len(), 5);
        assert!(a11_alleles.iter().all(|allele| {
            let kir_ligand_info = allele_factory.kir_ligand_map.get_allele_info(allele);
            let motif = kir_ligand_info[0].motif();

            *motif == A11
        }))
    }
    #[test]
    fn test_allele_kir_bound_factory_generation() {
        let kir = Kir::from_str("KIR2DL2").unwrap();
        let allele_factory = AlleleFactory::new().unwrap();
        let kir2dl2_alleles = allele_factory.gen_allele_with_kir_bound(&kir, 5).unwrap();

        assert_eq!(kir2dl2_alleles.len(), 5);
        assert!(kir2dl2_alleles.iter().all(|allele| {
            let kir_ligand_info = allele_factory.kir_ligand_map.get_allele_info(allele);
            let motif = kir_ligand_info[0].motif();
            let kir_bound = get_bound_kirs(&allele_factory.kir_binding, &motif);

            kir_bound.iter().any(|bound_kir| *bound_kir == kir)
        }))
    }
}
