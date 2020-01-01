use crate::netmhcpan::parser::{BindingInfo, PepInfo};
use crate::prelude::collections::{HashMap, HashSet};
/// # TODO
/// 1. Implement conversion from Vec<&str>
/// 2. Add addition of Ligand Info using some sort of matching
/// 3. Implement comparison between/Retreiving motifs as Vec of strings
use crate::prelude::fs_tool::{BindLevel, Peptide, PeptideIdentity, Proteome, HLA};
use crate::prelude::traits::FromStr;

pub struct StrongThreshold(f32);
pub struct WeakThreshold(f32);
pub struct Distance(f32);

#[derive(Debug)]
pub struct NearestNeighbour {
    index: HLA,
    distance: f32,
    nearest_neighbour: HLA,
}

impl PartialEq for NearestNeighbour {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index && self.nearest_neighbour == other.nearest_neighbour
    }
}

impl std::hash::Hash for NearestNeighbour {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.index.hash(state);
        self.nearest_neighbour.hash(state);
    }
}

impl Eq for NearestNeighbour {}

impl From<(HLA, f32, HLA)> for NearestNeighbour {
    fn from(t: (HLA, f32, HLA)) -> Self {
        Self {
            index: t.0,
            distance: t.1,
            nearest_neighbour: t.2,
        }
    }
}
#[derive(Debug)]
pub struct NetMHCpanRecord {
    peptide_identity: PeptideIdentity,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

impl NetMHCpanRecord {
    fn new(
        peptide_identity: PeptideIdentity,
        bind_data: (f32, Option<f32>, f32, BindLevel),
    ) -> Self {
        let (score, aff, rank, bind_level) = bind_data;

        Self {
            peptide_identity,
            score,
            aff,
            rank,
            bind_level,
        }
    }

    fn bound(&self) -> bool {
        match self.bind_level {
            BindLevel::NB => false,
            _ => true,
        }
    }

    fn bound_by_rank(&self, threshold: Option<f32>) -> bool {
        if let Some(threshold) = threshold {
            self.rank < threshold
        } else {
            self.bound()
        }
    }
}

#[derive(Debug)]
pub struct NetMHCpanSummary {
    pub alleles: HashSet<NearestNeighbour>,
    records: HashMap<HLA, Vec<NetMHCpanRecord>>,
    proteome: Proteome,
    peptides: HashMap<PeptideIdentity, Peptide>,
    pub weak_threshold: Option<f32>,
    pub strong_threshold: Option<f32>,
}

impl NetMHCpanSummary {
    pub(crate) fn new() -> Self {
        Self {
            alleles: HashSet::<NearestNeighbour>::new(),
            records: HashMap::<HLA, Vec<NetMHCpanRecord>>::new(),
            proteome: Proteome::new(),
            peptides: HashMap::<PeptideIdentity, Peptide>::new(),
            weak_threshold: None,
            strong_threshold: None,
        }
    }

    pub fn add_hla(&mut self, nn: NearestNeighbour) -> bool {
        self.alleles.insert(nn)
    }

    pub fn add_peptide(&mut self, pep: Peptide) -> Option<Peptide> {
        let peptide_identity = PeptideIdentity::from(&pep);
        self.peptides.insert(peptide_identity, pep)
    }

    pub fn add_sequence<T: AsRef<str>>(&mut self, id: T, pos: usize, seq: T) {
        self.proteome
            .add_peptide(id.as_ref().to_string(), pos, seq.as_ref().to_string())
    }

    pub fn is_threshold_set(&self) -> bool {
        self.weak_threshold.is_some() && self.strong_threshold.is_some()
    }

    pub fn add_threshold(&mut self, rank_type: &str, threshold: f32) {
        if rank_type.trim() == "Strong" {
            self.strong_threshold = Some(threshold);
        } else if rank_type.trim() == "Weak" {
            self.weak_threshold = Some(threshold)
        }
    }

    /* Need to deal with error */
    pub fn insert_hla_record(
        &mut self,
        peptide_identity: PeptideIdentity,
        BindingInfo(hla_id, score, aff, rank, bind_level): BindingInfo,
    ) {
        let hla = hla_id.parse::<HLA>().expect("could not parse hla");

        let netmhcpan_record =
            NetMHCpanRecord::new(peptide_identity, (score, aff, rank, bind_level));

        let records = self.records.entry(hla).or_insert(Vec::new());
        records.push(netmhcpan_record);
    }

    pub fn get_bound(&self, hla: &HLA, threshold: Option<f32>) -> Vec<&Peptide> {
        let mut peptides = Vec::<&Peptide>::new();

        if let Some(records) = self.records.get(hla) {
            records
                .iter()
                .filter(|record| record.bound_by_rank(threshold))
                .for_each(|record| {
                    if let Some(peptide) = self.peptides.get(&record.peptide_identity) {
                        peptides.push(peptide)
                    }
                });
        }
        peptides
    }

    pub fn get_motifs(&self, peptides: &Vec<&Peptide>, motif_pos: Vec<usize>) -> Vec<String> {
        let mut slices = Vec::<String>::new();
        peptides.iter().for_each(|pep| {
            if let Some(protein_seq) = self.proteome.proteins.get(&pep.protein) {
                let slice = protein_seq.as_bytes()[pep.pos..pep.pos + pep.length]
                    .iter()
                    .enumerate()
                    .filter(|(pep_pos, aa)| motif_pos.contains(&(pep_pos + 1)))
                    .map(|(_, aa)| *aa as char)
                    .collect::<String>();

                slices.push(slice);
            }
        });

        slices.sort();
        slices
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::fs_tool::*;

    const netmhcpan: &str =
        "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\
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
                ";

    #[test]
    fn test_get_motifs() {
        let mut netmhcpan_summary = read_netmhcpan(&netmhcpan).unwrap();
        let hla = "HLA-A03:01".parse::<HLA>().unwrap();

        let cd8 = vec![2, 3, 4, 5, 6, 9];
        let hla_bound = netmhcpan_summary.get_bound(&hla, Some(10_f32));
        dbg!(netmhcpan_summary.get_motifs(&hla_bound, cd8));
    }
}
