use crate::calc::LilrbScore;
use crate::cohort::CohortTemp;
use crate::error::Error;
use immunoprot::ig_like::kir::Kir;
use immunoprot::ig_like::kir_ligand::{KirLigandMap, LigandMotif};
use immunoprot::mhc::hla::ClassI;
use log::info;
use std::collections::HashMap;

pub fn read_project_ligand_info() -> Option<KirLigandMap> {
    let mut kir_ligand_map = None;
    if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
        let path = project_dir.data_dir();
        kir_ligand_map = KirLigandMap::from_path(path.join(crate::PROJECT_LIGAND_TABLE)).ok()
    }

    if kir_ligand_map.is_none() {
        info!("No local kir ligand info found. Defaulting to: https://raw.githubusercontent.com/bjohnnyd/fs-tool/dev_fs/immunoprot/src/resources/allele_motifs.tsv");
        kir_ligand_map = KirLigandMap::init().ok();
    }
    kir_ligand_map
}

pub fn read_temp_cohort<P>(path: P) -> Result<Vec<CohortTemp>, Error>
where
    P: AsRef<std::path::Path>,
{
    let mut cohort_temp = Vec::<CohortTemp>::new();

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(crate::DEFAULT_DELIM)
        .has_headers(true)
        .from_path(path)?;

    for entry in rdr.deserialize() {
        let temp: CohortTemp = entry?;
        cohort_temp.push(temp);
    }

    Ok(cohort_temp)
}

pub fn read_kir_motif_binding() -> HashMap<Kir, Vec<LigandMotif>> {
    crate::KIR_MOTIF_BINDING
        .lines()
        .fold(HashMap::new(), |mut motif_binding, line| {
            let kir_motif = line.split('\t').collect::<Vec<&str>>();

            let motifs = kir_motif[1]
                .split(';')
                .filter_map(|motif| motif.parse::<LigandMotif>().ok())
                .collect::<Vec<LigandMotif>>();

            if !motifs.is_empty() {
                if let Ok(kir) = kir_motif[0].parse::<Kir>() {
                    motif_binding.insert(kir, motifs);
                }
            }

            motif_binding
        })
}

// TODO: might be better with a proper reader and serializing
pub fn read_lilrb_scores() -> Vec<LilrbScore> {
    crate::LILRB_SIMSCORES
        .lines()
        .filter(|line| !line.contains("lilrb") && !line.starts_with('#'))
        .fold(Vec::new(), |mut scores, line| {
            let entry = line.split('\t').collect::<Vec<&str>>();

            let first_allele = entry[0].parse::<ClassI>();
            let second_allele = entry[1].parse::<ClassI>();
            let lilrb1_score = entry[2].parse::<f32>();
            let lilrb2_score = entry[3].parse::<f32>();

            if let (Ok(first_allele), Ok(second_allele), Ok(lilrb1_score), Ok(lilrb2_score)) =
                (first_allele, second_allele, lilrb1_score, lilrb2_score)
            {
                scores.push(LilrbScore {
                    first_allele,
                    second_allele,
                    lilrb1_score,
                    lilrb2_score,
                });
            }

            scores
        })
}

#[cfg(test)]
mod tests {
    use crate::io::reader::{read_kir_motif_binding, read_lilrb_scores};

    #[test]
    fn test_create_motif_binding() {
        let motif_binding = read_kir_motif_binding();
        assert_eq!(motif_binding.len(), 9);
    }
    #[test]
    fn test_create_lilrb_scores() {
        let lilrb_scores = read_lilrb_scores();
        dbg!(lilrb_scores);
    }
}
