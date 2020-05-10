use crate::cohort::CohortTemp;
use crate::error::Error;
use immunoprot::ig_like::kir_ligand::KirLigandMap;
use log::warn;

pub fn read_project_ligand_info() -> Option<KirLigandMap> {
    let mut kir_ligand_map = None;
    if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
        let path = project_dir.data_dir();
        kir_ligand_map = KirLigandMap::from_path(path.join(crate::PROJECT_LIGAND_TABLE)).ok()
    }

    if kir_ligand_map.is_none() {
        warn!("No local kir ligand info found. Defaulting to: https://raw.githubusercontent.com/bjohnnyd/fs-tool/dev_fs/immunoprot/src/resources/allele_motifs.tsv");
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
