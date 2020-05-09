use immunoprot::ig_like::kir_ligand::KirLigandMap;
use log::warn;

pub(crate) fn read_project_ligand_info() -> Option<KirLigandMap> {
    let mut kir_ligand_map = None;
    if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
        let path = project_dir.data_dir();
        kir_ligand_map = KirLigandMap::from_path(path.join(crate::PROJECT_LIGAND_TABLE)).ok()
    }

    if kir_ligand_map.is_none() {
        warn!("No local kir ligand info found. Defaulting to: https://raw.githubusercontent.com/bjohnnyd/fs-tool/dev_fs/immunoprot/src/resources/2019-12-29_lg.tsv");
        kir_ligand_map = KirLigandMap::init().ok();
    }
    kir_ligand_map
}
