use immunoprot::ig_like::kir_ligand::{KirLigandMap, IPD_KIR_URL};
use immunoprot::mhc::hla::ClassI;
use log::warn;

/// Writes ligand information to the global file
pub(crate) fn write_project_ligand_info(kir_ligand: &KirLigandMap) {
    let dt_now = chrono::Utc::now().format("%Y-%m-%d");
    let comment_line = format!("# Obtained from {} on {}", IPD_KIR_URL, dt_now);

    let mut alleles = kir_ligand.alleles.iter().collect::<Vec<&ClassI>>();
    alleles.sort();

    if let Some(project_dir) = directories::ProjectDirs::from("", "", "fstool") {
        if !project_dir.data_dir().exists() {
            std::fs::create_dir(project_dir.data_dir()).unwrap_or_else(|_| {
                warn!(
                    "Could not create project cache directory {}",
                    project_dir.data_dir().display()
                )
            });
        }

        let path = project_dir.data_dir().join(crate::PROJECT_LIGAND_TABLE);

        let wrt = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(path.to_path_buf());

        if let Ok(mut wrt) = wrt {
            match wrt.write_record(&vec![comment_line, "".to_string(), "".to_string()]) {
                Ok(_) => (),
                Err(_) => warn!("Could not write download date time information to ligand file"),
            }

            alleles.iter().for_each(|allele| {
                if let Some(info) = kir_ligand.cache.get(allele) {
                    let record = vec![
                        info.allele().to_string(),
                        info.motif().to_string(),
                        info.freq().to_string(),
                    ];

                    match wrt.write_record(&record) {
                        Ok(_) => (),
                        Err(_) => warn!(
                            "Could not write kir ligand information for allele {}",
                            allele.to_nomenclature_string()
                        ),
                    }
                }
            })
        } else {
            warn!(
                "Could not open kir ligand information file located at {}",
                path.display()
            )
        }
    } else {
        warn!("Could not access global project directory for this OS")
    }
}
