use crate::calc::{CalcFsResult, CohortResult};
use crate::error::Error;
use crate::meta::{AlleleMeta, BindingMeta};

use immunoprot::ig_like::kir_ligand::{KirLigandMap, IPD_KIR_URL};
use immunoprot::mhc::hla::ClassI;
use log::warn;

/// Writes ligand information to the global file
pub(crate) fn write_project_ligand_info(kir_ligand: &KirLigandMap) {
    let dt_now = chrono::Utc::now().format("%Y-%m-%d");
    let comment_line = format!("# Obtained from {} on {}", IPD_KIR_URL, dt_now);

    let mut alleles = kir_ligand.alleles.iter().collect::<Vec<&ClassI>>();
    alleles.sort();

    if let Some(project_dir) = directories::ProjectDirs::from("", "", crate::TOOL_NAME) {
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

pub struct OutputWriters {
    pub allele_meta: csv::Writer<std::fs::File>,
    pub binding_meta: csv::Writer<std::fs::File>,
    pub allele_fs_result: csv::Writer<std::fs::File>,
    pub cohort_result: Option<csv::Writer<std::fs::File>>,
}

impl OutputWriters {
    pub fn write_allele_meta(
        &mut self,
        metadata: &[AlleleMeta],
    ) -> std::result::Result<Vec<()>, Error> {
        let write_result = metadata
            .iter()
            .map(|meta| self.allele_meta.serialize(meta))
            .collect::<Result<Vec<_>, _>>();

        Ok(write_result.or_else(|_| Err(Error::CouldNotWriteAlleleMeta))?)
    }

    pub fn write_binding_meta(
        &mut self,
        metadata: &[BindingMeta],
    ) -> std::result::Result<Vec<()>, Error> {
        let write_result = metadata
            .iter()
            .map(|meta| self.binding_meta.serialize(meta))
            .collect::<Result<Vec<_>, _>>();

        Ok(write_result.or_else(|_| Err(Error::CouldNotWriteBindingMeta))?)
    }

    pub fn write_fs_result(
        &mut self,
        fs_results: &[CalcFsResult],
    ) -> std::result::Result<Vec<()>, Error> {
        let write_result = fs_results
            .iter()
            .map(|result| self.allele_fs_result.serialize(result))
            .collect::<Result<Vec<_>, _>>();

        Ok(write_result.or_else(|_| Err(Error::CouldNotWriteFsResult))?)
    }

    pub fn write_cohort_result(
        &mut self,
        cohort_results: &[CohortResult],
    ) -> std::result::Result<Vec<()>, Error> {
        if let Some(ref mut cohort_result_file) = self.cohort_result {
            let write_result = cohort_results
                .iter()
                .map(|result| cohort_result_file.serialize(result))
                .collect::<Result<Vec<_>, _>>();
            Ok(write_result.or_else(|_| Err(Error::CouldNotWriteFsResult))?)
        } else {
            Ok(Vec::new())
        }
    }
}
