use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::parser::*;
use crate::result::*;

use log::debug;

// TODO: Needs refactoring and error handling not finished converting from nom::Err
/// Reads a netmhcpan output file.  Not optimized to skip peptides already processed.
pub fn read_raw_netmhcpan<T>(paths: Vec<T>) -> Result<BindingData, Box<dyn std::error::Error>>
where
    T: AsRef<Path>,
{
    use RankThreshold::*;
    let mut binding_data = BindingData::new();
    #[allow(unused_mut)]
    let mut paths = paths
        .iter()
        .map(|p| p.as_ref().to_owned())
        .collect::<Vec<PathBuf>>();

    #[cfg(target_os = "windows")]
    {
        paths = paths
            .iter()
            .flat_map(|p| p.to_str())
            .flat_map(|p| glob::glob(p))
            .map(|p| p.flat_map(|p| p).collect::<Vec<PathBuf>>())
            .flatten()
            .collect::<Vec<PathBuf>>();
    }

    for path in &paths {
        let f = File::open(path)?;
        let rdr = BufReader::new(f);
        let mut netmhcpan_version = crate::NETMHCPAN_VERSION[0].to_string();

        for (n, line) in rdr.lines().enumerate() {
            debug!(
                "Parsing line number {} in netMHCpan output from file {}",
                n,
                path.display()
            );
            let line = line?;
            let (i, version_line) = is_version_line(&line).unwrap();

            if version_line.is_some() {
                let (_, version) = get_netmhcpan_version(i).unwrap();
                netmhcpan_version = version.to_string();
            } else {
                let (i, nn_line) = is_nn_line(i).unwrap();

                if nn_line.is_some() {
                    let (_, nn) = get_nn_info(i).unwrap();
                    binding_data.alleles.insert(nn);
                } else {
                    if binding_data.strong_threshold.is_none()
                        || binding_data.weak_threshold.is_none()
                    {
                        let (i, rank_line) = is_rank_line(i).unwrap();

                        if rank_line.is_some() {
                            let (_, rank) = get_rank_info(i).unwrap();
                            match rank {
                                Strong(threshold) => {
                                    binding_data.strong_threshold = Some(threshold)
                                }
                                Weak(threshold) => binding_data.weak_threshold = Some(threshold),
                            }
                        }
                    }
                    let (i, pep_line) = is_peptide_line(i).unwrap();

                    if pep_line {
                        let (i, (pos, allele, pep_seq)) = get_netmhc_entry_info(i).unwrap();
                        let (i, (alignment_mods, icore, identity)) =
                            get_netmhc_align_info(i).unwrap();

                        let protein = binding_data
                            .proteome
                            .entry(identity.to_string())
                            .or_insert_with(|| Protein::new(identity.to_string()));
                        protein
                            .add_sequence_at_pos(pos, pep_seq.to_string())
                            .unwrap();

                        let peptide = Peptide::new(
                            pos,
                            pep_seq.to_string(),
                            identity.to_string(),
                            icore.to_string(),
                            &alignment_mods,
                        );

                        binding_data.peptides.insert(peptide.clone());

                        let (_, binding_info) =
                            get_netmhc_binding_info(i, peptide, &netmhcpan_version).unwrap();

                        let allele_binding = binding_data.allele_binding.entry(allele).or_default();
                        allele_binding.push(binding_info);
                    }
                }
            }
        }
    }
    Ok(binding_data)
}

#[cfg(test)]
mod tests {
    use crate::reader::read_raw_netmhcpan;

    #[test]
    fn read_binding_protein() {
        let bd = read_raw_netmhcpan(vec!["tests/netmhcpan_wBA.txt"]).unwrap();
        let seq_expected = "TPQDLNTMLNTVGGHQAAMQMLKETINEEA";
        assert_eq!(bd.proteome.get("Gag_180_209").unwrap().seq(), seq_expected);
    }
    #[test]
    fn read_hla_bound_at_threshold() {
        let bd = read_raw_netmhcpan(vec!["tests/netmhcpan_wBA.txt"]).unwrap();
        let allele = "B27:05".parse().unwrap();

        let motif_pos = vec![2usize, 7usize, 8usize, 9usize]
            .iter()
            .map(|i| i - 1)
            .collect::<Vec<usize>>();
        let expected_motif = "PTML";

        let bound = bd.get_binding_info(&allele).unwrap();
        dbg!(&bound[0]);
        assert_eq!(bound[0].motif(&motif_pos), expected_motif.to_string());
    }
}
