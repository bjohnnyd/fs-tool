use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::error::Error;
use crate::parser::*;
use crate::result::*;

pub fn read_raw_netmhcpan<T>(path: T) -> Result<BindingData, Error>
where
    T: AsRef<Path>,
{
    use RankThreshold::*;

    let mut f = File::open(path)?;
    let mut raw_data = String::new();
    f.read_to_string(&mut raw_data).unwrap();

    let mut sequence_info_finished = false;

    Ok(raw_data
        .lines()
        .fold(BindingData::new(), |mut binding_data, line| {
            let (i, nn_line) = is_nn_line(line.as_bytes()).unwrap();

            if let Some(nn_line) = nn_line {
                let (_, nn) = get_nn_info(i).unwrap();
                binding_data.alleles.insert(nn);
            } else {
                if binding_data.strong_threshold.is_none() || binding_data.weak_threshold.is_none()
                {
                    let (i, rank_line) = is_rank_line(i).unwrap();

                    if let Some(rank_line) = rank_line {
                        let (_, rank) = get_rank_info(i).unwrap();
                        match rank {
                            Strong(threshold) => binding_data.strong_threshold = Some(threshold),
                            Weak(threshold) => binding_data.weak_threshold = Some(threshold),
                            _ => (),
                        }
                    }
                }
                let (i, pep_line) = is_peptide_line(i).unwrap();

                if pep_line {
                    let (i, (pos, allele, pep_seq)) = get_netmhc_entry_info(i).unwrap();
                    let (i, (alignment_mods, icore, identity)) = get_netmhc_align_info(i).unwrap();

                    let mut protein = binding_data
                        .proteome
                        .entry(identity.clone())
                        .or_insert_with(|| Protein::new(identity.clone()));
                    protein.add_sequence_at_pos(pos, pep_seq.clone());

                    let peptide = Peptide::new(
                        pos,
                        pep_seq.len(),
                        pep_seq,
                        identity.clone(),
                        icore,
                        &alignment_mods,
                    );

                    binding_data.peptides.insert(peptide.clone());

                    let (_, binding_info) = get_netmhc_binding_info(i, peptide).unwrap();

                    let allele_binding = binding_data.allele_binding.entry(allele).or_default();
                    allele_binding.push(binding_info);
                }
            }

            binding_data
        }))
}

#[cfg(test)]
mod tests {
    use crate::reader::read_raw_netmhcpan;

    #[test]
    fn read_binding_data_raw() {
        let bd = read_raw_netmhcpan("tests/netmhcpan_wBA.txt").unwrap();
        dbg!(&bd);
    }
}
