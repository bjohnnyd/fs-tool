use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use crate::error::Error;
use crate::parser::*;
use crate::result::*;

// TODO: Needs serious refactoring and error handling not finished covnerting from nom::Err
/// Reads a netmhcpan output file.  Not optimized to skip peptides already processed.
pub fn read_raw_netmhcpan<T>(path: T) -> Result<BindingData, Box<dyn std::error::Error>>
where
    T: AsRef<Path>,
{
    use RankThreshold::*;

    let mut f = File::open(path)?;
    // let mut raw_data = String::new();
    let rdr = BufReader::new(f);
    // f.read_to_string(&mut raw_data).unwrap();

    let mut binding_data = BindingData::new();

    for line in rdr.lines() {
        let line = line?;
        let (i, nn_line) = is_nn_line(line.as_str()).unwrap();

        if let Some(nn_line) = nn_line {
            let (_, nn) = get_nn_info(i).unwrap();
            binding_data.alleles.insert(nn);
        } else {
            if binding_data.strong_threshold.is_none() || binding_data.weak_threshold.is_none() {
                let (i, rank_line) = is_rank_line(i).unwrap();

                if let Some(rank_line) = rank_line {
                    let (_, rank) = get_rank_info(i).unwrap();
                    match rank {
                        Strong(threshold) => binding_data.strong_threshold = Some(threshold),
                        Weak(threshold) => binding_data.weak_threshold = Some(threshold),
                    }
                }
            }
            let (i, pep_line) = is_peptide_line(i).unwrap();

            if pep_line {
                let (i, (pos, allele, pep_seq)) = get_netmhc_entry_info(i).unwrap();
                let (i, (alignment_mods, icore, identity)) = get_netmhc_align_info(i).unwrap();

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

                let (_, binding_info) = get_netmhc_binding_info(i, peptide).unwrap();

                let allele_binding = binding_data.allele_binding.entry(allele).or_default();
                allele_binding.push(binding_info);
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
        let bd = read_raw_netmhcpan("tests/netmhcpan_wBA.txt").unwrap();
        let seq_expected = "TPQDLNTMLNTVGGHQAAMQMLKETINEEA";
        assert_eq!(bd.proteome.get("Gag_180_209").unwrap().seq(), seq_expected);
    }
    #[test]
    fn read_hla_bound_at_threshold() {
        use crate::result::WEAK_BOUND_THRESHOLD;

        let bd = read_raw_netmhcpan("tests/netmhcpan_wBA.txt").unwrap();
        let allele = "B27:05".parse().unwrap();

        let motif_pos = vec![2, 3, 4, 5, 7, 8];
        let expected_motif = "AAMQLK";

        let bound = bd.get_bound_info(&allele, WEAK_BOUND_THRESHOLD);
        assert_eq!(bound[0].motif(&motif_pos), expected_motif.to_string());
    }
}
