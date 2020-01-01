use crate::prelude::fs_tool::{BindingInfo, PepInfo, VariantInfo};
use nom::lib::std::collections::HashMap;

pub const MAX_PEPTIDE_LEN: usize = 12;

#[derive(Debug, PartialEq, Eq)]
pub struct Deletion(pub usize, pub usize);
#[derive(Debug, PartialEq, Eq)]
pub struct Insertion(pub usize, pub usize);

#[derive(Debug, Eq, PartialEq)]
pub enum BindLevel {
    SB,
    WB,
    NB,
}

impl From<Option<&str>> for BindLevel {
    fn from(bind_level: Option<&str>) -> Self {
        match bind_level {
            Some("WB") => BindLevel::WB,
            Some("SB") => BindLevel::SB,
            Some("NB") => BindLevel::NB,
            _ => BindLevel::NB,
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
pub enum PeptideDifference {
    Match(usize),
    Deletion(usize),
    Insertion(usize),
}

impl PeptideDifference {
    fn kind(&self) -> usize {
        match self {
            PeptideDifference::Match(_) => 0,
            PeptideDifference::Deletion(_) => 1,
            PeptideDifference::Insertion(_) => 2,
        }
    }
}

impl std::ops::AddAssign<usize> for PeptideDifference {
    fn add_assign(&mut self, rhs: usize) {
        match self {
            PeptideDifference::Match(val) => val.add_assign(rhs),
            PeptideDifference::Deletion(val) => val.add_assign(rhs),
            PeptideDifference::Insertion(val) => val.add_assign(rhs),
        }
    }
}

#[derive(Debug, PartialEq, Eq, Hash)]
pub struct PeptideIdentity {
    pub pos: usize,
    pub identity: String,
}

impl<'a> From<&'a Peptide> for PeptideIdentity {
    fn from(pep: &'a Peptide) -> Self {
        let pos = pep.pos;
        let identity = pep.protein.to_string();

        Self { pos, identity }
    }
}

impl<'a> From<&'a PepInfo<'a>> for PeptideIdentity {
    fn from(pep_info: &'a PepInfo<'a>) -> Self {
        let pos = pep_info.0;
        let identity = pep_info.4.to_string();

        Self { pos, identity }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Peptide {
    pub pos: usize,
    pub length: usize,
    core: Vec<PeptideDifference>,
    interaction_core: Vec<PeptideDifference>,
    core_start_offset: usize,
    deletion: Deletion,
    insertion: Insertion,
    pub protein: String,
}

impl<'a> From<(PepInfo<'a>, VariantInfo)> for Peptide {
    fn from(info: (PepInfo<'a>, VariantInfo)) -> Self {
        let PepInfo(pos, pep_seq, core_seq, icore_seq, protein) = info.0;
        let VariantInfo(core_start_offset, del_gp, del_gl, ins_gp, ins_gl) = info.1;

        let core = edit_distance(pep_seq, core_seq);
        let interaction_core = edit_distance(pep_seq, icore_seq);
        let length = pep_seq.len();
        let deletion = Deletion(del_gp, del_gl);
        let insertion = Insertion(ins_gp, ins_gl);

        Self {
            pos,
            length,
            core,
            interaction_core,
            core_start_offset,
            deletion,
            insertion,
            protein: protein.to_string(),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
pub struct Proteome {
    pub proteins: HashMap<String, String>,
}

impl Proteome {
    pub(crate) fn new() -> Self {
        let proteins = HashMap::<String, String>::new();
        Self { proteins }
    }

    pub(crate) fn add_peptide<T>(&mut self, identity: T, pos: usize, pep: T)
    where
        T: AsRef<str>,
    {
        let sequence = self
            .proteins
            .entry(identity.as_ref().to_string())
            .or_insert(String::new());
        if sequence.is_empty() {
            sequence.push_str(pep.as_ref());
        } else {
            if let Some(end) = pep.as_ref().chars().last() {
                sequence.push(end);
            }
        }
    }
}

fn pep_to_core<T: AsRef<str>>(steps: &[PeptideDifference], peptide: T) -> String {
    let mut peptide_seq = peptide.as_ref().chars();
    let mut core = String::new();
    steps.iter().for_each(|edit| match edit {
        PeptideDifference::Match(aa) => {
            core.push_str(peptide_seq.by_ref().take(*aa).collect::<String>().as_ref())
        }
        PeptideDifference::Deletion(aa) => {
            peptide_seq.by_ref().take(*aa).count();
        }
        PeptideDifference::Insertion(aa) => (0..*aa).for_each(|_| core.push('-')),
    });
    core
}

fn edit_distance<T: AsRef<str>>(peptide: T, core: T) -> Vec<PeptideDifference> {
    let mut core_chars = core.as_ref().chars().peekable();
    let peptide_chars = peptide.as_ref().chars();
    let mut result = Vec::<PeptideDifference>::with_capacity(MAX_PEPTIDE_LEN);

    for pep_aa in peptide_chars {
        if let Some(next_core_aa) = core_chars.peek() {
            let next_diff = match next_core_aa {
                _c if *_c == pep_aa => {
                    core_chars.next();
                    PeptideDifference::Match(1)
                }
                '-' => PeptideDifference::Insertion(1),
                _ => PeptideDifference::Deletion(1),
            };

            match result.last_mut() {
                None => result.push(next_diff),
                Some(current_diff) => {
                    if current_diff.kind() == next_diff.kind() {
                        *current_diff += 1;
                    } else {
                        result.push(next_diff)
                    }
                }
            }
        }
    }

    result
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_proteome() {
        let pep1 = (1, "TPQDLNTMLNT");
        let pep2 = (4, "DLNTMLNTVGG");
        let pep3 = (12, "VGGHQAAMQML");
        let identity = "Gag_180_209";

        let mut base_hashmap = HashMap::new();
        base_hashmap.insert(identity.to_string(), pep1.1.to_string());
        let expected = Proteome {
            proteins: base_hashmap.clone(),
        };

        let mut proteome = Proteome::new();

        proteome.add_peptide(identity, pep1.0, pep1.1);
        assert_eq!(proteome, expected);

        proteome.add_peptide(identity, pep2.0, pep2.1);
        assert_eq!(proteome, expected);

        let expected_seq = base_hashmap
            .entry(identity.to_string())
            .or_insert(String::new());
        expected_seq.push_str(pep3.1);
        let expected = Proteome {
            proteins: base_hashmap,
        };

        proteome.add_peptide(identity, pep3.0, pep3.1);
        assert_eq!(proteome, expected);
    }

    #[test]
    fn calculate_edits() {
        let peptide = "TVGGHQAAMQM";
        let core = "TVHQAAMQM";

        let result = edit_distance(peptide, core);

        let expected = vec![
            PeptideDifference::Match(2),
            PeptideDifference::Deletion(2),
            PeptideDifference::Match(7),
        ];

        assert_eq!(result, expected);
    }

    #[test]
    fn get_core() {
        let peptide = "TVGGHQAAMQM";
        let core = "TVHQAAMQM";

        let result = edit_distance(peptide, core);

        assert_eq!(pep_to_core(&result, peptide), core);
    }

    #[test]
    fn test_netmhcpan_line() {
        let example_wo_ba = vec![
            "3",
            "HLA-A*03:01",
            "QDLNTMLNTVG",
            "QDLNTNTVG",
            "0",
            "5",
            "2",
            "0",
            "0",
            "QDLNTMLNTVG",
            "Gag_180_209",
            "0.0000040",
            "96.0000",
        ];
    }
}
