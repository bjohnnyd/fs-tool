// # TODO
// 1. (MAYBE) Separate peptide into own module/file
// 2. Implement Vec<PeptideDifference>.to_core(&peptide) -> String
// 3. Start working on parser

pub const MAX_PEPTIDE_LEN: usize = 12;

#[derive(Debug)]
pub struct Deletion(usize, usize);
#[derive(Debug)]
pub struct Insertion(usize, usize);

#[derive(Debug, Eq, PartialEq)]
pub enum BindLevel {
    SB,
    WB,
    NB,
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

#[derive(Debug)]
pub struct Peptide<'a> {
    pos: usize,
    length: usize,
    core: Vec<PeptideDifference>,
    interaction_core: Vec<PeptideDifference>,
    core_start_offset: usize,
    deletion: Deletion,
    insertion: Insertion,
    protein: &'a Protein,
    score: f32,
    aff: Option<f32>,
    rank: f32,
    bind_level: BindLevel,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Protein {
    identity: String,
    sequence: String,
}

impl Protein {
    fn new(identity: String) -> Self {
        Self {
            identity,
            sequence: String::new(),
        }
    }

    fn add_peptide<T>(&mut self, sequence: T, pos: usize)
    where
        T: AsRef<str>,
    {
        if self.sequence.is_empty() || pos == self.sequence.len() + 1 {
            self.sequence.push_str(sequence.as_ref())
        }
    }
}

fn pep_to_core<T: AsRef<str>>(steps: &Vec<PeptideDifference>, peptide: T) -> String {
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
    let mut peptide_chars = peptide.as_ref().chars();
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
    fn build_peptide() {
        let mut protein = Protein::new("Gag_180_209".to_string());
        let first = ("TPQDLNTMLNT", 1);
        let second = ("NTVGGHQAAMQ", 10);
        let third = ("VGGHQAAMQML", 12);

        protein.add_peptide(first.0, first.1);
        assert_eq!(
            protein,
            Protein {
                identity: String::from("Gag_180_209"),
                sequence: String::from("TPQDLNTMLNT")
            }
        );

        protein.add_peptide(second.0, second.1);
        assert_eq!(
            protein,
            Protein {
                identity: String::from("Gag_180_209"),
                sequence: String::from("TPQDLNTMLNT")
            }
        );

        protein.add_peptide(third.0, third.1);
        assert_eq!(
            protein,
            Protein {
                identity: String::from("Gag_180_209"),
                sequence: String::from("TPQDLNTMLNTVGGHQAAMQML")
            }
        );
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
