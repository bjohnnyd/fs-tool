const IPD_URL: &str = "https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?";
// https://www.ebi.ac.uk/cgi-bin/ipd/kir/retrieve_ligands.cgi?C*01:02

use super::*;
use select::document::Document;
use select::predicate::Name;

enum P80 {
    I,
    T,
    N,
    K,
}

enum Bw4 {
    P80,
}

enum C1 {
    P80,
}

enum C2 {
    P80,
}

#[derive(Debug, Eq, PartialEq)]
pub(crate) enum LigandGroup {
    A3,
    A11,
    Bw4,
    Bw6,
    C1,
    C2,
    Unclassified,
    Unknown,
}


#[cfg(test)]
mod tests {
    use super::*;

}
