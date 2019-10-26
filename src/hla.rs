use nom::{
    IResult,
    bytes::complete::tag,
};

pub mod mhc {
    use super::*;
//    fn parse_hla_name(hla_name: &str) -> IResult<&str, &str> {
//
//    }

    pub struct HLA {
        gene: Locus,
        allele_group: u8,
        hla_protein: u8,
        cds_synonymous_sub: u8,
        non_coding_difference: u8,
        expression_change: ExpressionChange,
        ligand_group: LigandGroup,
        mhc_class: MHC,
    }

//    impl HLA {
//        pub fn new(name: String) -> HLA {}
//    }

    enum Locus {
        A,
        B,
        C,
        DP,
        DM,
        DO,
        DQ,
        DR,
    }

    enum MHC {
        I,
        II,
    }

    enum ExpressionChange {
        N,
        L,
        S,
        C,
        A,
        Q,
    }

    pub enum LigandGroup {
        I,
        II,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

//    #[test]
//    fn test_create_from_short_name() {
//        assert_eq!("HLA-A*03:01", )
//    }
}
