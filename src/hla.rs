mod mhc {

    pub struct HLA {
        gene: Locus,
        allele_group: u8,
        hla_protein: u8,
        cds_synonymous_sub: u8,
        non_coding_difference: u8,
        expression_change: ExpressionChange,
        ligand_group: LigandGroup,
        mhc_class: MHC
    }

    pub enum Locus {
        A,
        B,
        C,
        DP,
        DM,
        DO,
        DQ,
        DR
    }

    pub enum MHC {
        I,
        II
    }

    pub enum ExpressionChange {
        N,
        L,
        S,
        C,
        A,
        Q,
    }
    pub enum LigandGroup {
        I,
        II
    }
}
