pub mod error;

pub mod mhc {
    pub mod hla;
    pub mod mhc_I;
}

pub mod ig_like {
    pub mod kir;
    pub mod lilrb;
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_ligand_info() {
        let lg_info = include_str!("resources/2019-12-29_lg.tsv");

        lg_info
            .lines()
            .for_each(|l|{
                dbg!(l);
            });
    }
}
