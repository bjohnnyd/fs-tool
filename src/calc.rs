use crate::prelude::traits::FromStr;
#[derive(Debug, PartialEq, Eq)]
pub struct Measure {
    pub name: String,
    pub pos: Vec<usize>,
}

impl FromStr for Measure {
    type Err = &'static str;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let name_pos = s.split(':').collect::<Vec<&str>>();
        let mut name = String::new();
        let mut pos = Vec::<usize>::new();

        if name_pos.len() < 2 {
            return Err("Incorrect measure speicifier");
        } else {
            let name = name_pos[0].to_string();
            let pos = name_pos[1]
                .split(',')
                .filter_map(|digit| digit.parse::<usize>().ok())
                .collect::<Vec<usize>>();
            Ok(Self { name, pos })
        }
    }
}
pub fn intersection_count_sorted_motifs(a: &[String], b: &[String]) -> u32 {
    let mut count = 0;
    let mut b_iter = b.iter().map(String::as_bytes);

    if let Some(mut current_b) = b_iter.next() {
        for current_a in a.iter().map(String::as_bytes) {
            while current_b < current_a {
                current_b = match b_iter.next() {
                    Some(current_b) => current_b,
                    None => return count,
                };
            }
            if current_a == current_b {
                count += 1;
            }
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::netmhcpan::parser::{BindingInfo, PepInfo};
    use crate::prelude::collections::{HashMap, HashSet};
    use crate::prelude::fs_tool::*;
    use crate::prelude::traits::FromStr;

    const netmhcpan: &str =
        "HLA-A03:01 : Distance to training data  0.000 (using nearest neighbor HLA-A03:01)\n\
            \n\
            # Rank Threshold for Strong binding peptides   0.500\n\
            # Rank Threshold for Weak binding peptides   2.000\n\
            -----------------------------------------------------------------------------------\n\
              Pos          HLA         Peptide       Core Of Gp Gl Ip Il        Icore        Identity     Score   %Rank  BindLevel\n\
            -----------------------------------------------------------------------------------\n\
                1  HLA-A*03:01     TPQDLNTMLNT  TPQDTMLNT  0  4  2  0  0  TPQDLNTMLNT     Gag_180_209 0.0000160 83.3333\n\
                2  HLA-A*03:01     PQDLNTMLNTV  PQDLNTMLV  0  8  2  0  0  PQDLNTMLNTV     Gag_180_209 0.0000120 87.0000\n\
                3  HLA-A*03:01     QDLNTMLNTVG  QDLNTNTVG  0  5  2  0  0  QDLNTMLNTVG     Gag_180_209 0.0000040 96.0000\n\
                4  HLA-A*03:01     DLNTMLNTVGG  DLNTNTVGG  0  4  2  0  0  DLNTMLNTVGG     Gag_180_209 0.0000040 96.0000\n\
                5  HLA-A*03:01     LNTMLNTVGGH  LMLNTVGGH  0  1  2  0  0  LNTMLNTVGGH     Gag_180_209 0.0001090 48.8333\n\
                6  HLA-A*03:01     NTMLNTVGGHQ  NTMTVGGHQ  0  3  2  0  0  NTMLNTVGGHQ     Gag_180_209 0.0001260 46.1429\n\
                7  HLA-A*03:01     TMLNTVGGHQA  TMLNGGHQA  0  4  2  0  0  TMLNTVGGHQA     Gag_180_209 0.0001260 46.1429\n\
                8  HLA-A*03:01     MLNTVGGHQAA  MLNGGHQAA  0  3  2  0  0  MLNTVGGHQAA     Gag_180_209 0.0002300 36.4000\n\
                9  HLA-A*03:01     LNTVGGHQAAM  NTVGGHQAM  1  7  1  0  0   NTVGGHQAAM     Gag_180_209 0.0000530 62.5000\n\
               10  HLA-A*03:01     NTVGGHQAAMQ  NTVGGAAMQ  0  5  2  0  0  NTVGGHQAAMQ     Gag_180_209 0.0001420 44.1250\n\
               11  HLA-A*03:01     TVGGHQAAMQM  TVHQAAMQM  0  2  2  0  0  TVGGHQAAMQM     Gag_180_209 0.0004120 28.5000\n\
               12  HLA-A*03:01     VGGHQAAMQML  VGGAAMQML  0  3  2  0  0  VGGHQAAMQML     Gag_180_209 0.0000120 87.0000\n\
               13  HLA-A*03:01     GGHQAAMQMLK  GQAAMQMLK  0  1  2  0  0  GGHQAAMQMLK     Gag_180_209 0.0313010  3.8215\n\
               14  HLA-A*03:01     GHQAAMQMLKE  GQAAMQMLK  0  1  1  0  0   GHQAAMQMLK     Gag_180_209 0.0004080 28.6176\n\
               15  HLA-A*03:01     HQAAMQMLKET  HQAAMQMLK  0  0  0  0  0    HQAAMQMLK     Gag_180_209 0.0003110 32.0000\n\
               16  HLA-A*03:01     QAAMQMLKETI  QAAMQMLTI  0  7  2  0  0  QAAMQMLKETI     Gag_180_209 0.0000140 85.0000\n\
               17  HLA-A*03:01     AAMQMLKETIN  AAMQKETIN  0  4  2  0  0  AAMQMLKETIN     Gag_180_209 0.0000060 93.7500\n\
               18  HLA-A*03:01     AMQMLKETINE  AMLKETINE  0  1  2  0  0  AMQMLKETINE     Gag_180_209 0.0001620 41.9000\n\
               19  HLA-A*03:01     MQMLKETINEE  MLKETINEE  0  1  2  0  0  MQMLKETINEE     Gag_180_209 0.0000850 53.5417\n\
               20  HLA-A*03:01     QMLKETINEEA  QMLKETINA  0  8  2  0  0  QMLKETINEEA     Gag_180_209 0.0000410 67.2727\n\
               21  HLA-A*03:01       MLKETINEE  MLKETINEE  0  0  0  0  0    MLKETINEE     Gag_180_209 0.0079270  7.4157\n\
               22  HLA-A*03:01       LKETINEEA  LKETINEEA  0  0  0  0  0    LKETINEEA     Gag_180_209 0.0000450 65.4545\n\
                ";

    #[test]
    fn test_measure() {
        let input_measure = "CD8:2,3,4,5,6,9";
        let measure = input_measure.parse::<Measure>().unwrap();

        assert_eq!(
            measure,
            Measure {
                name: "CD8".to_string(),
                pos: vec![2, 3, 4, 5, 6, 9]
            }
        )
    }
    #[test]
    fn test_get_motifs() {
        let mut f = std::fs::File::open("resources/netmhcout.txt").unwrap();
        //        let mut netmhcpan_summary = read_netmhcpan(&netmhcpan).unwrap();
        let mut netmhcpan_summary = read_netmhcpan(f).unwrap();
        let hla1 = "HLA-C03:04".parse::<HLA>().unwrap();
        let hla2 = "HLA-C08:01".parse::<HLA>().unwrap();

        let cd8 = vec![2, 3, 4, 5, 6, 9];
        let hla1_bound = netmhcpan_summary.get_bound(&hla1, Some(10_f32));
        let hla2_bound = netmhcpan_summary.get_bound(&hla2, Some(10_f32));
        let motifs1 = netmhcpan_summary.get_motifs(&hla1_bound, &cd8);
        let motifs2 = netmhcpan_summary.get_motifs(&hla2_bound, &cd8);

        assert_eq!(
            intersection_count_sorted_motifs(&motifs1, &motifs2) as f32 / motifs1.len() as f32,
            1f32
        );
    }
}
