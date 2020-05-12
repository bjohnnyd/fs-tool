use serde::de::{self, Visitor};
use serde::{ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer};

use crate::meta::{AlleleMeta, LigandMeta};
use immunoprot::ig_like::kir_ligand::{KirLigandInfo, LigandMotif};

struct LigandMetaVisitor;

impl Serialize for LigandMeta {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut ligand_info = serializer.serialize_struct("LigandMeta", 3)?;
        ligand_info.serialize_field("kir_ligand_allele", &self.kir_ligand_allele.to_string())?;
        ligand_info.serialize_field("kir_ligand_motif", &self.kir_ligand_motif.to_string())?;
        ligand_info.serialize_field("kir_ligand_freq", &self.kir_ligand_allele_freq.to_string())?;
        ligand_info.end()
    }
}

impl<'de> Visitor<'de> for LigandMetaVisitor {
    type Value = LigandMeta;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("Kir Ligand Information")
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: serde::de::Error,
    {
        if let Ok(ligand_info) = value.parse::<KirLigandInfo>() {
            Ok(LigandMeta::new(&ligand_info))
        } else {
            Err(E::custom(format!(
                "Could not create ligand motif information from {}",
                value
            )))
        }
    }
}

impl<'de> Deserialize<'de> for LigandMeta {
    fn deserialize<D>(deserializer: D) -> Result<Self, <D as Deserializer<'de>>::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(LigandMetaVisitor)
    }
}

impl Serialize for AlleleMeta {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut allele_meta = serializer.serialize_struct("LigandMeta", 6)?;
        allele_meta.serialize_field("allele", &self.allele.to_string())?;
        allele_meta.serialize_field("netmhcpan_nn", &self.netmhcpan_nn.to_string())?;
        allele_meta.serialize_field(
            "netmhcpan_nn_distance",
            &self.netmhcpan_nn_distance.to_string(),
        )?;

        if let Some(ligand_meta) = &self.ligand_meta {
            allele_meta.serialize_field(
                "kir_ligand_allele",
                &ligand_meta.kir_ligand_allele.to_string(),
            )?;
            allele_meta.serialize_field(
                "kir_ligand_motif",
                &ligand_meta.kir_ligand_motif.to_string(),
            )?;
            allele_meta.serialize_field(
                "kir_ligand_freq",
                &ligand_meta.kir_ligand_allele_freq.to_string(),
            )?;
        } else {
            allele_meta.serialize_field("kir_ligand_allele", "NA")?;
            allele_meta.serialize_field("kir_ligand_motif", "NA")?;
            allele_meta.serialize_field("kir_ligand_freq", "NA")?;
        }
        allele_meta.end()
    }
}

pub fn optional_motif_serialize<S>(x: &Option<LigandMotif>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match x {
        Some(motif) => s.serialize_str(motif.to_string().as_str()),
        None => s.serialize_str("NA"),
    }
}

pub fn optional_motif_deserialize<'de, D>(deserializer: D) -> Result<Option<LigandMotif>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;

    Ok(s.parse::<LigandMotif>().ok())
}

pub fn optional_bool_serialize<S>(x: Option<bool>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match x {
        Some(x) => s.serialize_str(if x { "1" } else { "0" }),
        None => s.serialize_str("NA"),
    }
}

pub fn optional_bool_deserialize<'de, D>(deserializer: D) -> Result<Option<bool>, D::Error>
where
    D: Deserializer<'de>,
{
    let s: &str = Deserialize::deserialize(deserializer)?;

    match s {
        "1" | "Y" | "T" | "TRUE" => Ok(Some(true)),
        "0" | "N" | "F" | "FALSE" => Ok(Some(false)),
        s => Err(de::Error::custom(format!(
            "Could not deduce TRUE/FALSE from {}",
            s
        ))),
    }
}

pub fn optional_float_serialize<S>(x: &Option<f32>, s: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    match x {
        Some(x) => s.serialize_str(format!("{:.2}", x).as_str()),
        None => s.serialize_str("NA"),
    }
}
