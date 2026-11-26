use anyhow::{Context, Result};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

#[derive(Debug, Clone, Deserialize)]
pub struct PrioritizeConfig {
    pub af: Vec<String>,
    pub maxaf: f32,
    pub regions: Vec<String>,
    pub scores: BTreeMap<String, f32>,
    pub constraint: f32,
    #[serde(default)]
    pub more_regions: Vec<String>,
    #[serde(default)]
    pub more_values: BTreeMap<String, f32>,
}

impl PrioritizeConfig {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let p = path.as_ref();
        let file =
            File::open(p).with_context(|| format!("Failed to open config file: {}", p.display()))?;
        let reader = BufReader::new(file);
        serde_json::from_reader(reader)
            .with_context(|| format!("Failed to parse config JSON: {}", p.display()))
    }
}

#[derive(Debug, Clone, Copy, Deserialize)]
pub struct DbSchema {
    pub id: usize,
    pub stdtype: usize,
    pub dbsource: usize,
    pub constraint: usize,
    pub closest: usize,
    #[serde(rename = "closestProt")]
    pub closest_prot: usize,
    pub genes: usize,
}

impl DbSchema {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let p = path.as_ref();
        let file = File::open(p)
            .with_context(|| format!("Failed to open dbschema file: {}", p.display()))?;
        let reader = BufReader::new(file);
        serde_json::from_reader(reader)
            .with_context(|| format!("Failed to parse dbschema JSON: {}", p.display()))
    }

    pub fn max_col(&self) -> usize {
        self.id
            .max(self.stdtype)
            .max(self.dbsource)
            .max(self.constraint)
            .max(self.closest)
            .max(self.closest_prot)
            .max(self.genes)
    }
}
