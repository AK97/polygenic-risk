use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

const DEFAULT_GWAS: &str = "data/data.tsv";
const DEFAULT_GENOME: &str = "data/human1_hg38.bed";
const DEFAULT_P_THRESHOLD: f64 = 1e-2;

#[derive(Debug)]
struct GwasRecord {
    position: u64,
    effect_allele: String,
    beta: f64,
}

fn main() -> Result<(), Box<dyn Error>> {
    run_prs()?;
    Ok(())
}

fn run_prs() -> Result<(), Box<dyn Error>> {
    let gwas_path = PathBuf::from(DEFAULT_GWAS);
    let genome_path = PathBuf::from(DEFAULT_GENOME);
    let p_threshold = DEFAULT_P_THRESHOLD;

    let genome = read_genome(&genome_path)?;
    let gwas = read_gwas(&gwas_path, p_threshold)?;

    println!("Total # Signifiicant SNPs: {}", gwas.len());

    let (prs, overlap) = calculate_prs(&genome, &gwas);
    println!("Overlapping SNPs between genome and GWAS: {overlap}");
    println!("Polygenic Risk Score: {prs}");

    Ok(())
}

fn read_gwas(path: &Path, p_threshold: f64) -> io::Result<Vec<GwasRecord>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    let header_line = match lines.next() {
        Some(line) => line?,
        None => {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "GWAS file is empty; expected header row",
            ))
        }
    };

    let headers: Vec<&str> = header_line.trim_end().split('\t').collect();
    let idx_position = find_column(&headers, "base_pair_location")?;
    let idx_effect = find_column(&headers, "effect_allele")?;
    let idx_beta = find_column(&headers, "beta")?;
    let idx_pvalue = find_column(&headers, "p_value")?;

    let mut records = Vec::new();

    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= idx_pvalue
            || cols.len() <= idx_effect
            || cols.len() <= idx_beta
            || cols.len() <= idx_position
        {
            continue;
        }

        let position = match cols.get(idx_position).and_then(|s| s.trim().parse::<u64>().ok()) {
            Some(v) => v,
            None => continue,
        };
        let p_value = match cols
            .get(idx_pvalue)
            .and_then(|s| s.trim().parse::<f64>().ok())
        {
            Some(v) => v,
            None => continue,
        };
        if p_value >= p_threshold {
            continue;
        }

        let beta = match cols.get(idx_beta).and_then(|s| s.trim().parse::<f64>().ok()) {
            Some(v) => v,
            None => continue,
        };
        let effect_allele = match cols.get(idx_effect) {
            Some(v) => v.trim().to_ascii_uppercase(),
            None => continue,
        };
        if effect_allele.is_empty() {
            continue;
        }

        records.push(GwasRecord {
            position,
            effect_allele,
            beta,
        });
    }

    Ok(records)
}

fn find_column(headers: &[&str], name: &str) -> io::Result<usize> {
    headers
        .iter()
        .position(|h| h.eq_ignore_ascii_case(name))
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Missing column '{name}' in GWAS file header"),
            )
        })
}

fn read_genome(path: &Path) -> io::Result<HashMap<u64, String>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut genome = HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 4 {
            continue;
        }

        let start = match cols.get(1).and_then(|s| s.trim().parse::<u64>().ok()) {
            Some(v) => v,
            None => continue,
        };
        let genotype = cols.get(3).map(|s| s.trim().to_ascii_uppercase());
        if let Some(gt) = genotype {
            genome.insert(start, gt);
        }
    }

    Ok(genome)
}

fn calculate_prs(genome: &HashMap<u64, String>, gwas: &[GwasRecord]) -> (f64, usize) {
    let mut prs = 0.0;
    let gwas_positions: HashSet<u64> = gwas.iter().map(|r| r.position).collect();

    for record in gwas {
        if let Some(genotype) = genome.get(&record.position) {
            let count = count_effect_allele(genotype, &record.effect_allele);
            prs += (count as f64) * record.beta;
        }
    }

    let overlap = genome
        .keys()
        .filter(|pos| gwas_positions.contains(pos))
        .count();

    (prs, overlap)
}

fn count_effect_allele(genotype: &str, effect_allele: &str) -> u32 {
    let effect_char = match effect_allele.chars().next() {
        Some(c) => c,
        None => return 0,
    };

    if effect_allele.chars().count() != 1 {
        return 0;
    }

    genotype.chars().filter(|c| *c == effect_char).count() as u32
}
