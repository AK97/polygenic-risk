Polygenic Risk Scoring
======================

This repo contains a small, end-to-end script to compute a polygenic risk score (PRS) from GWAS summary statistics and a single individual's genotypes.

The data used is Genome-Wide Association Data from the UK BioBank and individual genomic data from the Personal Genome Project. UCSC's LiftOver web tool was used to convert between genomic builds.

Files
------------------
- `main.py` loads GWAS summary stats (`data/data.tsv`) and a sample genome in BED format (`data/human1_hg38.bed`), filters variants by p-value, and prints the PRS.
- `import_data.py` provides helpers that load `.tsv`/`.bed` files from `data/` and cache them as `.pk1` for faster re-runs.
- `prs.py` aligns GWAS and genome variants by position and sums `copies_of_effect_allele * beta` to produce the PRS.
- `convert.py` is a helper for turning raw 23andMe-style TSV data into BED for LiftOver, if you need to regenerate the sample genome file.

Running it
----------
1) Install Python 3.10+ and dependencies:
```
python -m venv venv
source venv/bin/activate
pip install pandas
```
2) Place inputs in `data/`:
- GWAS summary stats: `data/<name>.tsv` with at least `base_pair_location`, `effect_allele`, `beta`, `p_value` columns. The default file used is `data/data.tsv`.
- Genome data: `data/<name>.bed` with columns `chromosome`, `start`, `end`, `genotype`, `dummy`. The default file used is `data/human1_hg38.bed`.
3) Run the script:
```
python main.py
```
It will report how many SNPs pass the p-value threshold (default `1e-2`) and the resulting PRS.

Results
-----
`> python3 main.py`
```
Total # Signifiicant SNPs: 1021228
Overlapping SNPs between genome and GWAS: 321
Polygenic Risk Score: 292.17575389999996
```

Notes
-----
- The loader caches to `data/*.pk1`; delete those if you need to force a fresh read.
- To work with a different dataset, change the filenames passed to `import_tsv` and `import_bed` in `main.py`.
- `convert.py` expects a TSV with columns `chromosome`, `position`, and `genotype`, filters to chromosome 10, converts to 0-based BED, and writes `data/human1.bed` as an intermediate for LiftOver.
