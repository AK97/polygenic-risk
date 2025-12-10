import pandas as pd

def calculate_prs(genome:pd.DataFrame, gwas:pd.DataFrame) -> float:
    """
    Calculate the polygenic risk score (PRS) for a single genome.

    Parameters
    ----------
    genome : DataFrame
        Columns:
          - 'start'    : genomic position (same coordinate system as GWAS)
          - 'genotype' : two-letter string, e.g. 'AA', 'GC'
    gwas : DataFrame
        Columns:
          - 'base_pair_location' : genomic position matching 'start'
          - 'effect_allele'      : allele for which beta is reported
          - 'beta'               : effect size per copy of effect_allele

    Returns
    -------
    float
        The PRS = sum_over_SNPs( copies_of_effect_allele * beta ).
    """

    # Align variants by genomic position
    merged = genome.merge(
        gwas[['base_pair_location', 'effect_allele', 'beta']],
        left_on='start',
        right_on='base_pair_location',
        how='inner'
    )

    # Normalize to uppercase strings
    genotypes = merged['genotype'].astype(str).str.upper()
    effect_alleles = merged['effect_allele'].astype(str).str.upper()

    # Count how many copies (0, 1, 2) of the effect_allele are in the genotype string
    # e.g. genotype 'GC', effect_allele 'C' -> 1
    #      genotype 'AA', effect_allele 'A' -> 2
    effect_counts = [
        g.count(ea) if len(g) == 2 else None
        for g, ea in zip(genotypes, effect_alleles)
    ]

    merged['effect_allele_count'] = pd.to_numeric(effect_counts, errors='coerce')

    # Drop rows with missing counts or beta
    merged = merged.dropna(subset=['effect_allele_count', 'beta'])

    # PRS = sum( copies_of_effect_allele * beta )
    prs = (merged['effect_allele_count'] * merged['beta']).sum()

    # Count overlap
    overlap = genome['start'].isin(gwas['base_pair_location']).sum()
    print(f"Overlapping SNPs between genome and GWAS: {overlap}")

    return float(prs)
