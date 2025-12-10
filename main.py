from import_data import import_tsv, import_bed
from prs import calculate_prs

data = import_tsv("data")
genome_data = import_bed("human1_hg38")

# 102182724 count in gwas sample
# 1014604 count in human sample

# we care about risk allele and effect sizes
relevant_data = data[['base_pair_location', 'effect_allele', 'beta', 'p_value']]
relevant_genome_data = genome_data[['start', 'genotype']]

# only consider SNPs with statistically significant p-values
P_VAL_THRESHOLD = 1e-2
significant_snps = relevant_data[relevant_data['p_value'] < P_VAL_THRESHOLD][['base_pair_location', 'effect_allele', 'beta']]

print(f"Total # Signifiicant SNPs: {len(significant_snps)}")

prs = calculate_prs(relevant_genome_data, significant_snps)
print(f"Polygenic Risk Score: {prs}")
