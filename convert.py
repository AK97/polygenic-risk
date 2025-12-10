# Convert from TSV to BED format for LiftOver
# For crossmapping genome variant with GWAS data

import pandas as pd

def tsv_to_bed(filename:str) -> None:

    df = pd.read_csv(filename, delimiter = "\t")

    print(df.head())

    # GWAS only has chromosome 10 data; filter
    df = df[df['chromosome'] == 10]

    # BED uses 0-based start, 1-based end
    df["chrom"] = "chr" + df["chromosome"].astype(str)
    df["start"] = df["position"] - 1
    df["end"]   = df["position"]

    # Columns for BED: chrom, start, end, genotype
    bed = df[["chrom", "start", "end", "genotype"]]

    print(bed.head())

    # Write BED
    bed.to_csv("data/human1.bed", sep = "\t", header = False, index = False)

tsv_to_bed("data/human1_23andme.tsv")