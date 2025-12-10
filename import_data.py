import csv
import pandas as pd

def import_tsv(filename:str) -> pd.DataFrame:
    # pickling for faster subsequent loads
    try:
        data = pd.read_pickle(f"data/{filename}.pk1")
    except:
        data = pd.read_csv(f"data/{filename}.tsv", delimiter = '\t')
        data.to_pickle(f"data/{filename}.pk1")

    return data

def import_bed(filename:str) -> pd.DataFrame:
    # pickling for faster subsequent loads
    try:
        data = pd.read_pickle(f"data/{filename}.pk1")
    except:
        data = pd.read_csv(f"data/{filename}.bed", delimiter = '\t', names=['chromosome', 'start', 'end', 'genotype', 'dummy'])
        data.to_pickle(f"data/{filename}.pk1")

    return data