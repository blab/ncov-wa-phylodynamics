'''
This script subsamples 25 strains from each interval (default is 14 days) since the first human sequence of ncov and writes them to a fasta file.
Inputs are:
    --alignment
    --metadata
    --interval (optional)
    --output

'''
import argparse
from Bio import SeqIO
import pandas as pd
import datetime
import random

def create_strain_list(alignment):
    '''
    Returns list of strains in alignment.
    '''
    strains = []
    for record in SeqIO.parse(alignment, 'fasta'):
        strains.append(record.id)
    return strains

def load_metadata(metadata, strains):
    '''
    Load metadata into dataframe.
    '''
    with open(metadata) as meta:
        df_meta = pd.read_csv(meta, sep='\t')
    df_meta = df_meta[df_meta['strain'].isin(strains)]
    df_meta['date'] = pd.to_datetime(df_meta['date'], infer_datetime_format=True)
    return df_meta

def select_strains(df_meta, interval):
    '''
    Selects 25 strains from every two-week interval since the date of the first virus sampled.
    '''
    min_date = df_meta.date.min()
    max_date = df_meta.date.max()
    days = datetime.timedelta(days=interval)
    sample = []
    for interval_start in (min_date + days*n for n in range(1 + int((max_date - min_date).days / interval))):
        sequences = df_meta['strain'][(df_meta['date'] < (interval_start + days)) & (df_meta['date'] >= interval_start)].tolist()
        if len(sequences) > 25:
            sample += random.sample(sequences, 25)
        else:
            sample += sequences
        interval_start += days
    return sample

def write_fasta(alignment, sample, output):
    '''
    Returns list of strains in alignment.
    '''
    sequences = SeqIO.parse(alignment, 'fasta')
    subsampled = [record for record in sequences if record.id in sample]
    SeqIO.write(subsampled, output, "fasta")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Subsamples strains over time.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type=str, required=True, help='alignment of viruses')
    parser.add_argument('--metadata', type=str, required=True, help = 'metadata for aligned sequences')
    parser.add_argument('--interval', type=int, default=14, help = 'interval  over which to sample viruses')
    parser.add_argument('--output', type=str, required=True, help='location of output file')
    args = parser.parse_args()

#Creates list of strains in alignment
strains = create_strain_list(args.alignment)

#Loads metadata to dataframe and filters to include only metadata for strains in alignment
df_meta = load_metadata(args.metadata, strains)

#Subsamples 25 strains for each interval since first virus sampled
sample = select_strains(df_meta, args.interval)

#writes subsampled strains to fasta file
write_fasta(args.alignment, sample, args.output)
