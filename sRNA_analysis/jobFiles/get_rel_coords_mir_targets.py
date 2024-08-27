import pandas as pd
import re
import sys
import csv  # Import CSV module to control quoting
targetsFile = sys.argv[1]
inTargets = open(targetsFile,'r') #"mirbase_athaliana.uniq.names.no_ppt.targets.ath_deg.txt"

targets = []
genes = []

for i in inTargets:
    i = i.rstrip("\n")
    targets.append(i)
    x = i.split(".")[0]
    genes.append(x)

# Load the GTF file
fileName = sys.argv[2] #'Araport11_GTF_genes_transposons.current.gtf'
inGTF = pd.read_csv(fileName, sep='\t', header=None)
column_names = ['chr', 'araport', 'type', 'start', 'end', 'foo', 'strand', 'foo2', 'info']
inGTF.columns = column_names

# Functions to extract transcript_id and gene_id
def extract_transcript_id(text):
    match = re.search(r'transcript_id "([^"]+)"', text)
    return match.group(1) if match else None

def extract_gene_id(text):
    match = re.search(r'gene_id "([^"]+)"', text)
    return match.group(1) if match else None

# Filter rows based on conditions
def filter_rows(row):
    transcript_id = extract_transcript_id(row['info'])
    gene_id = extract_gene_id(row['info'])
    if transcript_id in targets:
        return True
    return False

filtered_df = inGTF[inGTF.apply(filter_rows, axis=1)]
filtered_df = filtered_df[filtered_df['type'] != 'protein']

# Function to process each group
def adjust_positions(group):
    transcript_id = extract_transcript_id(group.iloc[0]['info'])
    gene_id = extract_gene_id(group.iloc[0]['info'])

    if group['strand'].iloc[0] == '+':
        # Get the mRNA start position
        mrna_start = group.loc[group['type'] == 'mRNA', 'start'].values[0]
        
        # Subtract the mRNA start from start and stop columns
        group['relStart'] = (group['start'] - mrna_start) +1
        group['relEnd'] = (group['end'] - mrna_start) +1
        
    elif group['strand'].iloc[0] == '-':
        # Get the mRNA end position
        mrna_end = group.loc[group['type'] == 'mRNA', 'end'].values[0]
        
        # Adjust positions for the '-' strand
        group['relEnd'] = ((group['start'] - mrna_end) * -1)+1
        group['relStart'] = ((group['end'] - mrna_end) * -1)+1
        group['strand'] = "+"

    group['chr'] = transcript_id
    group['info'] = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
    return group

# Apply adjustments
df_adjusted = filtered_df.groupby('info').apply(adjust_positions)
df_adjusted.reset_index(drop=True, inplace=True)
df_adjusted = df_adjusted.drop(columns=['start', 'end'])

# Reorder columns
df_adjusted = df_adjusted[['chr', 'araport', 'type', 'relStart', 'relEnd', 'foo', 'strand', 'foo2', 'info']]

# Sort and save the output
df_adjusted = df_adjusted.sort_values(by=['info', 'relStart', 'type'], ascending=[True, True, False]).reset_index(drop=True)
outFile = targetsFile.split(".txt")[0]+"."+fileName.replace(".gtf", "_relCoords.gtf").split("/")[-1]
#outFile = fileName.replace(".gtf", "_relCoords.gtf")


df_adjusted.to_csv(outFile, index=False, sep="\t", header=False, quoting=csv.QUOTE_NONE, escapechar='\\')
