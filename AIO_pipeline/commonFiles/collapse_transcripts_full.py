#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import sys
import glob
import os

# Define paths from command line arguments
base_path = sys.argv[1]  # Path to the directory containing the files
cutoff = int(sys.argv[2])  # Cutoff value for proximity

# Read and merge Wt Always Red Files
wt_red_files = glob.glob(os.path.join(base_path, "*35S.1red.rRNA_free.mapped_to_targets.htseq.txt")) + \
               glob.glob(os.path.join(base_path, "*35S.2reddish.rRNA_free.mapped_to_targets.htseq.txt"))

# Initialize an empty DataFrame
wt_red = pd.DataFrame()
for i, file in enumerate(wt_red_files):
    df = pd.read_table(file, header=None, names=["Gene"] + [f"wt_red.{i+1}"])
    wt_red = pd.merge(wt_red, df, on="Gene", how="outer") if not wt_red.empty else df

# Fill NaNs with 0
wt_red = wt_red.fillna(0)

# Read and merge dcl1234 Always Red Files
dcl_red_files = glob.glob(os.path.join(base_path, "*35S_dcl1234.1red.rRNA_free.mapped_to_targets.htseq.txt"))

# Initialize an empty DataFrame
dcl_red = pd.DataFrame()
for i, file in enumerate(dcl_red_files):
    df = pd.read_table(file, header=None, names=["Gene"] + [f"dcl_red.{i+1}"])
    dcl_red = pd.merge(dcl_red, df, on="Gene", how="outer") if not dcl_red.empty else df

# Fill NaNs with 0
dcl_red = dcl_red.fillna(0)

# Read and merge Full Red Files
wt_fullred_files = glob.glob(os.path.join(base_path, "*35S.5fullRed.rRNA_free.mapped_to_targets.htseq.txt"))

# Initialize an empty DataFrame
wt_fullred = pd.DataFrame()
for i, file in enumerate(wt_fullred_files):
    df = pd.read_table(file, header=None, names=["Gene"] + [f"wt_fullRed.{i+1}"])
    wt_fullred = pd.merge(wt_fullred, df, on="Gene", how="outer") if not wt_fullred.empty else df

# Fill NaNs with 0
wt_fullred = wt_fullred.fillna(0)

# Combine dataframes for all phenotypes
count_data = wt_red.merge(dcl_red, on="Gene", how="inner").merge(wt_fullred, on="Gene", how="inner")

# Filter transcripts present in at least 2 replicates of at least 1 sample
def filter_transcripts(df):
    return df[
        (df.filter(like="wt_red.").gt(0).sum(axis=1) >= 2) |
        (df.filter(like="dcl_red.").gt(0).sum(axis=1) >= 2) |
        (df.filter(like="wt_fullRed.").gt(0).sum(axis=1) >= 2)
    ]

count_data_filtered = filter_transcripts(count_data)

# Ensure that we are working with a copy of the DataFrame
count_data_filtered = count_data_filtered.copy()

# Convert 'Gene' column to string (if it is not already)
count_data_filtered['Gene'] = count_data_filtered['Gene'].astype(str)

# Split 'Gene' column into 'chr', 'coords', and 'class'
split_gene = count_data_filtered['Gene'].str.split('.', expand=True)

# Ensure split_gene has the correct number of columns
expected_columns = ['chr', 'coords', 'class']
if len(split_gene.columns) == len(expected_columns):
    count_data_filtered[expected_columns] = split_gene
else:
    print(f"Unexpected number of columns in split_gene: {len(split_gene.columns)}")

# Split 'coords' into 'start' and 'end'
split_coords = count_data_filtered['coords'].str.split('_', expand=True)
count_data_filtered[['start', 'end']] = split_coords

# Split 'class' into 'type' and 'strand'
split_class = count_data_filtered['class'].str.split('_', expand=True)
count_data_filtered[['type', 'strand']] = split_class

# Drop the now redundant 'coords' and 'class' columns
count_data_filtered = count_data_filtered.drop(columns=['coords', 'class', 'Gene'])

# Reorder columns: Place 'chr', 'start', 'end', 'type', 'strand' first, then others
desired_order = ['chr', 'start', 'end', 'type', 'strand'] + [col for col in count_data_filtered.columns if col not in ['chr', 'start', 'end', 'type', 'strand']]
count_data_filtered = count_data_filtered[desired_order]

# Sort DataFrame by 'chr', 'end', 'start'
count_data_sorted = count_data_filtered.sort_values(by=['chr', 'end', 'start'], ascending=[True, True, True]).reset_index(drop=True)

def get_collapse_transcripts(df, cutoff):
    # Calculate next and previous positions
    df['next_start'] = df.groupby(['chr', 'strand', 'type'])['start'].shift(-1)
    df['next_end'] = df.groupby(['chr', 'strand', 'type'])['end'].shift(-1)
    df['prev_start'] = df.groupby(['chr', 'strand', 'type'])['start'].shift(1)
    df['prev_end'] = df.groupby(['chr', 'strand', 'type'])['end'].shift(1)

    # Ensure numeric columns are of type int or float
    df[['next_start', 'next_end', 'prev_start', 'prev_end', 'start', 'end']] = \
        df[['next_start', 'next_end', 'prev_start', 'prev_end', 'start', 'end']].apply(pd.to_numeric, errors='coerce')

    # Determine proximity to next and previous rows
    df['nextLine'] = ((abs(df['next_start'] - df['start']) < cutoff) & 
                      (abs(df['next_end'] - df['end']) < cutoff)).astype(int)
    df['prevLine'] = ((abs(df['start'] - df['prev_start']) < cutoff) & 
                      (abs(df['end'] - df['prev_end']) < cutoff)).astype(int)

    # Calculate the logical sum
    df['logicalSum'] = df['nextLine'] + df['prevLine']

    # Drop intermediate columns
    df = df.drop(columns=['next_start', 'next_end', 'prev_start', 'prev_end'])

    # Initialize grouping variables
    groups = [0] * len(df)
    current_group = 0
    current_chr = df.loc[0, 'chr']
    
    # Iterate through each row to assign group numbers
    for i in range(len(df)):
        # Reset group number for a new chromosome
        if df.loc[i, 'chr'] != current_chr:
            current_group = 0
            current_chr = df.loc[i, 'chr']
        
        # Assign group number based on previous row's logical sum and strand
        if i == 0 or df.at[i-1, 'chr'] != df.at[i, 'chr']:
            groups[i] = current_group
        else:
            if df.at[i-1, 'logicalSum'] == 0:
                current_group += 1
                groups[i] = current_group
            elif df.at[i-1, 'logicalSum'] > 0:
                if df.at[i, 'prevLine'] == 1 and df.at[i, 'strand'] == df.at[i-1, 'strand']:
                    groups[i] = current_group
                else:
                    current_group += 1
                    groups[i] = current_group

    # Assign the group values to the DataFrame
    df['group'] = groups
    df = df.drop(columns=['nextLine', 'prevLine', 'logicalSum'])

    # Specify columns to exclude from TotalSum calculation
    exclude_columns = ['chr', 'start', 'end', 'strand', 'type', 'group']
    value_columns = [col for col in df.columns if col not in exclude_columns]

    # Calculate TotalSum column
    df['TotalSum'] = df[value_columns].sum(axis=1)

    # Define aggregation functions
    agg_funcs = {
        'start': lambda x: x[df.loc[x.index, 'TotalSum'].idxmax()],
        'end': lambda x: x[df.loc[x.index, 'TotalSum'].idxmax()],
    }

    # Add sum aggregations for value columns
    agg_funcs.update({col: 'sum' for col in value_columns})

    # Aggregate using groupby
    result = df.groupby(['chr', 'strand', 'type', 'group']).agg(agg_funcs).reset_index()

    # Reorder columns
    result = result[['chr', 'start', 'end', 'type', 'group', 'strand'] + [col for col in result.columns if col not in ['chr', 'start', 'end', 'group', 'type', 'strand']]]

    return result

# Collapse transcripts with initial cutoff
results1 = get_collapse_transcripts(count_data_sorted, cutoff)

# Sort results from round 1
results1_sorted = results1.sort_values(by=['chr', 'start', 'end'], ascending=[True, True, True]).reset_index(drop=True)

# Collapse transcripts round 2

finalCollapse =get_collapse_transcripts(results1_sorted,cutoff)

# Step 1: Combine 'start' and 'end' into 'coords'
finalCollapse['coords'] = finalCollapse['start'].astype(str) + '_' + finalCollapse['end'].astype(str)

# Step 2: Combine 'type' and 'strand' into 'type'
finalCollapse['type'] = finalCollapse['type'] + '_' + finalCollapse['strand']

# Step 3: Combine 'chr', 'coords', and 'type' into 'Gene'
finalCollapse['Gene'] = finalCollapse['chr'] + '.' + finalCollapse['coords'] + '.' + finalCollapse['type'] + '.transcript' + finalCollapse['group'].astype(str) 

# Drop intermediate columns if needed
finalCollapse = finalCollapse.drop(columns=['chr','start','coords','type', 'end','group', 'strand'])

reordered_columns = ['Gene'] + [col for col in finalCollapse.columns if col != 'Gene']
finalCollapse = finalCollapse[reordered_columns]

#print(finalCollapse.to_string(index=False, header=True))

finalCollapse.to_csv('merged_input_deseq.txt', index=False,sep="\t",header=True )


