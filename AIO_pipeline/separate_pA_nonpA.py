#pip install cutadapt
#pip install pysam

import pysam
import subprocess
import os
import sys
import cutadapt
inAnnotation=open("/cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.tar
getChr.construct.transgenes.bed","r")

senseTargets = []
antiTargets=[]
for line in inAnnotation:
    data = line.rstrip().split('\t')
    if data[5] == "+" :
        senseTargets.append(data[0])
    elif data[5] == "-":
        antiTargets.append(data[0])



def extract_soft_clipped_sequences(bam_file, output_fasta_for,output_fasta_rev):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    with open(output_fasta_for, 'w') as fasta_file_for, open(output_fasta_rev, 'w') as fasta_file_rev:
        # Iterate over each read in the BAM file
        for read in bam:
            # Check if the read is mapped (flag not equal to 4)
            if not read.is_unmapped:
                # Check for soft clipping in the CIGAR string
                if read.cigartuples:
                    if read.is_forward:
                        # Find if there is a soft clip at the 3' end
                        if read.cigartuples[-1][0] == 4:  # Soft clipping operation at the end of the read
                            length = read.cigartuples[-1][1]
                            if read.query_sequence is not None:
                                sequence = read.query_sequence[-length:]
                                fasta_file_for.write(f">{read.query_name}\n{sequence}\n")
                    if read.is_reverse:
                        # Find if there is a soft clip at the 5' end
                        if read.cigartuples[0][0] == 4:  # Soft clipping operation at the beginning of the read
                            length = read.cigartuples[0][1]
                            if read.query_sequence is not None:
                                sequence = read.query_sequence[:length]
                                fasta_file_rev.write(f">{read.query_name}\n{sequence}\n")
    
    bam.close()

def check_soft_clips(bam_file, output_fasta_for,output_fasta_rev):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    with open(output_fasta_for, 'w') as fasta_file_for, open(output_fasta_rev, 'w') as fasta_file_rev:
        # Iterate over each read in the BAM file
        for read in bam:
            # Check if the read is mapped (flag not equal to 4)
            if not read.is_unmapped:
                # Check for soft clipping in the CIGAR string
                if read.cigartuples:
                    reference_name = bam.get_reference_name(read.reference_id)
                    # Get Reference name
                    if reference_name in senseTargets:
                        if read.is_reverse:
                            # Find if there is a soft clip at the 3' end
                            if read.cigartuples[-1][0] == 4:  # Soft clipping operation at the end of the read
                                length = read.cigartuples[-1][1]
                                if read.query_sequence is not None:
                                    sequence = read.query_sequence[-length:]
                                    fasta_file_for.write(f">{read.query_name}\n{sequence}\n")
                    else:
                        if read.is_forward:
                            # Find if there is a soft clip at the 3' end
                            if read.cigartuples[0][0] == 4:  # Soft clipping operation at the beginning of the read
                                length = read.cigartuples[0][1]
                                if read.query_sequence is not None:
                                    sequence = read.query_sequence[:length]
                                    fasta_file_rev.write(f">{read.query_name}\n{sequence}\n")
    
    bam.close()


def search_pA_for_with_cutadapt(input_fasta, search_output_fasta,untrim_out):
    # Run cutadapt to search for the specific sequence
    command = [
        sys.executable, "-m", "cutadapt",
        "--quiet",
        "-g", "XA{50};max_error_rate=0.10;min_overlap=4;",
        "-g", "A{20};max_error_rate=0.25;min_overlap=8;",
        "--action=lowercase",
        "-o", search_output_fasta,
        "--untrimmed-output",untrim_out,
        input_fasta
    ]
    subprocess.run(command, check=True)
    
def search_pA_rev_with_cutadapt(input_fasta, search_output_fasta,untrim_out):
    # Run cutadapt to search for the specific sequence
    command = [
        sys.executable, "-m", "cutadapt",
        "--quiet",
        "-a", "A{50}X;max_error_rate=0.10;min_overlap=4;",
        "-a", "A{20};max_error_rate=0.25;min_overlap=8;",
        "--action=lowercase",
        "-o", search_output_fasta,
        "--untrimmed-output",untrim_out,
        input_fasta
    ]
    subprocess.run(command, check=True)
    
def search_pT_for_with_cutadapt(input_fasta,search_output_fasta2,untrim_out):
    # Run cutadapt to search for the specific sequence
    command = [
        sys.executable, "-m", "cutadapt",
        "--quiet",
        "-g", "T{50};max_error_rate=0.10;min_overlap=10;",
        "--action=lowercase",
        "-o", search_output_fasta2,
        "--untrimmed-output",untrim_out,
        input_fasta
    ]
    subprocess.run(command, check=True)

def search_pT_rev_with_cutadapt(input_fasta,search_output_fasta2,untrim_out):
    # Run cutadapt to search for the specific sequence
    command = [
        sys.executable, "-m", "cutadapt",
        "--quiet",
        "-a", "T{50};max_error_rate=0.10;min_overlap=10;",
        "--action=lowercase",
        "-o", search_output_fasta2,
        "--untrimmed-output",untrim_out,
        input_fasta
    ]
    subprocess.run(command, check=True)

def search_3p_adapter_with_cutadapt(input_fasta,search_output_fasta3,untrim_out):
    # Run cutadapt to search for the specific sequence
    command = [
        sys.executable, "-m", "cutadapt",
        "--quiet",
        "-g", "CTGTAGGCACCATCAAT",
        "--action=lowercase",
        "-o", search_output_fasta3,
        "--untrimmed-output",untrim_out,
        input_fasta
    ]
    subprocess.run(command, check=True)


def process_bam_file_with_cutadapt(bam_file, output_bam_a, output_bam_b):
    extracted_fasta_for = "soft_clipped_sequences.for.tmp.fasta"
    extracted_fasta_rev = "soft_clipped_sequences.rev.tmp.fasta"
    search_pA_output_fasta_for = "cutadapt_output.for.pA.fasta"
    search_pA_output_fasta_rev = "cutadapt_output.rev.pA.fasta"
    search_pT_output_fasta_for = "cutadapt_output.for.dT.fasta"
    search_pT_output_fasta_rev = "cutadapt_output.rev.dT.fasta"
    untrim_out = "non_trim_output.fasta"

    # Extract 3' soft-clipped sequences to a FASTA file
    print(f"extracting soft-clipped sequences")
    extract_soft_clipped_sequences(bam_file, extracted_fasta_for, extracted_fasta_rev)

    # Search for the specific sequence in the extracted soft-clipped sequences using cutadapt
    print(f"looking for pA tails")
    search_pA_for_with_cutadapt(extracted_fasta_for, search_pA_output_fasta_for,untrim_out)
    search_pA_rev_with_cutadapt(extracted_fasta_rev, search_pA_output_fasta_rev,untrim_out)
    
    
    print(f"looking for pT tails")
    search_pT_for_with_cutadapt(extracted_fasta_for, search_pT_output_fasta_for,untrim_out)
    search_pT_rev_with_cutadapt(extracted_fasta_rev, search_pT_output_fasta_rev,untrim_out)
    
    
    # Read the names of the reads that contain the pA
    reads_with_pA = set()
    with open(search_pA_output_fasta_for, 'r') as pA_fasta_file_for, open(search_pA_output_fasta_rev, 'r') as pA_fasta_file_rev:
        for line in pA_fasta_file_for:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pA.add(read_name)
        for line in pA_fasta_file_rev:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pA.add(read_name)

    # Read the names of the reads that contain the pT
    reads_with_pT_for = set()
    reads_with_pT_rev = set()
    with open(search_pT_output_fasta_for, 'r') as pT_fasta_file_for,open(search_pT_output_fasta_rev, 'r') as pT_fasta_file_rev:
        for line in pT_fasta_file_for:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pT_for.add(read_name)
        for line in pT_fasta_file_rev:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pT_rev.add(read_name)

    # Open the BAM file for reading and create two output BAM files
    print(f"separating pA and non pA bam file")
    bam = pysam.AlignmentFile(bam_file, "rb")
    tmp_bam_a = "unsorted_" + output_bam_a
    tmp_bam_b = "unsorted_" + output_bam_b
    bam_a = pysam.AlignmentFile(tmp_bam_a, "wb", template=bam)
    bam_b = pysam.AlignmentFile(tmp_bam_b, "wb", template=bam)

    # Write the reads to the appropriate BAM files
    pAcount = 0
    nonpAcount = 0
    pTcount =0
    
    for read in bam:
        if read.query_name in reads_with_pA:
            bam_a.write(read)
            pAcount+=1
        elif read.query_name in reads_with_pT_rev:
            bam_a.write(read)
            pAcount+=1
        else:
            if read.query_name not in reads_with_pT_for:
                bam_b.write(read)
                nonpAcount+=1
    
    bam.close() 
    bam_a.close()
    bam_b.close()
    print(f"{pAcount} reads are pA+;\n{nonpAcount} reads are pA-;")

    # Sort and index the BAM files
    print(f"sorting")
    pysam.sort("-o", output_bam_a, tmp_bam_a)
    pysam.sort("-o", output_bam_b, tmp_bam_b)

    # Clean up temporary files
    os.remove(extracted_fasta_for)
    os.remove(extracted_fasta_rev)
    os.remove(search_pA_output_fasta_for)
    os.remove(search_pA_output_fasta_rev)
    os.remove(search_pT_output_fasta_for)
    os.remove(search_pT_output_fasta_rev)
    os.remove(untrim_out)
    os.remove(tmp_bam_a)
    os.remove(tmp_bam_b)

def check_non_pA_reads(pA_file, non_pA_file, output_bam_a, output_bam_b):
    extracted_fasta_for = "check_soft_clipped_sequences.for.tmp.fasta"
    extracted_fasta_rev = "check_soft_clipped_sequences.rev.tmp.fasta"
    search_pA_output_fasta_for = "check_cutadapt_output.for.pA.fasta"
    search_pT_output_fasta_rev = "check_cutadapt_output.rev.pT.fasta"
    search_3p_adapter_output_fasta_rev = "check_cutadapt_output.for.3p_adapter.fasta"
    untrim_out = "non_trim_output.fasta"

    # Extract 3' soft-clipped sequences to a FASTA file
    print(f"extracting soft-clipped sequences")
    check_soft_clips(non_pA_file, extracted_fasta_for, extracted_fasta_rev)

    # Search for the specific sequence in the extracted soft-clipped sequences using cutadapt
    print(f"looking for pA tails")
    search_pA_for_with_cutadapt(extracted_fasta_for, search_pA_output_fasta_for,untrim_out)

    
    print(f"looking for pT tails")
    search_pT_rev_with_cutadapt(extracted_fasta_rev, search_pT_output_fasta_rev,untrim_out)
    
    print(f"looking for 3' adapter")
    search_3p_adapter_with_cutadapt(extracted_fasta_for,search_3p_adapter_output_fasta_rev,untrim_out)
    
    # Read the names of the reads that contain the pA
    reads_with_pA = set()
    with open(search_pA_output_fasta_for, 'r') as pA_fasta_file_for, open(search_pT_output_fasta_rev, 'r') as pT_fasta_file_rev:
        for line in pA_fasta_file_for:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pA.add(read_name)
        for line in pT_fasta_file_rev:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_with_pA.add(read_name)

    reads_reorient = set()
    with open(search_3p_adapter_output_fasta_rev, 'r') as adapter_output:
        for line in adapter_output:
            if line.startswith('>'):
                read_name = line[1:].strip()
                reads_reorient.add(read_name)

    # Open the BAM file for reading and create two output BAM files
    print(f"separating pA and non pA bam file")
    bam = pysam.AlignmentFile(non_pA_file, "rb")
    tmp_bam_a = "unsorted_check_" + output_bam_a
    tmp_bam_b = "unsorted_check_" + output_bam_b
    bam_a = pysam.AlignmentFile(tmp_bam_a, "wb", template=bam)
    bam_b = pysam.AlignmentFile(tmp_bam_b, "wb", template=bam)

    # Write the reads to the appropriate BAM files
    pAcount = 0
    nonpAcount = 0
    
    for read in bam:
        if read.query_name in reads_with_pA:
            if read.is_reverse:
                read.is_reverse = False
            bam_a.write(read)
            pAcount+=1
        elif read.query_name in reads_reorient:
            if read.is_reverse:
                read.is_reverse = False
            bam_b.write(read)
            nonpAcount+=1
        else:
                bam_b.write(read)
                nonpAcount+=1
    
    bam.close() 
    bam_a.close()
    bam_b.close()
    print(f"{pAcount} reads are pA+;\n{nonpAcount} reads are pA-;")

    # Sort and index the BAM files
    print(f"merging pA bam files")
    tmp_merge_file = "merge_pA_file.tmp.bam"
    pysam.merge("-f",tmp_merge_file, pA_file, tmp_bam_a)

    print(f"sorting and indexing")
    pysam.sort("-o", output_bam_a, tmp_merge_file)
    pysam.sort("-o", output_bam_b, tmp_bam_b)
    pysam.index(output_bam_a)
    pysam.index(output_bam_b)

    # Clean up temporary files
    os.remove(extracted_fasta_for)
    os.remove(extracted_fasta_rev)
    os.remove(search_pA_output_fasta_for)
    os.remove(search_pT_output_fasta_rev)
    os.remove(search_3p_adapter_output_fasta_rev)
    os.remove(untrim_out)
    os.remove(tmp_bam_a)
    os.remove(tmp_bam_b)
    os.remove(tmp_merge_file)

# Usage example

bam_file = sys.argv[1] #"03_MK016_rep1.35S.1red.mod.non_pA.filtered.rRNA_free.bam"

tmp_output_bam_a = bam_file.split("bam")[0]+"pA.tmp.bam" #"output_with_sequence.bam"
tmp_output_bam_b = bam_file.split("bam")[0]+"non_pA.tmp.bam" #"output_without_sequence.bam"
print(f"Reading {bam_file}")
# Process the BAM file
process_bam_file_with_cutadapt(bam_file, tmp_output_bam_a, tmp_output_bam_b)


output_bam_a = bam_file.split("bam")[0]+"pA.bam" #"output_with_sequence.bam"
output_bam_b = bam_file.split("bam")[0]+"non_pA.bam" #"output_without_sequence.bam"
print(f"Checking {bam_file}")
check_non_pA_reads(tmp_output_bam_a, tmp_output_bam_b, output_bam_a, output_bam_b)

print(f"removing tmp files")
os.remove(tmp_output_bam_a)
os.remove(tmp_output_bam_b)

print(f"polyA+ reads saved to {output_bam_a}")
print(f"polyA- reads saved to {output_bam_b}")

######
## Example usage
# cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/2_map_to_capture/output/mapped_to_targets

# for file in *targets.bam; do python /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/separate_pA_nonpA.py $file ; done &

