from Bio import SeqIO
import pandas as pd
import numpy as np
from io import StringIO 
import pyranges as pr
from numpy import array
import pybedtools
import time
from Bio.SeqUtils import gc_fraction
from itertools import product
import os

####################################################################################################################


#count di neucleotide count
def Counting_di(seq):
    sequence = seq
    di_nucleotides = [''.join(x) for x in product('ACGT', repeat=2)]
    # Generate all possible di-nucleotides
    di_count = {di: 0 for di in di_nucleotides} # Initialize dictionary with all di-nucleotides
  
    for i in range(len(sequence) - 1):
        di = sequence[i:i+2].upper() # Extract current di-nucleotide
        if di in di_count:
            di_count[di] += 1 # Increment count if di-nucleotide is valid
    
    #print(di_count)

    return di_count



#count tri neucleotide count
def Counting_tri(seq):
    sequence = seq
    tri_nucleotides = [''.join(x) for x in product('ACGT', repeat=3)] # Generate all possible tri-nucleotides
    tri_count = {tri: 0 for tri in tri_nucleotides} # Initialize dictionary with all tri-nucleotides
    
    for i in range(len(sequence) - 2):
        tri = sequence[i:i+3].upper() # Extract current tri-nucleotide
        if tri in tri_count:
            tri_count[tri] += 1 # Increment count if tri-nucleotide is valid
        

    return tri_count


#count tetra neucleotide count
def Counting_tetra(seq):
    sequence = seq
    tetra_nucleotides = [''.join(x) for x in product('ACGT', repeat=4)] # Generate all possible terai-nucleotides
    tetra_count = {tetra: 0 for tetra in tetra_nucleotides} # Initialize dictionary with all tri-nucleotides
    
    for i in range(len(sequence) - 3):
        tetra = sequence[i:i+4].upper() # Extract current tri-nucleotide
        if tetra in tetra_count :
            tetra_count[tetra] += 1 # Increment count if tri-nucleotide is valid
        

    return tetra_count




####################################################################################################################


start= time.time()

round_dir='/shares/kslotkin_share/private/tratnayake/rubyRF/Round2/'
fasta_file_list='/shares/kslotkin_share/private/tratnayake/rubyRF/Round2/fa_files.txt'


with open(fasta_file_list, 'r') as file:
    for line in file:
        
        fa_file=line.strip()
        infa_pth=os.path.join(round_dir ,fa_file)
        print(infa_pth)
        psplit=os.path.split(fa_file)[1]
        
        outdr=os.path.split(fa_file)[0]
        print(psplit)
        outfname=psplit.replace('_3p5ptrim_ruby.fasta', '')

        output_file= str(outfname)+'R2_3p5ptr_features_out.txt'
        output_pth=os.path.join(outdr,output_file)
        print(output_pth)




        first_iteration = True

        for seq_record in SeqIO.parse(infa_pth, 'fasta'):
        
            readid=seq_record.id
            print(readid)
            #print(repr(seq_record.seq))
            sequence = str(seq_record.seq)
            #print(len(sequence))
        
        
            gc_content=gc_fraction(sequence)
            #print('GC content ' + seq_record.id + ' ' + str(gc_content) + '%')
        
        
            
            read_dic= {'seq_id':seq_record.id,'length': len(sequence) , 'gc_cont': gc_content }
            #print(read_dic)
        
            di_dic=Counting_di(sequence)
            #print(di_dic)
            read_dic.update(di_dic)
        
            tri_dic=Counting_tri(sequence)
            #print(tri_dic)
            read_dic.update(tri_dic)

            tetra_dic=Counting_tetra(sequence)
            #print(tetra_dic)
            read_dic.update(tetra_dic)

            #print('final >>>>:' )

            #print(read_dic)

            ddf=pd.DataFrame.from_dict(read_dic,orient='index',columns=['Count'])

            dfn = (ddf.T).reset_index(drop=True)


        

            dfn.to_csv(output_pth, mode='a', header=first_iteration, index=False,sep='\t')
            #dfn.to_csv(output_file, mode='a', header=first_iteration, index=False)

            if first_iteration:
                first_iteration = False
    


etime=(time.time()-start)
print('Execution time : ' + str(etime) + 'Seconds')