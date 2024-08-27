#! /usr/bin/python2.7

import sys
import argparse
import subprocess
import re
import random

parser = argparse.ArgumentParser(description='Calculate basewise coverage over BED intervals using BIGWIG format as an intermediate')
parser.add_argument('in_ref', help='BED12 reference file')
parser.add_argument('out_bed', help='BED12 files with comma seperated list of scores in column 13')
parser.add_argument('--temp_directory_path', '-tdp', action='store', help='specify temporary directory path')
parser.add_argument('--bigWig_plus', '-plus', action='store', help="Plus strand bigwig file for stranded info")
parser.add_argument('--bigWig_minus', '-minus', action='store', help="Minus strand bigwig file for stranded info")
parser.add_argument('--bigWig_unstranded', '-unstr', action='store', help="bigWig input for unstranded data")
parser.add_argument('--NA_fill_value', '-NAval', action='store', help="Value to insert for positions with no value in the bigwig files. Must be numeric")
args = parser.parse_args()

# Set bigWig inputs
assert not ( args.bigWig_unstranded and (args.bigWig_plus or args.bigWig_minus) ) and ( (args.bigWig_plus and args.bigWig_minus) or (not args.bigWig_plus and not args.bigWi
g_minus) ), "must either input bigwigs for both strands using '-plus' and '-minus' or input a single bigwig using '-unstr'"

if args.bigWig_unstranded:
	bwPlus_file = args.bigWig_unstranded
	bwMinus_file = args.bigWig_unstranded

else:
	bwPlus_file = args.bigWig_plus
	bwMinus_file = args.bigWig_minus

# Set NA fill
if args.NA_fill_value:
	fill = ["-fill=" + args.NA_fill_value]
else:
	fill = []

# Generate random tag and temporary directory
def set_temp_dir(out_fn,tmp_tag):
	rTag = str(random.randint(1000000,9999999))
	l_slash = out_fn.rfind('/')
	tmpPth=out_fn[:l_slash+1] + tmp_tag + rTag + '/'
	subprocess.call(['mkdir','-p',tmpPth])
	return [rTag,tmpPth]

(rTag,tmpPth) = set_temp_dir(args.out_bed, 'tmp_bwscore_')

# reformat bed file reference
openRef = open(args.in_ref,'r')
tmpRef_file = tmpPth+"tmp_ref_"+rTag+".txt"
tmpRef_open = open(tmpRef_file,'w')

print >> sys.stderr, "reformating reference"
for line in openRef:
	j = line.rstrip().split("\t")
	tmpRef_open.write("\t".join([j[0],j[1],j[2],")@@(".join(j[3:len(j)])]) + "\n")

print >> sys.stderr, "done"

openRef.close()
tmpRef_open.close()

# transfer data to BED reference
bedPlus_file = tmpPth+"tmp_bedP_"+rTag+".txt"
bedMinus_file = tmpPth+"tmp_bedM_"+rTag+".txt"
print >> sys.stderr, "extracting coverage to BED format"
bedP = subprocess.Popen(["bwtool","extract","bed",tmpRef_file,bwPlus_file,bedPlus_file] + fill)
bedM = subprocess.Popen(["bwtool","extract","bed",tmpRef_file,bwMinus_file,bedMinus_file] + fill)
bedP.wait()
bedM.wait()
print >> sys.stderr, "done"

# filter strand data
bedAll_file = tmpPth+"tmp_bedA"+rTag+".txt"
bedAll_open = open(bedAll_file, 'w')

bedPlus_open = open(bedPlus_file,'r')
bedMinus_open = open(bedMinus_file,'r')

print >> sys.stderr, "filtering strand data"

for line in bedPlus_open:
	j = line.rstrip().split(")@@(")
	k = "\t".join(j)
	l = k.split("\t")
	go_line = l[:len(l)-2] + l[len(l)-1:]
	if j[2] == "+":
		bedAll_open.write("\t".join(go_line) + "\n")

bedPlus_open.close()

for line in bedMinus_open:
	j = line.rstrip().split(")@@(")
	k = "\t".join(j)
	l = k.split("\t")
	go_line = l[:len(l)-2] + l[len(l)-1:]
	if j[2] == "-":
		bedAll_open.write("\t".join(go_line) + "\n")

bedMinus_open.close()
bedAll_open.close()

print >> sys.stderr, "done"

# Create final output by bed sort
final_out = open(args.out_bed,'w')

print >> sys.stderr, "eliminating and sorting output"
inFile = open(bedAll_file,'r')
slimFile = tmpPth+"tmp_bedB"+rTag+".txt"
slimOut = open(slimFile,'w')
for line in inFile:
	data = line.rstrip().split('\t')
	score = sum([float(i) for i in data[12].split(',')])
	if score > 0.0:
		outLine = data[:len(data)]
		slimOut.write('\t'.join(outLine)+'\n')
inFile.close()
slimOut.close()
sortFile = subprocess.check_call(["sort","-k1,1","-k2,2n","-i",slimFile],stdout=final_out)

final_out.close()
print >> sys.stderr, "done"

# remove temporary directory
subprocess.check_call(["rm","-r",tmpPth])
