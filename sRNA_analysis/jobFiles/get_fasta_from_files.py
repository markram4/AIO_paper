import sys

inFile = open(sys.argv[1],'r') #"ath.deg"
for line in inFile:
    line = line.rstrip()
    #print(line)
    if line.startswith("ath"):
        (mirna,target,score,pval) = line.split("\t")
        print(">"+mirna+"-"+target)
    elif line.startswith("tar"):
        (tar,p5,seq,p3)=line.split(" ")
        newseq=seq.replace("U", "T")

        print(newseq)
