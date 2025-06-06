## Download Araport11 annotation from Ensembl https://plants.ensembl.org/info/data/ftp/index.html
#ensembl file has the protein coding genes gtf that can be successfully converted to bed12
bedparse gtf2bed --extraFields gene_biotype,gene_name /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.gtf | sort -k1,1 -k2,2n | awk 'BEGIN{FS="\t"}{OFS="\t"}{n=split($4,a,".")}{if ($13 != "miRNA" && a[1] != "AT3G17185") print $1,$2+1,$3,$4,a[1],$6,$7+1,$8,$9,$10,$11,$12,$13,$14; else if (a[1] == "AT3G17185") print $1,$2+1,$3,$4,a[1],$6,$2+1,$3,$9,$10,$11,$12,$13,$14}' > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.bed12

## Download Araport11 annotation from https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FAraport11_genome_release
# Araport has transposon information
# The ATNTE### genes sometimes have several smaller orfs annotated as transposable_element_gene (see AT2TE02465 on jbrowse for an example)
awk 'BEGIN{FS="\t"}{OFS="\t"}{n=split($1,a,"Chr")}{n=split($9,b,"\"")}{if ($3 == "transposable_element") print a[2], $4,$5,b[2],b[2],$7,$4,$5, "0", "1",$5-$4+1",","0,", $3, b[2]}'  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.gtf > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element.bed12
awk 'BEGIN{FS="\t"}{OFS="\t"}{n=split($1,a,"Chr")}{n=split($9,b,"\"")}{if ($3 == "transposable_element_gene") print a[2], $4,$5,b[2],b[2],$7,$4,$5, "0", "1",$5-$4+1",","0,", $3, b[2]}'  /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.gtf > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element_gene.bed12


## If transposable element gene overlaps 100% with a transposable element
intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element_gene.bed12 -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element.bed12 -wb -F 1 | awk '{OFS="\t"}{print $15,$16,$17,$18,$4,$20,$21,$22,$23,$24,$25,$26,$27,$28}' > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.TE_to_TEgene.bed12
intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element.bed12 -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.transposable_element_gene.bed12 -wb -F 1 -v > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.TE.bed12

cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.TE.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Araport_annotations/Araport11_GTF_genes_transposons.current.TE_to_TEgene.bed12 | sort -k1,1 -k2,2n > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.bed12

################################################
## Pull out genome coordinates for genome file
################################################
grep ">" /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa| awk '{OFS="\t"}{n=split($1,a,">")}{n=split(a[2],b,"_")}{n=split(b[2],c,"Chr")}{if($1 ~ "_") print c[2], b[3], b[4],b[1], a[2], "."}' | sort -k1,1 -k2,2n > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa.not_construct.bed
grep ">" /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa | awk '{OFS="\t"}{n=split($1,a,">")}{n=split(a[2],b,"_")}{if($1 ~ "_") print b[2], b[3], b[4],b[1], a[2], "."}' | sort -k1,1 -k2,2n > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa.not_construct.chr.bed

## Get gene information and annotation for genes in target genome file
intersectBed -a /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa.not_construct.bed -b /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.bed12 -nonamecheck -wao -F 0.9| awk '{OFS="\t"}{n=split($5,a,"_")}{if ($19 ~ /transposable/ || a[1] != $11) print $5, ($8-$2)+1, $9-$2,$10,$11,$12,($13-$2)+1,$14-$2,$15,$16,$17,$18,$19,a[1];else print $5, ($8-$2)+1, $9-$2,$10,$11,$12,($13-$2)+1,$14-$2,$15,$16,$17,$18,$19,$20} ' | sort -k1,1 -k2,2n  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.bed12 

# Make bed12 for construct sequences
awk -f /cluster/pixstor/slotkinr-lab/mkramer/annotations/soy_williams82/fastaToTbl /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa | awk 'BEGIN {FS="\t"}{OFS="\t"}{n=split($1,a," ")}{if ($1 !~ "_") print a[1],1,length($2)-1,a[1],a[1],"+",1,length($2)-1,0,1,length($2)-1",","0,","construct",a[1]}' > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.tbl.construct.bed12

## Combine with construct and transgenes
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.tbl.construct.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/35S_RUBY_transgene.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ZmUBQ_RUBY_transgene.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/AtUBQ_RUBY_transgene.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/LsUBQ_RUBY_transgene.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/35S_RUBY_transgene_mut.bed12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/35S_RUBY_transgene_ins.bed12 > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed12

cut -f1,2,3,4,5,6,7,8,9,10,11,12 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed12 > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.2.bed12

cut -f 1,2,3,4,14,6 /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed12 | awk '{OFS="\t"}{if ($6 != "") print $1,$2,$3,$4,$6,$5;else print $1,$2,$3,$4,$4,$5}' > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.Arabidopsis_thaliana.TAIR10.56.Araport11_transposable_element.current.targetChr.construct.transgenes.bed


################################################
## Create annotation for Riboseq
################################################
awk '{OFS=FS="\t"}{if ($1 != "Mt" && $1 != "Pt") print}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.gtf  > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.gtf

## Filter for only protein-coding genes in gtf file
awk '{OFS=FS="\t"}{if ($0 ~ "#") print; else if ($9 ~ /protein_coding/ || $1 == "35S_RUBY_transgene") print}' /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.gtf > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/ensembl_downloads/Arabidopsis_thaliana.TAIR10.56.35S_ruby_transgene.noMtPt.pcg_only.gtf
