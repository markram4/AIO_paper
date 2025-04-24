##########################################
## 5 Map to capture sequences
##########################################
## Wt 35S::RUBY
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_35s.transgene.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa

## This file will map then split into unmapped and mapped directories and will index bam files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt ; done 

bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round2 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round2/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa

########
## ZmUBQ::RUBY in Arabidopsis
########
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_atubq.transgene.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.AtUBQ_RUBY.fa

## This file will map then split into unmapped and mapped directories and will index bam files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round3/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round3/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt ; done 

sbatch /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round3 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round3/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.ZmUBQ_RUBY.fa

########
## 35S::RUBY , AtUBQ::RUBY, LsUBQ::RUBY in lettuce
########
cat At_array.v2.fa ruby_atubq.transgene.fa > At_array.v2.AtUBQ_RUBY.fa
cat At_array.v2.fa ruby_LsUBQ.transgene.fa > At_array.v2.LsUBQ_RUBY.fa
cat At_array.v2.fa ruby_35s.transgene.fa > At_array.v2.RUBY.fa

## This file will map then split into unmapped and mapped directories and will index bam files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*35S*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.35S.txt ; done 
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*At*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.AtUBQ.txt ; done 
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*LsU*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.LsUBQ.txt ; done 

bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round4 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.35S.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa
bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round4 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.AtUBQ.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.AtUBQ_RUBY.fa
bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round4 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round4/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.LsUBQ.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.LsUBQ_RUBY.fa

########
## 35S::RUBY Mutants
########
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_35s.transgene.pMCK007.HtoQ.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.pMCK007.HtoQ.fa
cat /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.fa /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/ruby_35s.transgene.pMCK010.2bpins.fa > /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.pMCK010.2bpins.fa

## This file will map then split into unmapped and mapped directories and will index bam files
for file in /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/*fastq; do echo $(basename $file .fastq) >> /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.txt ; done 

## separate by which transgene version they should map to
cd /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq

grep SL2621 samples.txt > samples.SL2621.txt
grep SL2622 samples.txt > samples.SL2622.txt
grep MK samples.txt > samples.Wt.txt

bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round6 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.Wt.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.fa
bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round6 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.SL2621.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.pMCK007.HtoQ.fa
bash /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/commonFiles/generate_target_mapping_scripts.sh ruby_round6 /cluster/pixstor/slotkinr-lab/mkramer/projects/target_capture/ruby_round6/data/6_minimap/1_map_to_rRNA_tRNA/output/rRNA_free_fastq/samples.SL2622.txt /cluster/pixstor/slotkinr-lab/mkramer/annotations/arabidopsis/At_v2/At_array.v2.RUBY.pMCK010.2bpins.fa
