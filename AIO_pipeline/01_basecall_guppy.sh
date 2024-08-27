######################################
## 1 Base call raw Nanopore data using guppy basecaller
######################################
/apps/src/ont-guppy_5.0.11/bin/guppy_basecaller -i /shares/slotkin/users/mkramer/ruby_round2/0_fast5/fc3/fc3 -s /scratch/TEMP_wkdir/TEMP_outf5_dir --device "cuda:all:100%" \
  --flowcell FLO-MIN106 \
  --kit SQK-LSK110 -q 0 --min_qscore 7 \
  --recursive --fast5_out --chunks_per_runner 1024 \
   --chunk_size 1024 --gpu_runners_per_device 1 --num_callers 1
