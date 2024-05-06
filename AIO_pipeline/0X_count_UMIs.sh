######################################
## Count number of unique UMIs in each sample for each target 
######################################
## Orient reads by the 3' adapter 

mkdir /home/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p_noTrim
mkdir /home/mkramer/projects/target_capture/ruby_round2/data/4_cutadapt_trim_5p_3p_adapters/3p_noTrim/logFiles

condor_submit /home/mkramer/projects/target_capture/ruby_round2/jobFiles/3_cutadapt_trim_3p_adapters.noTrim.job

