# Step 1 Obtain clean data for each sample

{pipeline_path}/scripts/fastp -i {sample_raw_R1_data} -I {sample_raw_R2_data} -o {sample}_R1.fastq.gz -O {sample}_R2.fastq.gz -w 4 -j {sample}.json -h {sample}.html
