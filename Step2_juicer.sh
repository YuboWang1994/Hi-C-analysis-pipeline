mkdir -p {PATH}/01_juicer
cd {PATH}/01_juicer

ln -s {genome} genome.fasta
{bwa-mem2_path} index -p genome.fasta genome.fasta
{samtools_path} faidx genome.fasta
python2 {juicer_pipe}/generate_site_positions.py {enzyme} genome genome.fasta

cat genome.gtf |awk '$3=="transcript"'|awk '{print $1"\t"$4"\t"$5"\t"$12}'|sed 's/"//g'|sed 's/;//g' > genome.gene.bed


# 每个样本执行以下操作
mkdir -p {PATH}/01_juicer/{sample}
cd {PATH}/01_juicer/{sample}

mkdir -p tmp
python {juicer_pipe}/juicer_pipe.py {sample}.R1.clean.fq.gz {sample}.R2.clean.fq.gz {sample} {PATH}/01_juicer/genome.fasta {PATH}/01_juicer/genome_MboI.txt {PATH}/01_juicer/{sample} --threads 20

mkdir -p ref_split
{script}/faSplit byname {PATH}/01_juicer/genome.fasta ref_split/

python3 {HIC_PIPE}/N03_map/src/juicer/mnd2hiclib1.py {name}_mq30.txt {name}-MboI_refined.frag genome.fasta.fai

/public/home/zyy2020/miniconda3/envs/py2/bin/python2 {HIC_PIPE}/N03_map/src/icepipe/get_hm_bychr.py -rf {name}-MboI_refined.frag -hm 100,200,400 -rs ref_split -e MboI -g {name} -s {name}
/public/home/zyy2020/miniconda3/envs/py2/bin/python2 {HIC_PIPE}/N03_map/src/icepipe/04_iterativeCorrection.py ref_split 100,200,400 {name} MboI
/public/home/zyy2020/miniconda3/envs/py2/bin/python2 {HIC_PIPE}/N03_map/src/icepipe/scale_plot.py --ref-split ref_split {name} --plot-reslu 100,200,400

python3 {HIC_PIPE}/N03_map/src/icepipe/get_hm_bychr.py -rf {name}-MboI_refined.frag -b 10,20,50,100,200 -rs ref_split  -e MboI -g {name} -s {name}
python3 {HIC_PIPE}/N03_map/src/icepipe/bychr_plot_new.py --bychr-dir {name} --plot-reslu 10,20,50,100,200 --fai genome.fasta.fai --outdir bychr
python3 {HIC_PIPE}/N03_map/src/frag2distance.py {name}-MboI_refined.frag 200,400,800,1200,1600,3200 {name} {name}/picture/distance
python3 {HIC_PIPE}/N03_map/src/frag2bin_depth.py {name}-MboI_refined.frag genome.fasta.fai {name} {name}/picture/resolution

