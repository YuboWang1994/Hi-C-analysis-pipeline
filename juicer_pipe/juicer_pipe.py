import os
import fire
from multiprocessing import Process
import time

seqkit = '{juicer_pipe}/seqkit'
bwa_mem = '{bwa-mem2_path}'
chimeric_blacklist='{juicer_pipe}/chimeric_blacklist.awk'
fragment = '{juicer_pipe}/fragment.pl'
dups = '{juicer_pipe}/dups.awk'
samtools ='{samtools_path}'
perl = '{perl_path}'
pigz = '{pigz_path}'

def main(fq1,fq2,sample,index,enzyme_txt,out,threads=10,cpus=2):
    start = time.time()
    out = os.path.abspath(out)
    fq1,fq2 = os.path.abspath(fq1),os.path.abspath(fq2)
    tmp = f'{out}/tmp'
    mkdir(out,tmp)
    
    split_num = threads//cpus
    minilist = sub_fq(fq1,fq2,split_num,tmp,threads)
    process = []
    fraglist = []
    for i in range(len(minilist)):
        tmp_fq1,tmp_fq2 = minilist[i].split()[0],minilist[i].split()[1]
        tmp_sample = f'{sample}_{i}'
        tmp_out = f'{out}/{tmp_sample}'
        mkdir(tmp_out)
        fraglist.append(f'{tmp_out}/{tmp_sample}_norm_frag.txt')
        t = Process(target=juicer_split,args=(tmp_fq1,tmp_fq2,tmp_sample,index,enzyme_txt,tmp_out,cpus))
        process.append(t)
    for i in process:
        i.start()
    for i in process:
        i.join()
    juicer_merge_dedup(fraglist,sample,out)
    print(time.time()-start)

def juicer_split(fq1,fq2,sample,index,enzyme_txt,out,threads=1):
    os.chdir(out)
    print(f'{bwa_mem} mem -SP5M -t {threads} {index} {fq1} {fq2} |{samtools}  view -bS -@ {threads} > {out}/{sample}.bam')
    os.system(f'{bwa_mem} mem -SP5M -t {threads} {index} {fq1} {fq2} |{samtools}  view -bS -@ {threads} > {out}/{sample}.bam')
    print(f'{samtools} view -@ {threads} {out}/{sample}.bam |LC_ALL=C awk -v "fname1"="{out}/{sample}_norm.txt" -v  "fname2"="{out}/{sample}_abnrom.sam" -v "fname3"="{out}/{sample}_unmapped.sam"  -f {chimeric_blacklist}')
    os.system(f'{samtools} view -@ {threads} {out}/{sample}.bam |LC_ALL=C awk -v "fname1"="{out}/{sample}_norm.txt" -v  "fname2"="{out}/{sample}_abnrom.sam" -v "fname3"="{out}/{sample}_unmapped.sam"  -f {chimeric_blacklist}')
    print(f'{perl} {fragment} {out}/{sample}_norm.txt {out}/{sample}_norm_frag.txt {enzyme_txt}')
    os.system(f'{perl} {fragment} {out}/{sample}_norm.txt {out}/{sample}_norm_frag.txt {enzyme_txt}')
    return f'{out}/{sample}_norm_frag.txt'

def juicer_merge_dedup(fraglist,sample,out):
    os.chdir(out)
    print(f'LC_ALL=C sort -T {out}/tmp  -k2,2V -k6,6V -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n -m {" ".join(fraglist)} > {sample}_sort.txt')
    os.system(f'LC_ALL=C sort -T {out}/tmp  -k2,2V -k6,6V -k4,4n -k8,8n -k1,1n -k5,5n -k3,3n -m {" ".join(fraglist)} > {sample}_sort.txt')
    print(f'LC_ALL=C awk -f {dups} -v name={out}/ {out}/{sample}_sort.txt')
    os.system(f'LC_ALL=C awk -f {dups} -v name={out}/{sample}_ {out}/{sample}_sort.txt')

def sub_fq(fq1,fq2,split_num,out,thread):
    #os.system(f'{seqkit} split2  -p {split_num} -1 {fq1} -2 {fq2}  -j {thread} -w 150 -O {out} -f ')
    prefix1 = os.path.basename(fq1).split('.')[0]
    prefix2 = os.path.basename(fq2).split('.')[0]
    tail = os.path.basename(fq1).split('.')[1]+'.'+os.path.basename(fq1).split('.')[2]
    minilist = []
    for i in range(1,split_num+1):
        minilist.append(f'{out}/{prefix1}.part_{str(i).rjust(3,"0")}.{tail} {out}/{prefix2}.part_{str(i).rjust(3,"0")}.{tail}')
    return minilist

def mkdir(*args):
    for i in args:
        if os.path.exists(i) == False:
            os.makedirs(i)

if __name__ == '__main__':
    fire.Fire(main) 
