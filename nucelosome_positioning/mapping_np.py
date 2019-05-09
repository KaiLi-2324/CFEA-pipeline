import re
import os
import glob
import argparse
import subprocess
import multiprocessing

"""
@Author: likai
@Email: likai@wibe.ac.cn
This script is developed for mapping of nucleosome positioning data, be sure that 
picard samtools, fastqc and multiqc are in your local path, one way is to use
conda. 
conda install -n py2 picard, samtools, fastqc, multiqc, bwa
source activate /share/pub/lik/soft/miniconda3/envs/py2
"""


def get_files(path_fastq):
    """
    get all raw fastq files
    :param path_fastq: directory to raw fastq files
    :return: list with all fastq files
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq.gz")
    return fastq_files


def mapping_bwa(fastq_dir, sample, bwa_index, mapping_out_dir, bam_dir, log_dir):
    # we use bwa to map the raw fastq sequences onto hg19.fa. After alignment, samtools
    # was used to convert sam files into bam files. Finally, picard was adopted to remove
    # duplicated sequences from the aligned bam files
    fastq1 = f"{fastq_dir}/{sample}_1.fastq.gz"
    fastq2 = f"{fastq_dir}/{sample}_2.fastq.gz"
    this_sample_mapping_out = f"{bam_dir}/{sample}"
    bam1 = f"{this_sample_mapping_out}/{sample}_1.sai"
    bam2 = f"{this_sample_mapping_out}/{sample}_2.sai"
    sam = f"{this_sample_mapping_out}/{sample}.sam"
    sort_bam = f"{this_sample_mapping_out}/{sample}.sort.bam"
    dedup_bam_file = f"{this_sample_mapping_out}/{sample}.dedup.bam"
    metrics_picard_file = f"{mapping_out_dir}/picard_out/{sample}.dedup.txt"
    mapping1_log = f"{log_dir}/{sample}.bwa1.log"
    mapping2_log = f"{log_dir}/{sample}.bwa2.log"
    sam_log = f"{log_dir}/{sample}.sam.log"
    this_sample_picard_log = f"{log_dir}/{sample}.picard.log"
    this_sample_error_log = f"{log_dir}/{sample}.bwa.error"
    cmd1 = f"bwa aln -t 40 -f {bam1} {bwa_index} {fastq1} > {mapping1_log} 2>&1"
    cmd2 = f"bwa aln -t 40 -f {bam2} {bwa_index} {fastq2} > {mapping2_log} 2>&1"
    cmd_sam = f"bwa sampe -f {sam} {bwa_index} {bam1} {bam2} {fastq1} {fastq2} > {sam_log} 2>&1"
    cmd_sort = f"samtools view -h {sam} | samtools sort -@ 10 -o {sort_bam}"
    cmd_picard_dedup = f"picard MarkDuplicates INPUT={sort_bam} OUTPUT={dedup_bam_file} " \
                       f"METRICS_FILE={metrics_picard_file} VALIDATION_STRINGENCY=LENIENT " \
                       f"ASSUME_SORTED=true REMOVE_DUPLICATES=true > {this_sample_picard_log} 2>&1"
    try:
        if not os.path.exists(this_sample_mapping_out):
            os.makedirs(this_sample_mapping_out)
        if not os.path.exists(f"{mapping_out_dir}/picard_out"):
            os.makedirs(f"{mapping_out_dir}/picard_out")
        print(f"Aligning:  {sample}")
        subprocess.check_output(cmd1, shell=True)
        subprocess.check_output(cmd2, shell=True)
        subprocess.check_output(cmd_sam, shell=True)
        print(f"Sorting:  {sample}")
        subprocess.check_output(cmd_sort, shell=True)
        print(f"Removing duplicates:  {sample}")
        subprocess.check_output(cmd_picard_dedup, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "w+") as log:
            log.write(str(e) + "\n")


def multi_process_run(fastq_files, fastq_dir, bwa_index, mapping_out_dir, bam_dir, log_dir, num_process):
    # To run faster, we use the multiprocessing module to create a processing pool
    # you can control the number of processes in the pool with -p parameter
    pool = multiprocessing.Pool(processes=int(num_process))
    total_samples = []
    for each_fastq in fastq_files:
        total_samples.append(re.sub("_\d", "", os.path.basename(each_fastq).split(".")[0]))
    total_samples = list(set(total_samples))
    for sample in total_samples:
        pool.apply_async(mapping_bwa, (fastq_dir, sample, bwa_index, mapping_out_dir, bam_dir, log_dir,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="mapping_out_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="bam_dir", help="path to put bam files")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bwa indexes")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_out_dir,
                results.log_dir, results.bam_dir, results.path_index]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_out_dir = results.mapping_out_dir if not results.mapping_out_dir.endswith(
        "/") else results.mapping_out_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    bam_dir = results.bam_dir if not results.bam_dir.endswith("/") else results.bam_dir.rstrip("/")
    path_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    for each_dir in [mapping_out_dir, log_dir, bam_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, path_fastq, path_index, mapping_out_dir, bam_dir, log_dir, num_process)


if __name__ == '__main__':
    main()
