import os
import glob
import argparse
import subprocess
import multiprocessing

"""
This script is developed for mapping of 5hmc data, be sure that picard
samtools, fastqc and multiqc are in your local path, one way is to use
conda. 
conda install -n py2 picard samtools fastqc multiqc
source activate /share/pub/lik/soft/miniconda3/envs/py2
__author__ = likai
__email__ = likai@wibe.ac.cn
"""


def get_files(path_fastq):
    """
    get all raw fastq files
    :param path_fastq: directory to raw fastq files
    :return: list with all fastq files
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq")
    return fastq_files


def mapping_bowtie2(fastq, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir):
    """
    we first align the fastq files with bowtie2 aligner and then remove duplicates with
    picard tool set, after that, we rerun fastqc to get the quality of the deduplicated data
    :param fastq: path to each fastq file
    :param bowtie2_path: path to the bowtie2 aligner
    :param path_index: path to index files created by bowtie2 aligner
    :param mapping_out_dir: path to put the qc results of aligned bam files
    :param bam_dir: path to put aligned bam files
    :param log_dir: path to put log files in each step
    :return: none
    """
    prefix = os.path.basename(fastq).split(".")[0]
    fastq_file_1 = f"{os.path.dirname(fastq)}/{prefix}.sra_1.fastq"
    fastq_file_2 = f"{os.path.dirname(fastq)}/{prefix}.sra_2.fastq"
    sam_file = f"{bam_dir}/{prefix}.sam"
    this_sample_log = f"{log_dir}/{prefix}.mapping.error.log"
    this_sample_picard_log = f"{log_dir}/{prefix}.picard.log"
    bowtie2_out_result = f"{mapping_out_dir}/bowtie2_out/{prefix}.txt"
    if not os.path.exists(f"{mapping_out_dir}/bowtie2_out"):
        os.makedirs(f"{mapping_out_dir}/bowtie2_out")
    cmd_bowtie2 = f"{bowtie2_path} -p 30 --end-to-end -q " \
                  f"-x {path_index} -1 {fastq_file_1} -2 {fastq_file_2} -S {sam_file} > {bowtie2_out_result} 2>&1"
    sort_bam_file = f"{bam_dir}/{prefix}.sort.bam"
    cmd_sort_samtools = f"samtools view -h {sam_file} | samtools sort -@ 30 -o {sort_bam_file}"
    dedup_bam_file = f"{bam_dir}/{prefix}.dedup.bam"
    metrics_picard_file = f"{mapping_out_dir}/picard_out/{prefix}.dedup.txt"
    if not os.path.exists(f"{mapping_out_dir}/picard_out"):
        os.makedirs(f"{mapping_out_dir}/picard_out")
    cmd_picard_dedup = f"picard MarkDuplicates INPUT={sort_bam_file} OUTPUT={dedup_bam_file} " \
                       f"METRICS_FILE={metrics_picard_file} VALIDATION_STRINGENCY=LENIENT " \
                       f"ASSUME_SORTED=true REMOVE_DUPLICATES=true > {this_sample_picard_log} 2>&1"
    fastqc_dedup_out_dir = f"{mapping_out_dir}/fastqc_dedup/{prefix}"
    if not os.path.exists(fastqc_dedup_out_dir):
        os.makedirs(fastqc_dedup_out_dir)
    cmd_fastqc_dedup_bam = f"fastqc -q --extract -o {fastqc_dedup_out_dir} {dedup_bam_file}" # run fastqc for bam files
    try:
        print(f"  Aligning: {prefix}")
        subprocess.check_output(cmd_bowtie2, shell=True)
        print(f"  Sorting: {prefix}")
        subprocess.check_output(cmd_sort_samtools, shell=True)
        print(f"  Removing duplicates: {prefix}")
        subprocess.check_output(cmd_picard_dedup, shell=True)
        print(f"  Dedup QC: {prefix}")
        subprocess.check_output(cmd_fastqc_dedup_bam, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_log, "w+") as out:
            out.write(str(e) + "\n")
    

def multi_process_run(fastq_files, num_process, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir):
    """
    to run faster, we create a processing pool with a certain number of processes, the params
    are the same to function mapping_bowtie2
    :param fastq_files:
    :param num_process:
    :param bowtie2_path:
    :param path_index:
    :param mapping_out_dir:
    :param bam_dir:
    :param log_dir:
    :return: none
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for fastq in fastq_files:
        if not fastq.endswith(".sra_2.fastq"):
            pool.apply_async(mapping_bowtie2, (fastq, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir,))
    pool.close()
    pool.join()


def get_multiqc_file(mapping_out_dir):
    """
    we use multiQC software to collect the qc results for the aligned bam files and deduplicated bam files
    this function is designed to get all the aligned and deduplicated bam files from the mapping directory
    :param mapping_out_dir: directory to the qc results of the aligned bam files and deduplicated bam files
    :return: path to multiQC result file of aligned bam files and deduplicated bam files
    """
    print(f"Start reading from {mapping_out_dir}")
    bowtie2_out_result_dir = f"{mapping_out_dir}/bowtie2_out"
    bowtie2_out_files = glob.glob(f"{bowtie2_out_result_dir}/*.txt")
    bowtie2_out_multiqc_path = f"{mapping_out_dir}/bowtie2_multiqc_sample.txt"
    dedup_qc_dir = f"{mapping_out_dir}/fastqc_dedup"
    dedup_qc_zip_files = glob.glob(f"{dedup_qc_dir}/*/*.zip")
    dedup_qc_multiqc_path = f"{mapping_out_dir}/dedup_multiqc_sample.txt"
    with open(bowtie2_out_multiqc_path, "w+") as out:
        out.write("\n".join(bowtie2_out_files) + "\n")
    with open(dedup_qc_multiqc_path, "w+") as out:
        out.write("\n".join(dedup_qc_zip_files) + "\n")
    return bowtie2_out_multiqc_path, dedup_qc_multiqc_path


def multiqc(bowtie2_out_multiqc_path, dedup_qc_multiqc_path, n, mapping_out_dir):
    """
    run multiQC for the deduplicated bam files by picard
    :param bowtie2_out_multiqc_path:
    :param dedup_qc_multiqc_path:
    :param n:
    :param mapping_out_dir:
    :return:
    """
    print("Start multiqc")
    n_bowtie2 = f"{n}_bowtie2"
    # run multiQC for aligned bam files
    cmd_bowtie2_multiqc = f"multiqc -n {n_bowtie2} -s -o {mapping_out_dir} -l {bowtie2_out_multiqc_path}"
    n_dedup_qc = f"{n}_dedup_qc"
    # run multiQC for deduplicated bam files
    cmd_dedup_multiqc = f"multiqc -n {n_dedup_qc} -s -o {mapping_out_dir} -l {dedup_qc_multiqc_path}"
    try:
        subprocess.check_output(cmd_bowtie2_multiqc, shell=True)
        subprocess.check_output(cmd_dedup_multiqc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="mapping_out_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="bam_dir", help="path to put bam files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-t", action="store", dest="bowtie2_path", help="path to bowtie2 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bowtie2 indexes")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_out_dir, results.log_dir, 
                results.n, results.bam_dir, results.bowtie2_path, results.path_index]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_out_dir = results.mapping_out_dir if not results.mapping_out_dir.endswith("/") else results.mapping_out_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    bam_dir = results.bam_dir if not results.bam_dir.endswith("/") else results.bam_dir.rstrip("/")
    n = results.n
    bowtie2_path = results.bowtie2_path if not results.bowtie2_path.endswith("/") else results.bowtie2_path.rstrip("/")
    path_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    for each_dir in [mapping_out_dir, log_dir, bam_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, bowtie2_path, path_index, mapping_out_dir, bam_dir, log_dir)
    bowtie2_out_multiqc_path, dedup_qc_multiqc_path = get_multiqc_file(mapping_out_dir)
    multiqc(bowtie2_out_multiqc_path, dedup_qc_multiqc_path, n, mapping_out_dir)

    
if __name__ == '__main__':
    main()
