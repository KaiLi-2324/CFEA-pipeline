import os
import glob
import argparse
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is written to align fastq files by MCTA-seq. Before aligning two reads separately, 
we use three custom perl scripts to trim fastq files.
"""


def get_files(path_fastq):
    """
    get all the fastq files from MCTA-seq, these files need to be ended with .gz
    :param path_fastq:
    :return: python list with all gzipped fastq files
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*/*.fastq.gz")
    return fastq_files


def trim_map_file(fastq, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index, bowtie1_path, mapping_result_dir):
    """
    we first use QCnolane.pl to trim both fastq files, then use rm_firstx_leny.pl
    to remove forward xbp, backward zbp, remain more than ybp in the trimmed fastq
    files. Bismark is used to align read1 and read2 separately.
    :param fastq: path to the input gzipped fastq file
    :param path_fastq: directory to all fastq files
    :param clean_dir: directory to put trimmed fastq files
    :param log_dir: directory to put log files
    :param mapping_dir: directory to put aligned bam files
    :param bismark_index: path to bismark index files
    :param bowtie1_path: path to bowtie1
    :param mapping_result_dir: directory to put bismark log files
    :return: none
    """
    prefix = os.path.basename(fastq).split("_")[0]
    this_sample_error_log = f"{log_dir}/{prefix}.error.log"
    cmd_trim = f"perl QCnolane.pl --indir {path_fastq} --outdir {clean_dir} --sample {prefix}"
    clean_fq_1 = f"{clean_dir}/{prefix}/{prefix}.1.clean.fq.gz"
    clean_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.clean.fq.gz"
    filter_fq_1 = f"{clean_dir}/{prefix}/{prefix}.1.filter.fq.gz"
    filter_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.filter.fq.gz"
    nonCG_fq_2 = f"{clean_dir}/{prefix}/{prefix}.2.nonCG.fq.gz"
    bismark_output_1 = f"{mapping_dir}/{prefix}_1"
    bismark_output_2 = f"{mapping_dir}/{prefix}_2"
    if not os.path.exists(bismark_output_1):
        os.makedirs(bismark_output_1)
    if not os.path.exists(bismark_output_2):
        os.makedirs(bismark_output_2)
    cmd_filter_2 = f"perl rm_firstx_leny.pl {clean_fq_2} {filter_fq_2} 6 20 12"
    cmd_filter_1 = f"perl rm_firstx_leny.pl {clean_fq_1} {filter_fq_1} 12 20 6"
    bismark_log_1 = f"{mapping_result_dir}/{prefix}/{prefix}.1.mapping.txt"
    bismark_log_2 = f"{mapping_result_dir}/{prefix}/{prefix}.2.mapping.txt"
    if not os.path.exists(f"{mapping_result_dir}/{prefix}"):
        os.makedirs(f"{mapping_result_dir}/{prefix}")
    cmd_ch3 = f"perl ch3deleate.pl {filter_fq_2} {nonCG_fq_2}"
    # we use bismark to align each read separately
    cmd_bismark_1 = f"bismark --bowtie1 --non_directional --fastq --phred33-quals --temp_dir {log_dir} " \
                    f"--path_to_bowtie {bowtie1_path} --output_dir {bismark_output_1} {bismark_index}" \
                    f" {filter_fq_1} > {bismark_log_1} 2>&1"
    cmd_bismark_2 = f"bismark --bowtie1 --non_directional --fastq --phred33-quals --temp_dir {log_dir} " \
                    f"--path_to_bowtie {bowtie1_path} --output_dir {bismark_output_2} {bismark_index} " \
                    f"{nonCG_fq_2} > {bismark_log_2} 2>&1"
    bam1 = f"{bismark_output_1}/{prefix}.1.filter_bismark.bam"
    bam2 = f"{bismark_output_2}/{prefix}.2.nonCG_bismark.bam"
    cmd_samtools_sort_1 = f"samtools sort -@ 20 {bam1} -o {os.path.dirname(bam1)}/{prefix}.1.sort.bam"
    cmd_samtools_sort_2 = f"samtools sort -@ 20 {bam2} -o {os.path.dirname(bam2)}/{prefix}.2.sort.bam"
    cmd_index_1 = f"samtools index {os.path.dirname(bam1)}/{prefix}.1.sort.bam"
    cmd_index_2 = f"samtools index {os.path.dirname(bam2)}/{prefix}.2.sort.bam"
    try:
        print(f"    Trimming:  {prefix}")
        subprocess.check_output(cmd_trim, shell=True)
        print(f"    Filter fq 2:  {prefix}")
        subprocess.check_output(cmd_filter_2, shell=True)
        print(f"    Filter fq 1:  {prefix}")
        subprocess.check_output(cmd_filter_1, shell=True)
        print(f"    Filter out CG of fq 2:  {prefix}")
        subprocess.check_output(cmd_ch3, shell=True)
        print(f"    Aligning fq 1:  {prefix}")
        subprocess.check_output(cmd_bismark_1, shell=True)
        print(f"    Aligning fq 2:  {prefix}")
        subprocess.check_output(cmd_bismark_2, shell=True)
        print(f"    Sorting bam1:  {prefix}")
        subprocess.check_output(cmd_samtools_sort_1, shell=True)
        print(f"    Sorting bam2:  {prefix}")
        subprocess.check_output(cmd_samtools_sort_2, shell=True)
        print(f"    Indexing bam1:  {prefix}")
        subprocess.check_output(cmd_index_1, shell=True)
        print(f"    Indexing bam2:  {prefix}")
        subprocess.check_output(cmd_index_2, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "a") as out:
            out.write(str(e) + "\n")


def multi_process_run(fastq_files, num_process, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index,
                      bowtie1_path, mapping_result_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    print("Start checking files")
    for each_fastq in fastq_files:
        prefix = os.path.basename(each_fastq).split("_")[0]
        if not os.path.exists(f"{os.path.dirname(each_fastq)}/{prefix}_R1.fastq.gz") or not os.path.exists(
                f"{os.path.dirname(each_fastq)}/{prefix}_R2.fastq.gz"):
            print(f"{prefix} not completed!")
            exit(-1)
    for fastq in fastq_files:
        if not fastq.endswith("_R2.fastq.gz"):
            pool.apply_async(trim_map_file, (
                fastq, path_fastq, clean_dir, log_dir, mapping_dir, bismark_index, bowtie1_path, mapping_result_dir))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq for MCTA-seq")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=40)
    parser.add_argument("-m", action="store", dest="mapping_result_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="mapping_dir", help="path to put bam files")
    parser.add_argument("-t", action="store", dest="bowtie1_path", help="path to bowtie1 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bismark/bowtie1 indexes")
    parser.add_argument("-c", action="store", dest="clean_dir", help="path to put trimmed fastq files")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_result_dir, results.log_dir,
                results.mapping_dir, results.bowtie1_path, results.path_index, results.clean_dir]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_result_dir = results.mapping_result_dir if not results.mapping_result_dir.endswith(
        "/") else results.mapping_result_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    bowtie1_path = results.bowtie1_path if not results.bowtie1_path.endswith("/") else results.bowtie1_path.rstrip("/")
    bismark_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    clean_dir = results.clean_dir if not results.clean_dir.endswith("/") else results.clean_dir.rstrip("/")
    for each_dir in [mapping_result_dir, log_dir, mapping_dir, clean_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, path_fastq, clean_dir, log_dir,
                      mapping_dir, bismark_index, bowtie1_path, mapping_result_dir)


if __name__ == '__main__':
    main()
