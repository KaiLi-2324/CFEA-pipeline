import os
import glob
import argparse
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is designed to trim and align RRBS/WGBS fastq files
with trim_galore and bismark
"""


def get_files(path_fastq):
    """
    directory to fastq files
    :param path_fastq:
    :return:
    """
    print(f"Start reading from {path_fastq}")
    fastq_files = glob.glob(f"{path_fastq}/*.fastq")
    return fastq_files


def trim_map_file(fastq, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir):
    """
    trim fastq files with trim_galore, and align them with bismark
    :param fastq: directory to fastq files
    :param clean_dir: directory to put trimmed fastq files
    :param log_dir: directory to put log files
    :param bismark_index: directory to bismark index files
    :param bowtie2_path: directory to bowtie2 executable
    :param mapping_dir: directory to put aligned bam files
    :param mapping_result_dir: directory to put bismark log files
    :return: none
    """
    prefix = os.path.basename(fastq).split(".")[0]
    this_sample_clean_out = f"{clean_dir}/{prefix}"
    this_sample_mapping_out = f"{mapping_dir}/{prefix}"
    this_sample_trim_log = f"{log_dir}/{prefix}.trim.log"
    this_sample_bismark_log = f"{mapping_result_dir}/{prefix}/{prefix}.bismark.txt"
    this_sample_error_log = f"{log_dir}/{prefix}.trim.error.log"
    if not os.path.exists(f"{mapping_result_dir}/{prefix}"):
        os.makedirs(f"{mapping_result_dir}/{prefix}")
    if not os.path.exists(this_sample_clean_out):
        os.makedirs(this_sample_clean_out)
    if not os.path.exists(this_sample_mapping_out):
        os.makedirs(this_sample_mapping_out)
    if fastq.endswith(".sra_1.fastq"):
        fastq_file_1 = f"{os.path.dirname(fastq)}/{prefix}.sra_1.fastq"
        fastq_file_2 = f"{os.path.dirname(fastq)}/{prefix}.sra_2.fastq"
        cmd_trim = f"trim_galore --rrbs --illumina --phred33 " \
                   f"--paired {fastq_file_1} {fastq_file_2} -o {this_sample_clean_out} > {this_sample_trim_log} 2>&1"
        clean_file_1 = f"{this_sample_clean_out}/{prefix}.sra_1_val_1.fq"
        clean_file_2 = f"{this_sample_clean_out}/{prefix}.sra_2_val_2.fq"
        cmd_bismark = f"bismark {bismark_index} --path_to_bowtie {bowtie2_path} --bowtie2 " \
                      f"-1 {clean_file_1} -2 {clean_file_2} -o {this_sample_mapping_out} " \
                      f"--temp_dir {log_dir} > {this_sample_bismark_log} 2>&1"
        bam_file = f"{this_sample_mapping_out}/{prefix}.sra_1_val_1_bismark_bt2_pe.bam"
    else:
        cmd_trim = f"trim_galore --rrbs --illumina --phred33 " \
                   f"{fastq} -o {this_sample_clean_out} > {this_sample_trim_log} 2>&1"
        clean_file = f"{this_sample_clean_out}/{prefix}.sra_trimmed.fq"
        cmd_bismark = f"bismark {bismark_index} --path_to_bowtie {bowtie2_path} " \
                      f"--bowtie2 {clean_file} -o {this_sample_mapping_out} " \
                      f"--temp_dir {log_dir} > {this_sample_bismark_log} 2>&1"
        bam_file = f"{this_sample_mapping_out}/{prefix}.sra_trimmed_bismark_bt2.bam"
    bam_qc_dir = f"{mapping_result_dir}/bam_qc/{prefix}"
    if not os.path.exists(bam_qc_dir):
        os.makedirs(bam_qc_dir)
    cmd_bam_qc = f"fastqc -o {bam_qc_dir} -q {bam_file}"
    try:
        print(f"  Trimming: {prefix}")
        subprocess.check_output(cmd_trim, shell=True)
        print(f"  Aligning: {prefix}")
        subprocess.check_output(cmd_bismark, shell=True)
        print(f"  Bam qc: {prefix}")
        subprocess.check_output(cmd_bam_qc, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(fastq_files, num_process, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir):
    """
    to run faster, we create a multiprocessing pool here
    :param fastq_files: same as above, same below
    :param num_process:
    :param clean_dir:
    :param log_dir:
    :param bismark_index:
    :param bowtie2_path:
    :param mapping_dir:
    :param mapping_result_dir:
    :return:
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for fastq in fastq_files:
        if not fastq.endswith(".sra_2.fastq"):
            pool.apply_async(trim_map_file, (fastq, clean_dir, log_dir, bismark_index, bowtie2_path, mapping_dir, mapping_result_dir,))
    pool.close()
    pool.join()


def get_multiqc_file(mapping_dir, mapping_result_dir):
    """
    collect qc file for multiQC
    :param mapping_dir:
    :param mapping_result_dir:
    :return:
    """
    print(f"Start reading from {mapping_dir}")
    bismark_report_files = glob.glob(f"{mapping_dir}/*/*_report.txt")
    bismark_report_multiqc_path = f"{mapping_result_dir}/mapping_multiqc_sample.txt"
    print(f"Start reading from {mapping_result_dir}/bam_qc")
    bam_qc_files = glob.glob(f"{mapping_result_dir}/bam_qc/*/*.zip")
    bam_qc_multiqc_path = f"{mapping_result_dir}/bam_qc_multiqc_sample.txt"
    with open(bismark_report_multiqc_path, "w+") as out:
        out.write("\n".join(bismark_report_files) + "\n")
    with open(bam_qc_multiqc_path, "w+") as out:
        out.write("\n".join(bam_qc_files) + "\n")
    return bismark_report_multiqc_path, bam_qc_multiqc_path


def multiqc(bismark_report_multiqc_path, bam_qc_multiqc_path, n, mapping_result_dir):
    """
    run multiQC
    :param bismark_report_multiqc_path:
    :param bam_qc_multiqc_path:
    :param n:
    :param mapping_result_dir:
    :return:
    """
    print("Start multiqc")
    n_bismark = f"{n}_bismark"
    n_bam_qc = f"{n}_bam_qc"
    cmd_bismark = f"multiqc -n {n_bismark} -s -o {mapping_result_dir} -l {bismark_report_multiqc_path}"
    cmd_bam_qc = f"multiqc -n {n_bam_qc} -s -o {mapping_result_dir} -l {bam_qc_multiqc_path}"
    try:
        subprocess.check_output(cmd_bismark, shell=True)
        subprocess.check_output(cmd_bam_qc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="path_fastq")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="mapping_result_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-b", action="store", dest="mapping_dir", help="path to put bam files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-t", action="store", dest="bowtie2_path", help="path to bowtie2 executable")
    parser.add_argument("-d", action="store", dest="path_index", help="path to bismark/bowtie2 indexes")
    parser.add_argument("-c", action="store", dest="clean_dir", help="path to put trimmed fastq files")
    results = parser.parse_args()
    if not all([results.path_fastq, results.num_process, results.mapping_result_dir, results.log_dir,
                results.n, results.mapping_dir, results.bowtie2_path, results.path_index, results.clean_dir]):
        print("too few arguments, type -h for more information!")
        exit(-1)
    path_fastq = results.path_fastq if not results.path_fastq.endswith("/") else results.path_fastq.rstrip("/")
    num_process = results.num_process
    mapping_result_dir = results.mapping_result_dir if not results.mapping_result_dir.endswith(
        "/") else results.mapping_result_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    n = results.n
    bowtie2_path = results.bowtie2_path if not results.bowtie2_path.endswith("/") else results.bowtie2_path.rstrip("/")
    path_index = results.path_index if not results.path_index.endswith("/") else results.path_index.rstrip("/")
    clean_dir = results.clean_dir if not results.clean_dir.endswith("/") else results.clean_dir.rstrip("/")
    for each_dir in [mapping_result_dir, log_dir, mapping_dir, clean_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    fastq_files = get_files(path_fastq)
    multi_process_run(fastq_files, num_process, clean_dir, log_dir,
                      path_index, bowtie2_path, mapping_dir, mapping_result_dir)
    bismark_report_multiqc_path, bam_qc_multiqc_path = get_multiqc_file(mapping_dir, mapping_result_dir)
    multiqc(bismark_report_multiqc_path, bam_qc_multiqc_path, n, mapping_result_dir)


if __name__ == '__main__':
    main()
