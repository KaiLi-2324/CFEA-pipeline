import os
import glob
import argparse
import subprocess
import multiprocessing


"""
this script is written to call peaks from the deduplicated bam files with macs2 software
__author__ = likai
__email__ = likai@wibe.ac.cn
"""


def get_bam_files(mapping_dir):
    """
    get all the deduplicated bam file from the mapping directory if there are any
    :param mapping_dir: directory to the deduplicated bam files
    :return: python list with all targeted bam files
    """
    print(f"Start reading from {mapping_dir}")
    bam_files = glob.glob(f"{mapping_dir}/*.dedup.bam")
    if not bam_files:
        bam_files = glob.glob(f"{mapping_dir}/*.sort.bam")
    return bam_files


def call_peak(bam, peak_dir, log_dir):
    """
    we use macs2 to call peak from deduplicated bam files
    :param bam: path to each deduplicated bam files
    :param peak_dir: directory to put results of macs2
    :param log_dir: directory to put log files in each step
    :return: none
    """
    prefix = os.path.basename(bam).split(".")[0]
    print(f"    Calling Peak:  {prefix}")
    this_sample_peak_out = f"{peak_dir}/{prefix}"
    peak_log = f"{log_dir}/{prefix}.macs2.log"
    error_log = f"{log_dir}/{prefix}.macs2.error.log"
    if not os.path.exists(this_sample_peak_out):
        os.makedirs(this_sample_peak_out)
    cmd_macs2 = f"macs2 callpeak -t {bam} -g hs --nomodel " \
                f"--extsize 200 -n {prefix} --outdir {this_sample_peak_out} --call-summits > {peak_log} 2>&1"
    try:
        subprocess.check_output(cmd_macs2, shell=True)
    except Exception as e:
        print(e)
        with open(error_log, "w+") as out:
            out.write(str(e) + "\n")


def multi_process_run(bam_files, num_process, peak_dir, log_dir):
    """
    to run faster, we create a processing pool with a certain number of processes, the params
    are the same to function call_peak
    :param bam_files:
    :param num_process:
    :param peak_dir:
    :param log_dir:
    :return:
    """
    pool = multiprocessing.Pool(processes=int(num_process))
    for bam in bam_files:
        pool.apply_async(call_peak, (bam, peak_dir, log_dir,))
    pool.close()
    pool.join()


def get_multiqc_files(peak_dir, peak_out_dir):
    """
    we use multiQC software to collect the qc results for each peak files by macs2
    :param peak_dir: directory to put results of macs2
    :param peak_out_dir:
    :return:
    """
    print(f"Start reading from {peak_dir}")
    peak_files = glob.glob(f"{peak_dir}/*/*_peaks.xls")
    peak_qc_path = f"{peak_out_dir}/peaks_multiqc_sample.txt"
    with open(peak_qc_path, "w+") as out:
        out.write("\n".join(peak_files) + "\n")
    return peak_qc_path


def multiqc(peak_qc_path, peak_out_dir, n):
    """
    run multiQC to collect statics for peak files of each deduplicated bam files by macs2
    :param peak_qc_path: directory to the peak file called by macs2
    :param peak_out_dir: directory to put multiqc result
    :param n:
    :return:
    """
    print("Start multiqc")
    n_peak = f"{n}_peak"
    cmd_peak_multiqc = f"multiqc -n {n_peak} -s -o {peak_out_dir} -l {peak_qc_path}"
    try:
        subprocess.check_output(cmd_peak_multiqc, shell=True)
    except Exception as e:
        print(e)


def main():
    parser = argparse.ArgumentParser(description="A program to call peaks from 5hmc bam files")
    parser.add_argument("-i", action="store", dest="mapping_dir")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="peak_dir", help="path to put called peaks in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    parser.add_argument("-n", action="store", dest="n", help="project name to use in MultiQC output[bioproject.13]")
    parser.add_argument("-o", action="store", dest="peak_out_dir", help="directory to put peak qc results")
    results = parser.parse_args()
    if not all([results.mapping_dir, results.num_process,
                results.peak_dir, results.log_dir, results.n, results.peak_out_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    num_process = results.num_process
    peak_dir = results.peak_dir if not results.peak_dir.endswith("/") else results.peak_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    n = results.n
    peak_out_dir = results.peak_out_dir if not results.peak_out_dir.endswith("/") else results.peak_out_dir.rstrip("/")
    for each_dir in [peak_dir, log_dir, peak_out_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    bam_files = get_bam_files(mapping_dir)
    multi_process_run(bam_files, num_process, peak_dir, log_dir)
    peak_qc_path = get_multiqc_files(peak_dir, peak_out_dir)
    multiqc(peak_qc_path, peak_out_dir, n)


if __name__ == '__main__':
    main()
