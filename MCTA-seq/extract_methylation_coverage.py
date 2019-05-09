import os
import glob
import argparse
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
"""


def get_bam_files(path_bam):
    print(f"Start reading from {path_bam}")
    bam_files = glob.glob(f"{path_bam}/*/*.sort.bam")
    return bam_files


def bam2cov(each_bam, cov_dir, log_dir):
    sample_name = f'{os.path.basename(each_bam).split(".")[0]}_{os.path.basename(each_bam).split(".")[1]}'
    print(f"Start parsing {sample_name}")
    this_log_dir = f"{log_dir}/{sample_name}.bismark_extractor.log"
    this_sample_out = f"{cov_dir}/{sample_name}"
    if not os.path.exists(this_sample_out):
        os.makedirs(this_sample_out)
    cmd = f"bismark_methylation_extractor -s {each_bam} --bedGraph --counts -o {this_sample_out} > {this_log_dir} 2>&1"
    try:
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)


def multi_process_run(bam_files, num_process, cov_dir, log_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    for each_bam in bam_files:
        pool.apply_async(bam2cov, (each_bam, cov_dir, log_dir,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to map raw fastq and remove duplicates")
    parser.add_argument("-i", action="store", dest="mapping_dir")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-c", action="store", dest="cov_dir", help="path to put mapping result in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    results = parser.parse_args()
    if not all([results.mapping_dir, results.num_process, results.cov_dir, results.log_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    num_process = results.num_process
    cov_dir = results.cov_dir if not results.cov_dir.endswith("/") else results.cov_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    for each_dir in [cov_dir, log_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    bam_files = get_bam_files(mapping_dir)
    multi_process_run(bam_files, num_process, cov_dir, log_dir)


if __name__ == '__main__':
    main()
