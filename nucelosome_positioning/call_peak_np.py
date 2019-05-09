import os
import glob
import argparse
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is written to call peaks from the deduplicated bam files with danpos software
"""


def get_bam_files(path_bam):
    # get all the bam files
    print(f"Start reading from {path_bam}")
    bam_files = glob.glob(f"{path_bam}/*.dedup.bam")
    return bam_files


def call_peak(bam, peak_dir, log_dir):
    # call peak from deduplicated bam files with danpos
    # you have to add danpos into your local environment
    prefix = os.path.basename(bam).split(".")[0]
    this_sample_peak_out = f"{peak_dir}/{prefix}"
    this_sample_log = f"{log_dir}/{prefix}.danpos.log"
    this_sample_error_log = f"{log_dir}/{prefix}.danpos.log.error"
    print(f" Calling peak:  {prefix}")
    try:
        cmd = f"danpos.py dpos {bam} -o {this_sample_peak_out} > {this_sample_log} 2>&1"
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)
        with open(this_sample_error_log, "W+") as out:
            out.write(str(e) + "\n")


def multi_process_run(bam_files, num_process, peak_dir, log_dir):
    pool = multiprocessing.Pool(processes=int(num_process))
    for bam in bam_files:
        pool.apply_async(call_peak, (bam, peak_dir, log_dir,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to call peaks from nucleosome positioning bam files")
    parser.add_argument("-i", action="store", dest="mapping_dir", help="directory to bam files")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-m", action="store", dest="peak_dir", help="path to put called peaks in")
    parser.add_argument("-l", action="store", dest="log_dir", help="path to put log files")
    results = parser.parse_args()
    if not all([results.mapping_dir, results.num_process, results.peak_dir, results.log_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    mapping_dir = results.mapping_dir if not results.mapping_dir.endswith("/") else results.mapping_dir.rstrip("/")
    num_process = results.num_process
    peak_dir = results.peak_dir if not results.peak_dir.endswith("/") else results.peak_dir.rstrip("/")
    log_dir = results.log_dir if not results.log_dir.endswith("/") else results.log_dir.rstrip("/")
    for each_dir in [peak_dir, log_dir]:
        if not os.path.exists(each_dir):
            os.makedirs(each_dir)
    bam_files = get_bam_files(mapping_dir)
    multi_process_run(bam_files, num_process, peak_dir, log_dir)


if __name__ == '__main__':
    main()
