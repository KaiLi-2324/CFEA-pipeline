import os
import sys
import glob
import subprocess
import multiprocessing


"""
this script is desgned to convert 5hmc bam files into wiggle format
if you want to run this script, you have to install MEDIPS R package
in your R executable
> source("https://bioconductor.org/biocLite.R")
> biocLite("MEDIPS")
we've provide a R script along with this python script, and you have to put 
these two files together in the same directory
"""


def get_bam_files(path_bam):
    # get all the deduplicated bam files
    print(f"Start reading from {path_bam}")
    bam_files = glob.glob(f"{path_bam}/*.dedup.bam")
    if not bam_files:
        bam_files = glob.glob(f"{path_bam}/*.sort.bam")
    return bam_files


def bam_to_wig(bam, wig_dir, log_dir):
    # we run run_medips.R through python to convert bam into wig
    print(f"Start parsing {bam}")
    sample = os.path.basename(bam).split(".")[0]
    wig = f"{wig_dir}/{sample}.wig"
    log_file = f"{log_dir}/{sample}.medips.log"
    error_log = f"{log_dir}/{sample}.medips.log.error"
    cmd = f"Rscript run_medips.R {bam} {wig} {sample} > {log_file} 2>&1"  # the R command
    try:
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)
        with open(error_log, "w+") as log:
            log.write(str(e) + "\n")


def multi_process_run(bam_files, wig_dir, log_dir):
    pool = multiprocessing.Pool(processes=40)
    for each_bam in bam_files:
        pool.apply_async(bam_to_wig, (each_bam, wig_dir, log_dir,))
    pool.close()
    pool.join()


def main():
    if len(sys.argv) != 4:
        print("Usage: python bam2wig.py [bam dir] [wig dir] [log dir]")
        exit(-1)
    path_bam = sys.argv[1]
    wig_dir = sys.argv[2]
    log_dir = sys.argv[3]
    if not os.path.exists(wig_dir):
        os.makedirs(wig_dir)
    if not os.path.exists("run_medips.R"):
        raise FileNotFoundError("run_medips.R not in the current directory!")
    bam_files = get_bam_files(path_bam)
    multi_process_run(bam_files, wig_dir, log_dir)


if __name__ == '__main__':
    main()
