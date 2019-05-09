import os
import sys
import glob
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
"""


def get_wig_files(path_wigs):
    print(f"Start reading from {path_wigs}")
    wig_files = glob.glob(f"{path_wigs}/*.wig")
    return wig_files


def wig_to_bw(path_wig, path_chrom_size):
    try:
        sample_name = os.path.basename(path_wig).split(".")[0]
        bw = f"{os.path.dirname(path_wig)}/{sample_name}.bw"
        if os.path.exists(bw):
            print(f"    {bw} already exists!")
        else:
            print(f"Start parsing {path_wig}")
            cmd = f"wigToBigWig {path_wig} {path_chrom_size} {bw}"
            subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)


def multi_process_run(wig_files, path_chrom_size):
    pool = multiprocessing.Pool(processes=40)
    for each_wig in wig_files:
        pool.apply_async(wig_to_bw, (each_wig, path_chrom_size,))
    pool.close()
    pool.join()


def main():
    path_wigs = sys.argv[1]
    path_chrom_size = sys.argv[2]
    wig_files = get_wig_files(path_wigs)
    multi_process_run(wig_files, path_chrom_size)


if __name__ == '__main__':
    main()
