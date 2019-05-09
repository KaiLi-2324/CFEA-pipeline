import os
import glob
import argparse
import subprocess
import multiprocessing


"""
@Author: likai
@Email: likai@wibe.ac.cn
"""


def get_peak_files(path_peak):
    print(f"Start reading from {path_peak}")
    peak_files = glob.glob(f"{path_peak}/*/*/*.xls")
    return peak_files


def peak_to_bed(peak, path_bed):
    print(f"Converting to bed:  {peak}")
    try:
        sample_name = os.path.basename(peak).split("_")[-1].split(".")[0]
        this_sample_out = f"{path_bed}/{sample_name}.bed"
        cmd = f"sed -n '2,$p' {peak} | grep -P '^chr(\d+|M|X|y)\s+' | cut -f 1,2,3 > {this_sample_out}"
        subprocess.check_output(cmd, shell=True)
    except Exception as e:
        print(e)


def multi_process_run(peak_files, num_process, path_bed):
    pool = multiprocessing.Pool(processes=int(num_process))
    for peak in peak_files:
        pool.apply_async(peak_to_bed, (peak, path_bed,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to convert nucleosome positioning peaks into bed format")
    parser.add_argument("-i", action="store", dest="peak_dir", help="directory to peak files")
    parser.add_argument("-p", action="store", dest="num_process", help="processes to use in the program", default=10)
    parser.add_argument("-b", action="store", dest="bed_dir", help="path to put bed files in")
    results = parser.parse_args()
    if not all([results.peak_dir, results.num_process, results.bed_dir]):
        print("too few arguments, type -h for more information")
        exit(-1)
    peak_dir = results.peak_dir
    num_process = results.num_process
    bed_dir = results.bed_dir
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    peak_files = get_peak_files(peak_dir)
    multi_process_run(peak_files, num_process, bed_dir)


if __name__ == '__main__':
    main()
