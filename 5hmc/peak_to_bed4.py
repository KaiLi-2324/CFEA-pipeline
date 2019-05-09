import re
import os
import glob
import argparse
import multiprocessing


def get_peak_files(path_peaks):
    print(f"Start parsing {path_peaks}")
    peak_files = glob.glob(f"{path_peaks}/*/*_peaks.xls")
    if not peak_files:
        peak_files = glob.glob(f"{path_peaks}/*_peaks.xls")
    return peak_files


def peak_to_bed4(peak_file, bed_dir):
    print(f"Start parsing {peak_file}")
    try:
        sample = os.path.basename(peak_file).split("_")[0]
        this_sample_bed = f"{bed_dir}/{sample}.bed"
        out_bed = open(this_sample_bed, "w+")
        with open(peak_file) as file:
            for line in file:
                if re.search(r"^#", line) or re.search(r"^chr\s+", line) or re.search(r"^\n$", line):
                    continue
                line_list = line.strip().split("\t")
                out_bed.write(line_list[0] + "\t" + line_list[1] + "\t"
                              + line_list[2] + "\t" + line_list[9] + "\n")
        out_bed.close()
    except Exception as e:
        print(e)


def multi_process_run(peak_files, num_process, bed_dir):
    pool = multiprocessing.Pool(processes=num_process)
    for each_peak in peak_files:
        pool.apply_async(peak_to_bed4, (each_peak, bed_dir,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to convert 5hmc peak files into bed4 format")
    parser.add_argument("-i", action="store", dest="peak_dir", help="directory to your peak files")
    parser.add_argument("-b", action="store", dest="bed_dir", help="directory to put bed4 files")
    parser.add_argument("-p", action="store", dest="num_process", help="number of processed to use in the program")
    results = parser.parse_args()
    if not all([results.peak_dir, results.bed_dir, results.num_process]):
        print("too few arguments, type -h for more information")
        exit(-1)
    peak_dir = results.peak_dir
    bed_dir = results.bed_dir
    num_process = int(results.num_process)
    if not os.path.exists(bed_dir):
        os.makedirs(bed_dir)
    peak_files = get_peak_files(peak_dir)
    multi_process_run(peak_files, num_process, bed_dir)


if __name__ == '__main__':
    main()
