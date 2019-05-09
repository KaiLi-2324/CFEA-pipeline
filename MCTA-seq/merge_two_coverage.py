import os
import argparse
import multiprocessing
from collections import defaultdict


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is designed to merge two coverage files by bismark_methylation_extractor
"""


def get_files(path_met_cov):
    # we have to collect MCTA samples first
    samples = os.listdir(path_met_cov)
    samples = list(map(lambda k: k.split(".")[0], samples))
    samples = list(set(samples))
    return samples


def merge_two_cov(each_sample, path_met_cov, path_cov_out):
    # since we align each read separately and extract methylation coverage separately, we
    # have to merge them into one coverage file here
    try:
        print(f"Start parsing {each_sample}")
        read1 = f"{path_met_cov}/{each_sample}.1.sort.bismark.met.cov"
        read2 = f"{path_met_cov}/{each_sample}.2.sort.bismark.met.cov"
        this_sample_out = f"{path_cov_out}/{each_sample}.cov"
        this_sample_cov = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
        for each_read in [read1, read2]:
            with open(each_read) as f:
                for line in f:
                    line_list = line.strip().split("\t")
                    this_sample_cov[line_list[0]][line_list[1]]["met"].append(line_list[3])
                    this_sample_cov[line_list[0]][line_list[1]]["unmet"].append(line_list[4])
        print(f"Start writing to {this_sample_out}")
        with open(this_sample_out, "w+") as out:
            for each_chr in this_sample_cov:
                this_chr_value = this_sample_cov[each_chr]
                for each_pos in sorted(this_chr_value, key=lambda k: int(k)):
                    this_pos_met_reads = sum(list(map(lambda k: int(k), this_chr_value[each_pos]["met"])))
                    this_pos_unmet_reads = sum(list(map(lambda k: int(k), this_chr_value[each_pos]["unmet"])))
                    this_pos_met_value = this_pos_met_reads / (this_pos_met_reads + this_pos_unmet_reads)
                    out.write(each_chr + "\t" + each_pos + "\t" +
                              each_pos + "\t" + str(this_pos_met_value) + "\t" +
                              str(this_pos_met_reads) + "\t" + str(this_pos_unmet_reads) + "\n")
    except Exception as e:
        print(e)


def multi_process_run(samples, path_met_cov, path_cov_out):
    # to run faster, we create a multiprocessing pool here with 20 processes, you can modify this number
    # according to your data
    pool = multiprocessing.Pool(processes=20)
    for each_sample in samples:
        pool.apply_async(merge_two_cov, (each_sample, path_met_cov, path_cov_out,))
    pool.close()
    pool.join()


def main():
    parser = argparse.ArgumentParser(description="A program to extract CPG from bismark coverage files")
    parser.add_argument("-i", action="store", dest="path_met_cov", help="directory to met coverage files")
    parser.add_argument("-o", action="store", dest="path_cov_out", help="path to put merged coverage files")
    results = parser.parse_args()
    if not all([results.path_met_cov, results.path_cov_out]):
        print("too few arguments, type -h for more information")
        exit(-1)
    path_met_cov = results.path_met_cov
    path_cov_out = results.path_cov_out
    if not os.path.exists(path_cov_out):
        os.makedirs(path_cov_out)
    samples = get_files(path_met_cov)
    multi_process_run(samples, path_met_cov, path_cov_out)


if __name__ == '__main__':
    main()
