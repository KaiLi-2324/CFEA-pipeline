import re
import sys
import glob
import gzip
import multiprocessing
from collections import defaultdict


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is designed to convert WGBS/RRBS coverage files into wig file
"""


def get_files(path_cov):
    print(f"Start reading from {path_cov}")
    cov_files = glob.glob(f"{path_cov}/*.cov.gz")
    if not cov_files:
        cov_files = glob.glob(f"{cov_files}/*.cov")
    if not cov_files:
        raise ValueError("Empty cov file dir!")
    return cov_files


def cov2wig(each_cov):
    print(f"Start parsing {each_cov}")
    total_sample_cov = defaultdict(dict)
    try:
        if each_cov.endswith(".gz"):
            f = gzip.open(each_cov)
            cov_out = re.sub(".cov.gz", ".wig.gz", each_cov)
        else:
            f = open(each_cov)
            cov_out = re.sub(".cov", ".wig.gz", each_cov)
        for line in f:
            if isinstance(line, bytes):
                line = line.decode()
            line_list = line.strip().split("\t")
            total_sample_cov[line_list[0]][line_list[1]] = line_list[3]
        print(f"Start writing to {cov_out}")
        # to save space, we output gzip format
        with gzip.open(cov_out, "wb") as out:
            for each_chr in total_sample_cov:
                out.write(f"variableStep chrom={each_chr}\n".encode())
                this_chr_values = total_sample_cov[each_chr]
                for each_pos in sorted(this_chr_values, key=lambda k: int(k)):
                    out.write(each_pos.encode() + "\t".encode() + this_chr_values[each_pos].encode() + "\n".encode())
        print(f"Finishing writing to {cov_out}")
    except Exception as e:
        print(e)


def multi_process_run(cov_files):
    pool = multiprocessing.Pool(processes=50)
    for each_file in cov_files:
        pool.apply_async(cov2wig, (each_file,))
    pool.close()
    pool.join()


def main():
    path = sys.argv[1]  # path to cov files
    cov_files = get_files(path)
    multi_process_run(cov_files)


if __name__ == '__main__':
    main()
