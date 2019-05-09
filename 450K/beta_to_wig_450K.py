import re
import sys
import glob
import multiprocessing
from collections import defaultdict


"""
this script is designed to convert 450K beta values into wig file
__author__ = likai
"""


def get_files(path_files):
    print(f"Start reading from {path_files}")
    files = glob.glob(f"{path_files}/*.beta_value.txt")
    return files


def parse_ref(path_ref):
    print(f"Start reading from {path_ref}")
    cg_pos = dict()
    with open(path_ref) as f:
        for line in f:
            line_list = line.strip().split("\t")
            cg_pos[line_list[0]] = ([line_list[1]], line_list[2])
    return cg_pos


def beta2wig(each_beta, cg_pos):
    print(f"Start parsing {each_beta}")
    beta_chr_pos = defaultdict(dict)
    beta_out = re.sub(".txt", ".wig", each_beta)
    try:
        with open(each_beta) as f:
            for line in f:
                if line.startswith("pos"):
                    continue
                line_list = line.strip().split("\t")
                if line_list[0] in cg_pos:
                    this_cg_values = cg_pos[line_list[0]]
                    if not this_cg_values[0].startswith("chr"):
                        beta_chr_pos[f"chr{this_cg_values[0]}"][this_cg_values[1]] = line_list[1]
                    else:
                        beta_chr_pos[this_cg_values[0]][this_cg_values[1]] = line_list[1]
                else:
                    pass
                    # print(f"no pos information found for {line_list[0]}")
        print(f"Start writing to {beta_out}")
        with open(beta_out, "w+") as out:
            for each_chr in beta_chr_pos:
                out.write(f"variableStep chrom={each_chr}\n")
                this_chr_values = beta_chr_pos[each_chr]
                for each_pos in sorted(this_chr_values, key=lambda k: int(k)):
                    out.write(each_pos + "\t" + this_chr_values[each_pos] + "\n")
    except Exception as e:
        print(e)


def multi_process_run(files, cg_pos):
    pool = multiprocessing.Pool(processes=50)
    for each_file in files:
        pool.apply_async(beta2wig, (each_file, cg_pos,))
    pool.close()
    pool.join()


def main():
    path_files = sys.argv[1]
    path_ref = sys.argv[2]
    files = get_files(path_files)
    cg_pos = parse_ref(path_ref)
    multi_process_run(files, cg_pos)


if __name__ == '__main__':
    main()
