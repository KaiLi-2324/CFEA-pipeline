import os
import re
import sys
import subprocess
from itertools import chain
from collections import defaultdict


"""
@Author: likai
@Email: likai@wibe.ac.cn
this script is written to merge different runs from same samples
"""


def get_run_files(path_meta, pmid):
    """
    
    :param path_meta:
    :param pmid:
    :return:
    """
    print(f"Start reading from {path_meta}")
    sample_run = dict()
    with open(path_meta) as f:
        for line in f:
            line_list = line.strip().split("\t")
            if line_list[2] == pmid:
                urls = line_list[7]
                urls_list = list(map(lambda k: re.sub(".sra", "", os.path.basename(k)), urls.strip().split(",")))
                sample_run[line_list[0]] = urls_list
    return sample_run


def check_run_file(sample_run, path_fastq):
    """
    collect run information for each sample
    :param sample_run: run information for each sample
    :param path_fastq: directory to each fastq file
    :return:
    """
    sample_run_out = defaultdict(lambda: defaultdict(list))
    print("Start checking run files")
    for each_sample in sample_run:
        print(f" Start parsing {each_sample}")
        this_sample_runs = list(map(lambda k: [f"{path_fastq}/{k}.sra.fastq",
                                               f"{path_fastq}/{k}.sra_1.fastq",
                                               f"{path_fastq}/{k}.sra_2.fastq"], sample_run[each_sample]))
        this_sample_runs = list(chain(*this_sample_runs))
        this_sample_fastq = [each_fastq for each_fastq in this_sample_runs if os.path.exists(each_fastq)]
        if not this_sample_fastq:
            raise ValueError(f"{each_sample} has no fastq files!")
        for each_fastq in this_sample_fastq:
            if each_fastq.endswith(".sra.fastq"):
                sample_run_out[each_sample]["read"].append(each_fastq)
            elif each_fastq.endswith(".sra_1.fastq"):
                sample_run_out[each_sample]["read1"].append(each_fastq)
            elif each_fastq.endswith(".sra_2.fastq"):
                sample_run_out[each_sample]["read2"].append(each_fastq)
            else:
                raise ValueError(f"sample {each_sample} has invalid fastq file types!")
    return sample_run_out


def get_cmd_paired(sample_run_out, path_out, each_sample):
    """
    shell command to merge the paired fastq files
    :param sample_run_out:
    :param path_out:
    :param each_sample:
    :return:
    """
    this_sample_read1_files = sample_run_out[each_sample]["read1"]
    this_sample_read2_files = sample_run_out[each_sample]["read2"]
    this_sample_read1_out = f"{path_out}/{each_sample}.sra_1.fastq"
    this_sample_read2_out = f"{path_out}/{each_sample}.sra_2.fastq"
    cmd_read1 = f"cat {' '.join(this_sample_read1_files)} > {this_sample_read1_out}"
    cmd_read2 = f"cat {' '.join(this_sample_read2_files)} > {this_sample_read2_out}"
    return cmd_read1, cmd_read2


def get_cmd_single(sample_run_out, path_out, each_sample):
    """
    shell command to merge the single end fastq files
    :param sample_run_out:
    :param path_out:
    :param each_sample:
    :return:
    """
    this_sample_read_files = sample_run_out[each_sample]["read"]
    this_sample_read_out = f"{path_out}/{each_sample}.sra.fastq"
    cmd_read = f"cat {' '.join(this_sample_read_files)} > {this_sample_read_out}"
    return cmd_read


def merge_fastq(sample_run_out, path_out, log_file):
    """
    rum the merge shell command to merge different run files in each sample
    :param sample_run_out:
    :param path_out:
    :param log_file:
    :return:
    """
    try:
        log = open(log_file, "w+")
        for each_sample in sample_run_out:
            print(f"  Start merging sample {each_sample}")
            this_sample_read_types = list(sample_run_out[each_sample].keys())
            if len(this_sample_read_types) == 3:
                print(f"   sample {each_sample} has both single end and paired end files!")
                log.write(f"{each_sample}")
                read_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read"]))
                read1_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read1"]))
                read2_fastq_sizes = list(map(lambda k: int(os.path.getsize(k)), sample_run_out[each_sample]["read2"]))
                if sum(read1_fastq_sizes) + sum(read2_fastq_sizes) >= sum(read_fastq_sizes):
                    cmd_read1, cmd_read2 = get_cmd_paired(sample_run_out, path_out, each_sample)
                    print("  Using paired fastq files instead!")
                    subprocess.check_output(cmd_read1, shell=True)
                    subprocess.check_output(cmd_read2, shell=True)
                else:
                    cmd_read = get_cmd_single(sample_run_out, path_out, each_sample)
                    print("  Using single fastq files instead!")
                    subprocess.check_output(cmd_read, shell=True)
            if this_sample_read_types == ["read1", "read2"]:
                cmd_read1, cmd_read2 = get_cmd_paired(sample_run_out, path_out, each_sample)
                print(cmd_read1)
                subprocess.check_output(cmd_read1, shell=True)
                print(cmd_read2)
                subprocess.check_output(cmd_read2, shell=True)
            elif this_sample_read_types == ["read"]:
                cmd_read = get_cmd_single(sample_run_out, path_out, each_sample)
                print(cmd_read)
                subprocess.check_output(cmd_read, shell=True)
            else:
                print(f"Invalid fastq file types {each_sample}")
                log.write(f"{each_sample} invalid fastq types")
        log.close()
    except Exception as e:
        print(e)


def main():
    path_meta = sys.argv[1]  # cfDNA_samples_meta_plasma.xls
    path_fastq = sys.argv[2]  # /share/pub/lik/cfDNA_database/projects/liquid.20/fastq
    path_out = sys.argv[3]  # /share/pub/lik/cfDNA_database/projects/liquid.20/fastq_merge
    pmid = sys.argv[4]  # 28263317 PMID number
    log_file = sys.argv[5]  # /share/pub/lik/cfDNA_database/projects/liquid.20/log/invalid_fastq_sample.txt
    if not os.path.exists(path_out):
        os.makedirs(path_out)
    sample_run = get_run_files(path_meta, pmid)
    sample_run_out = check_run_file(sample_run, path_fastq)
    merge_fastq(sample_run_out, path_out, log_file)


if __name__ == '__main__':
    main()
