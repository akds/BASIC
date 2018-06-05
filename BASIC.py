from multiprocessing import Process
import os.path
import itertools
import re
import sys
import time
import subprocess
import math
from collections import defaultdict
import argparse


def count_nt(seq, vec, w):
    basemap = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    letters = list(seq)
    x = [basemap[base] for base in letters if base != 'N']
    y = [base[0] for base in enumerate(letters) if base[1] != 'N']
    vec[(x,y)]+=w
    return vec


def index2seq(vec):
    basemap = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    x = list(vec)
    x = [basemap[base] for base in x]
    return ''.join(x)


basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
def complement(s):
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


def revcom(s):
    return complement(s[::-1])


def splitread(s,rl):
    table = []
    rl = min(rl,len(s))
    for x in range(rl - 1, 14, -1):
        table.append(s[-x:])
    return table


def extend(s, d, verb, rl, reverse_comp=False):
    while 1:
        table = splitread(s, rl)
        edges = [sys.maxsize]
        best_c = 0
        best_k = ''
        try:
            cache_dict = {k: v for (k, v) in d.iteritems() if table[-1] in k}
        except AttributeError:
            cache_dict = {k: v for (k, v) in d.items() if table[-1] in k}

        for hc_start_block in table:
            try:
                match_dict = {k: v for (k, v) in cache_dict.iteritems() if k.startswith(hc_start_block)}
            except AttributeError:
                match_dict = {k: v for (k, v) in cache_dict.items() if k.startswith(hc_start_block)}

            if len(match_dict) == 0:
                continue
            n = float(sum(match_dict.values()) + .5)
            v = list(match_dict.values())
            p = [(float(x) / n) for x in v]
            entropy = sum([(-1 * x * math.log(x, 2)) for x in p])
            k = list(match_dict.keys())
            if (entropy <= min(edges)) and (max(v) > 1):
                best_k = k[v.index(max(v))]
                best_rl = len(best_k)
                best_c = -(best_rl - len(hc_start_block))
                edges.append(entropy)
            else:
                continue
        if len(edges) == 1:
            return s
        if best_k.strip('N') in s:
            return s
        s += best_k[best_c:]
        s = s.strip('N')
        if verb:
            if reverse_comp:
                print(revcom(s))
            else:
                print(s)
    return s


def stitch(chain_end, chain_start, const, variable, verb, read_length, cellid, output_location, chain_name):
    if verb: print("Stitching {} chain sequence (5' <--- 3') ...".format(chain_name))
    tmp = revcom(extend(chain_end, const, verb, read_length, reverse_comp=True))
    if verb: print("Stitching {} chain sequence (5' ---> 3') ...".format(chain_name))
    const_token = extend(tmp, variable, verb, read_length)

    if (chain_start in const_token) or (revcom(chain_start) in const_token):
        if verb: print("Path found!")
        result = ">cell_id={};{}_chain;BASIC\n{}\n".format(cellid, chain_name, const_token)
    else:
        if verb: print("Re-stitching {} chain sequence (5' ---> 3') ...".format(chain_name))
        tmp = extend(chain_start, variable, verb, read_length)
        if verb: print("Re-stitching {} chain sequence (5' <--- 3') ...".format(chain_name))
        var_token = revcom(extend(revcom(tmp), const, verb, read_length, reverse_comp=True))
        if (revcom(chain_end) in var_token) or (chain_end in var_token):
            if verb: print("Path found!")
            result = ">cell_id={};{}_chain;BASIC\n{}\n".format(cellid, chain_name, var_token)
        else:
            if verb: print("Path not found!")
            result = ">cell_id={};{}_chain[variable_region_contig];BASIC\n{}\n".format(cellid, chain_name, var_token)
            result +=">cell_id={};{}_chain[constant_region_contig];BASIC\n{}\n".format(cellid, chain_name, const_token)

    cmd = output_location + "/" + cellid + ".fasta"
    with open(cmd, "a") as text_file:
        text_file.write(result)


def find_anchor(mapped_output):
    """ find the sequence anchor within processed bowtie2 output """

    max_read_length = 0

    unmapped_ptrn = "\t*\t*\t"

    all_reads = defaultdict(int)
    mapped_reads = defaultdict(int)
    with open(mapped_output) as f:
        for line in f:
            ax, bx, cx, dx = line.split('\t')
            read_length = len(ax)
            max_read_length = max(max_read_length, read_length)
            if max([float(ax.count(base))/read_length for base in ['A', 'T', 'G', 'C']]) >= .5:
                continue
            if read_length < 16:
                continue
            all_reads[ax] += 1
            all_reads[revcom(ax)] += 1
            if unmapped_ptrn not in line:
                mapped_reads[ax] += 1
    try:
        anchor = max(mapped_reads, key=mapped_reads.get)
        anchor = anchor.strip('N')
    except (ValueError, TypeError):
        print('Error: reads did not map to {}'.format(mapped_output.split('.')[-1]))
        print('Program terminated.')
        exit(0)

    return anchor, all_reads, max_read_length


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', action='store', dest='type',
                        default='BCR',
                        help='BCR or TCR (default: BCR)')

    parser.add_argument('-p', action='store', dest='num_threads',
                        default='2',
                        help='Launch p > 2 threads that will run on separate processors/cores (default: 2)')

    parser.add_argument('-n', action='store', dest='name',
                        default='result',
                        help='Name of output file (default: result)')

    parser.add_argument('-SE', action='store', type=str, dest='FASTQ',
                        default='',
                        help='Single end FASTQ file. (example: se.fastq)')

    parser.add_argument('-PE_1', action='store', type=str, dest='LEFT',
                        default='',
                        help='Paired end (left) FASTQ file. -PE_2 is required and pairs must match order. (example: pe_1.fastq)')

    parser.add_argument('-PE_2', action='store', type=str, dest='RIGHT',
                        default='',
                        help='Paired end (right) FASTQ files. (example: pe_2.fastq)')

    parser.add_argument('-g', action='store', dest='genome',
                        default='human',
                        help='human or mouse (default: human)')

    parser.add_argument('-b', action='store', dest='bowtie',
                        default='./',
                        help='Absolute path to directory that contains the bowtie2 executable')

    parser.add_argument('-o', action='store', dest='output_location',
                        default='./',
                        help='Output dir (default: none -- current working directory)')

    parser.add_argument('-v', action='store_true', dest='VERBOSE',
                        default=False,
                        help='Turns on verbosity (more details)')

    parser.add_argument('--version', action='version', version='%(prog)s 1.2 (beta)')

    return parser.parse_args()


def main():
    single = 0
    paired = 0

    max_read_length = 0;

    results = parse_args()

    output_location = "./" + str(re.sub('\W+', '', results.output_location))
    output_file = str(re.sub(r'\W+', '', results.name))
    database = os.path.dirname(os.path.realpath(sys.argv[0])) + "/db"

    if results.VERBOSE: print(('Run ID: ' + results.name ))

    if (len(results.LEFT) == 0) | (len(results.RIGHT) == 0):
        if len(results.FASTQ) == 0:
            print('Sequencing data missing')
            print('Program terminated.')
            exit(0)
        else:
            if results.VERBOSE: print(('List of single end (SE) sequencing inputs: ' + results.FASTQ + ' ...'))
            single = 1
    elif (len(results.LEFT) != 0) & (len(results.RIGHT) != 0):
        if results.VERBOSE: print(('List of paired end (PE_1) sequencing inputs (left): ' + results.LEFT + ' ...'))
        if results.VERBOSE: print(('List of paired end (PE_2) sequencing inputs (right): ' + results.RIGHT + ' ...'))
        paired = 1
    else:
        print('Both data left and right paired end sequencing must be included')
        print('Program terminated.')
        exit(0)

    if single == 1:
        f_list = results.FASTQ.split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if v == False:
                print((s + ' not found'))
                print('Program terminated.')
                exit(0)
    elif paired == 1:
        f_list = results.LEFT.split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if v == False:
                print((s + ' not found'))
                print('Program terminated.')
                exit(0)

        f_list = results.RIGHT.split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if v == False:
                print((s + ' not found'))
                print('Program terminated.')
                exit(0)
    else:
        print('Sequencing data missing')
        print('Program terminated.')
        exit(0)


    try:
        os.makedirs(output_location)
    except OSError:
        if not os.path.isdir(output_location):
            print(('Output directory: ' + output_location + ' is not writeable -- check permissions or try different directory'))
            print('Program terminated.')
            exit(0)

    if results.type not in ['BCR', 'TCR']:
        print('Currently only BCR or TCR assembly is supported')
        print('Program terminated.')
        exit(0)
    else:
        if results.VERBOSE: print(('Assembling ' + results.type + ' sequences ...'))

    valid_genomes = ['human', 'mouse']
    if not (results.genome in valid_genomes):
        print('Unsupported genome')
        print('Program terminated.')
        exit(0)
    else:
        db_path = '{}/{}/{}'.format(database, results.genome, results.type)
        check_read = os.access('{}/hc'.format(db_path), os.R_OK)
        if check_read:
            if results.VERBOSE: print(('Using ' +  db_path + ' genome files ...'))
        else:
            print('Genome index files not found: {}'.format(db_path))
            print('Check database path (-d) and that index files are readable')
            print('Program terminated.')
            exit(0)


    # Bowtie2 for mapping reads
    if results.VERBOSE: print('Using Bowtie2 to find initial seeds using {} threads:'.format(results.num_threads))

    bowtie_base_cmd = "{}/bowtie2 --very-sensitive --quiet \
        --threads {} --no-hd".format(results.bowtie, results.num_threads)

    # awk extracts read sequence (column 10), name of ref sequence aligned to (3)
    # CIGAR representation of alignment (6), and 1-based offset into the forward
    # reference strand where leftmost character of the alignment occurs (4)
    awk_cmd = "awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'"
    bowtie_options = {'hv': '--norc',
                      'hc': '--nofw',
                      'lv': '--norc',
                      'lc': '--nofw'}

    if single == 1:
        seq_options = "-U ".format(results.FASTQ)
    elif paired == 1:
        seq_options = "-U {},{}".format(results.LEFT, results.RIGHT)

    # chain types- heavy variable (hv), heavy constant (hc) and likewise for light
    chains = ['hv', 'hc', 'lv', 'lc']
    for chain_type in chains:

        output_path = '{}/{}.{}'.format(output_location, results.name, chain_type)
        if results.VERBOSE:
            print('Mapping reads to {}...'.format(chain_type))

        try:
            cmd = "{} {} -x {}/{} {} | {} > {} ".format(bowtie_base_cmd,
                                                        bowtie_options[chain_type],
                                                        db_path,
                                                        chain_type,
                                                        seq_options,
                                                        awk_cmd,
                                                        output_path)

            if results.VERBOSE: print("Calling bowtie with command:\n{}".format(cmd))
            proc = subprocess.check_call(cmd, shell=True)
        except:
            print('Error Bowtie2 not installed or not found in path')
            print('Program terminated.')
            exit(0)

    # find anchors and collapse reads into dictionary
    anchors_str = {}
    anchors_dict = {}
    for chain_type in chains:
        if results.VERBOSE: print("Searching for an anchor in {} ...".format(chain_type))

        output_path = "{}/{}.{}".format(output_location, results.name, chain_type)
        anchors_str[chain_type], anchors_dict[chain_type], max_rl = find_anchor(output_path)
        max_read_length = max(max_read_length, max_rl)

        if results.VERBOSE:
            if chain_type.endswith('v'):
                print("anchor: {}".format(anchors_str[chain_type]))
            else:
                print("anchor: {}".format(revcom(anchors_str[chain_type])))

    # stitch ends of each chain together
    p1 = Process(target=stitch,
                 args=(anchors_str['hc'], anchors_str['hv'], anchors_dict['hc'], anchors_dict['hv'],
                       results.VERBOSE, max_read_length, results.name, output_location, 'heavy'))
    p2 = Process(target=stitch,
                 args=(anchors_str['lc'], anchors_str['lv'], anchors_dict['lc'], anchors_dict['lv'],
                       results.VERBOSE, max_read_length, results.name, output_location, 'light'))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    # remove intermediate bowtie outputs
    for chain_type in chains:
        output_path = "{}/{}.{}".format(output_location, results.name, chain_type)
        os.remove(output_path)


if __name__ == '__main__':
    start_time = time.time()
    main()

    print(("Done. Total run time: %s seconds" % (time.time() - start_time)))
