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


def complement(s):
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
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
                match_dict = {k: v for (k, v) in cache_dict.iteritems() if re.match(hc_start_block, k) is not None}
            except AttributeError:
                match_dict = {k: v for (k, v) in cache_dict.items() if re.match(hc_start_block, k) is not None}

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
    database = os.path.dirname(os.path.realpath(sys.argv[0])) + "/db/"

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


    if results.VERBOSE: print('Using Bowtie2 to find initial seeds using {} threads:'.format(results.num_threads))

    bowtie_base_cmd = "{}/bowtie2 --very-sensitive --quiet \
        --threads {} --no-hd".format(results.bowtie, results.num_threads)

    awk_cmd = "awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'"
    bowtie_options = {'hv': '--norc',
                      'hc': '--nofw',
                      'lv': '--norc',
                      'lc': '--nofw'}

    if single == 1:
        seq_options = "-U ".format(results.FASTQ)
    elif paired == 1:
        seq_options = "-U {},{}".format(results.LEFT, results.RIGHT)

    for chain_type in ['hv', 'hc', 'lv', 'lc']:

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
            proc = subprocess.check_call(cmd, shell=True)
        except:
            print('Error Bowtie2 not installed or not found in path')
            print('Program terminated.')
            exit(0)

    ptrn = output_location + "/" + results.name + ".hv"
    with open(ptrn, 'r') as f:
        line = f.readline()
    ax, bx, cx, dx = line.split('\t')
    read_length = len(ax)
    max_read_length = max(max_read_length,read_length)
    ptrn = "\t*\t*\t"
    d = defaultdict(int)

    # anchor start of heavy chain
    hv = defaultdict(int)
    if results.VERBOSE: print("Searching for an anchor in heavy chain variable region ...")
    cmd = output_location + "/" + results.name + ".hv"
    with open(cmd) as f:
        for line in f:
            ax, bx, cx, dx = line.split('\t')
            read_length = len(ax)
            max_read_length = max(max_read_length,read_length)
            if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
                continue
            if read_length < 16:
                continue
            hv[ax] += 1
            hv[revcom(ax)] += 1
            if ptrn not in line:
                d[ax] += 1
    try:
        hc_start = max(d, key=d.get)
        hc_start = hc_start.strip('N')
    except (ValueError, TypeError):
        print(('Error: ' + results.name + ' did not map to heavy chain variable region.'))
        print('Program terminated.')
        exit(0)
    if results.VERBOSE: print(("anchor: " + hc_start))
    d.clear()

    # anchor end of heavy chain
    hc = defaultdict(int)
    if results.VERBOSE: print("Searching for an anchor in heavy chain constant region ...")
    cmd = output_location + "/" + results.name + ".hc"
    with open(cmd) as f:
        for line in f:
            ax, bx, cx, dx = line.split('\t')
            read_length = len(ax)
            max_read_length = max(max_read_length,read_length)
            if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
                continue
            if read_length < 16:
                continue
            hc[ax] += 1
            hc[revcom(ax)] += 1
            if ptrn not in line:
                d[ax] += 1
    try:
        hc_end = max(d, key=d.get)
        hc_end = hc_end.strip('N')
    except (ValueError, TypeError):
        print(('Error: ' + results.name + ' did not map to heavy chain constant region.'))
        print('Program terminated.')
        exit(0)
    if results.VERBOSE: print(("anchor: " + revcom(hc_end)))
    d.clear()

    # anchor start of light chain
    lv = defaultdict(int)
    if results.VERBOSE: print("Searching for an anchor in light chain variable region ...")
    cmd = output_location + "/" + results.name + ".lv"
    with open(cmd) as f:
        for line in f:
            ax, bx, cx, dx = line.split('\t')
            read_length = len(ax)
            max_read_length = max(max_read_length,read_length)
            if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
                continue
            if read_length < 16:
                continue
            lv[ax] += 1
            lv[revcom(ax)] += 1
            if ptrn not in line:
                d[ax] += 1
    try:
        lc_start = max(d, key=d.get)
        lc_start = lc_start.strip('N')
    except (ValueError, TypeError):
        print(('Error: ' + results.name + ' did not map to light chain variable region.'))
        print('Program terminated.')
        exit(0)
    if results.VERBOSE: print(("anchor: " + lc_start))
    d.clear()

    # anchor end of light chain
    lc = defaultdict(int)
    if results.VERBOSE: print("Searching for an anchor in light chain constant region ...")
    cmd = output_location + "/" + results.name + ".lc"
    with open(cmd) as f:
        for line in f:
            ax, bx, cx, dx = line.split('\t')
            read_length = len(ax)
            max_read_length = max(max_read_length,read_length)
            if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
                continue
            if read_length < 16:
                continue
            lc[ax] += 1
            lc[revcom(ax)] += 1
            if ptrn not in line:
                d[ax] += 1
    try:
        lc_end = max(d, key=d.get)
        lc_end = lc_end.strip('N')
    except (ValueError, TypeError):
        print(('Error: ' + results.name + ' did not map to light chain constant region.'))
        print('Program terminated.')
        exit(0)
    if results.VERBOSE: print(("anchor: " + revcom(lc_end)))
    d.clear()
    ### End of single end sequencing processing

    p1 = Process(target=stitch,
                 args=(hc_end, hc_start, hc, hv, results.VERBOSE, max_read_length, results.name, output_location, 'heavy'))
    p2 = Process(target=stitch,
                 args=(lc_end, lc_start, lc, lv, results.VERBOSE, max_read_length, results.name, output_location, 'light'))
    p1.start()
    p2.start()
    p1.join()
    p2.join()

    ptrn = output_location + "/" + results.name + ".hv"
    os.remove(ptrn)
    ptrn = output_location + "/" + results.name + ".hc"
    os.remove(ptrn)
    ptrn = output_location + "/" + results.name + ".lv"
    os.remove(ptrn)
    ptrn = output_location + "/" + results.name + ".lc"
    os.remove(ptrn)


if __name__ == '__main__':
    start_time = time.time()
    main()

    print(("Done. Total run time: %s seconds" % (time.time() - start_time)))
