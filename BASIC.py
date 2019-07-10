#!/usr/bin/env python
from multiprocessing import Process, Queue
import os.path
import re
import sys
import time
import subprocess
import math
from collections import defaultdict
import argparse


basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}


def complement(s):
    letters = [basecomplement[base] for base in s]
    return ''.join(letters)


def revcom(s):
    return complement(s[::-1])


def splitread(s, rl):
    table = []
    rl = min(rl, len(s))
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
            cache = {k: v for (k, v) in d.iteritems() if table[-1] in k}
        except AttributeError:
            cache = {k: v for (k, v) in d.items() if table[-1] in k}

        for hc_start_block in table:
            try:
                match_dict = {k: v for (k, v) in cache.iteritems() if k.startswith(hc_start_block)}
            except AttributeError:
                match_dict = {k: v for (k, v) in cache.items() if k.startswith(hc_start_block)}

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


def stitch(chain_end, chain_start, const, variable, verb, read_length, cellid,
           chain_name, result_queue):
    if chain_end:
        if verb:
            print("Stitching {} chain sequence "
                  "(5' <--- 3') ...".format(chain_name))
        tmp = revcom(extend(chain_end, const, verb, read_length,
                            reverse_comp=True))
        if verb:
            print("Stitching {} chain sequence "
                  "(5' ---> 3') ...".format(chain_name))
        const_token = extend(tmp, variable, verb, read_length)
    else:
        const_token = ''

    if chain_start and ((chain_start in const_token) or (revcom(chain_start) in const_token)):
        if verb:
            print("Path found for the {} chain!".format(chain_name))
        result = ">cell_id={};{}_chain;BASIC\n{}\n".format(cellid, chain_name,
                                                           const_token)
    else:
        if chain_start:
            if verb:
                print("Re-stitching {} chain sequence "
                      "(5' ---> 3') ...".format(chain_name))
            tmp = extend(chain_start, variable, verb, read_length)
            if verb:
                print("Re-stitching {} chain sequence "
                      "(5' <--- 3') ...".format(chain_name))
            var_token = revcom(extend(revcom(tmp), const, verb, read_length,
                                      reverse_comp=True))
        else:
            var_token = ''
        if chain_end and ((revcom(chain_end) in var_token) or (chain_end in var_token)):
            if verb:
                print("Path found for the {} chain!".format(chain_name))
            result = ">cell_id={};{}_chain;BASIC\n{}\n".format(cellid,
                                                               chain_name,
                                                               var_token)
        else:
            if verb:
                print("Path not found for the {} chain!".format(chain_name))
            result = ''
            if var_token:
                if verb:
                    print("However, found partial {} chain variable region "
                          "contig".format(chain_name))
                result += ">cell_id={};{}_chain[variable_region_contig];BASIC\n{}\n".format(cellid, chain_name, var_token)
            else:
                if verb:
                    print("No {} chain variable region contig "
                          "assembled".format(chain_name))
            if const_token:
                if verb:
                    print("However, found partial {} chain constant region "
                          "contig".format(chain_name))
                result += ">cell_id={};{}_chain[constant_region_contig];BASIC\n{}\n".format(cellid, chain_name, const_token)
            else:
                if verb:
                    print("No {} chain constant region contig "
                          "assembled".format(chain_name))

    result_queue.put({chain_name: result})


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
            if max([float(ax.count(base))/read_length for base in [
                    'A', 'T', 'G', 'C']]) >= .5:
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
        anchor = ''

    return anchor, all_reads, max_read_length


def write_fasta_str(str_out, output_file):
    with open(output_file, 'a') as out:
        out.write(str_out)


def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', action='store', dest='type',
                        default='BCR',
                        help='Type of receptor. Choices: "BCR" or "TCR" '
                        '(default: BCR)')

    parser.add_argument('-p', action='store', type=int, dest='num_threads',
                        default='2',
                        help='Launch p > 2 threads that will run on separate '
                        'processors/cores (default: 2)')

    parser.add_argument('-n', action='store', dest='name',
                        default='result',
                        help='Name of output file. Note: a ".fasta" extension '
                        'will be added (default: result.fasta)')

    parser.add_argument('-SE', action='store', type=str, dest='FASTQ',
                        default='',
                        help='Single end FASTQ file (optionally gzipped). '
                        '(example: se.fastq)')

    parser.add_argument('-PE_1', action='store', type=str, dest='LEFT',
                        default='',
                        help='Paired end (left) FASTQ file (optionally '
                        'gzipped). -PE_2 is required and pairs must match '
                        'order. (example: pe_1.fastq)')

    parser.add_argument('-PE_2', action='store', type=str, dest='RIGHT',
                        default='',
                        help='Paired end (right) FASTQ file (optionally '
                        'gzipped). (example: pe_2.fastq)')

    parser.add_argument('-g', action='store', dest='genome',
                        default='human',
                        help='Options: "human" or "mouse" '
                        '(default: human). Note: other species are possible '
                        'by adding the appropriate Bowtie2 indices and '
                        'following the existing db/ directory structure')

    parser.add_argument('-b', action='store', dest='bowtie',
                        default='./',
                        help='Absolute path to bowtie2 executable or '
                        'directory containing it')

    parser.add_argument('-t', action='store', dest='tmpdir',
                        default='./',
                        help='Path to directory for writing intermediate '
                        'files. (default: current working directory)')

    parser.add_argument('-o', action='store', dest='output_location',
                        default='./',
                        help='Output directory (default: current working '
                        'directory)')

    parser.add_argument('-a', action='store_true', dest='allow_partial',
                        default=False,
                        help='Allow for partial reconstruction i.e. do not '
                        'terminate if one or both chains do not assemble.')

    parser.add_argument('-v', action='store_true', dest='VERBOSE',
                        default=False,
                        help='Turns on verbosity for more details.')

    parser.add_argument('-d', action='store', type=str, dest='database_path',
                        default='',
                        help='Path to database directory. Defaults to '
                        '<path of BASIC.py>/db')

    parser.add_argument('--version', action='version',
                        version='%(prog)s 1.4.1')

    return vars(parser.parse_args())


def config_args(arg_dict):
    """ Validate and configure arguments """

    if arg_dict['VERBOSE']:
        print(('Run ID: ' + arg_dict['name']))

    # configure single vs. paired end reads
    if (len(arg_dict['LEFT']) == 0) | (len(arg_dict['RIGHT']) == 0):
        if len(arg_dict['FASTQ']) == 0:
            print('Sequencing data missing')
            print('Program terminated.')
            exit(1)
        else:
            if arg_dict['VERBOSE']:
                print('List of single end (SE) sequencing inputs: '
                      '{} ...'.format(arg_dict['FASTQ']))
            arg_dict['read_type'] = 'single'
    elif (len(arg_dict['LEFT']) != 0) & (len(arg_dict['RIGHT']) != 0):
        if arg_dict['VERBOSE']:
            print('List of paired end (PE_1) sequencing inputs '
                  '(left): {} ...'.format(arg_dict['LEFT']))
        if arg_dict['VERBOSE']:
            print('List of paired end (PE_2) sequencing inputs '
                  ' (right): {} ...'.format(arg_dict['RIGHT']))
        arg_dict['read_type'] = 'paired'
    else:
        print('Both left and right paired end sequencing must be included')
        print('Program terminated.')
        exit(1)
    if arg_dict['read_type'] == 'single':
        f_list = arg_dict['FASTQ'].split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if not v:
                print((s + ' not found'))
                print('Program terminated.')
                exit(1)
    elif arg_dict['read_type'] == 'paired':
        f_list = arg_dict['LEFT'].split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if not v:
                print((s + ' not found'))
                print('Program terminated.')
                exit(1)

        f_list = arg_dict['RIGHT'].split(',')
        for s in f_list:
            v = os.path.isfile(s)
            if not v:
                print((s + ' not found'))
                print('Program terminated.')
                exit(1)
    else:
        print('Sequencing data missing')
        print('Program terminated.')
        exit(1)

    # setup output directory
    output_location = str(re.sub('[^a-zA-Z0-9_. /-]',
                                 '',
                                 arg_dict['output_location']))
    if output_location != arg_dict['output_location']:
        print('Warning: the output directory (-o) contained unexpected '
              'characters, using the following instead: '
              '{}'.format(output_location))
    arg_dict['output_location'] = output_location
    try:
        os.makedirs(output_location)
    except OSError:
        if not os.path.isdir(output_location):
            print(('Output directory: {} is not writeable -- '
                   ' check permissions or try different '
                   'directory'.format(arg_dict['output_location'])))
            print('Program terminated.')
            exit(1)

    # validate receptor type
    if arg_dict['type'] not in ['BCR', 'TCR']:
        print('Currently only BCR or TCR assembly is supported')
        print('Program terminated.')
        exit(1)
    else:
        if arg_dict['VERBOSE']:
            print(('Assembling ' + arg_dict['type'] + ' sequences ...'))

    # validate genomes
    valid_genomes = ['human', 'mouse']
    if not (arg_dict['genome'] in valid_genomes):
        print('Unsupported genome')
        print('Program terminated.')
        exit(1)

    # validate database path
    if not arg_dict['database_path']:
        # default to <path to BASIC.py>/db
        arg_dict['database_path'] = os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "db")

    arg_dict['db_path'] = os.path.join(arg_dict['database_path'],
                                       arg_dict['genome'],
                                       arg_dict['type'])
    check_read = os.access(os.path.join(arg_dict['db_path'], 'hc'), os.R_OK)
    if check_read:
        if arg_dict['VERBOSE']:
            print(('Using ' + arg_dict['db_path'] + ' genome files ...'))
    else:
        print('Genome index files not found: {}'.format(arg_dict['db_path']))
        print('Check database path (-d) and that index files are readable')
        print('Program terminated.')
        exit(1)

    # validate bowtie2 as executable
    # accepts directory containing bowtie2 or the executable itself
    if os.path.basename(arg_dict['bowtie']) != 'bowtie2':
        if os.path.isfile(os.path.join(arg_dict['bowtie'], 'bowtie2')):
            arg_dict['bowtie'] = os.path.join(arg_dict['bowtie'], 'bowtie2')

    if os.path.isfile(arg_dict['bowtie']) and os.access(arg_dict['bowtie'],
                                                        os.X_OK):
        try:
            subprocess.check_call("{} --version".format(arg_dict['bowtie']),
                                  shell=True)
        except subprocess.CalledProcessError:
            print("Error: Bowtie2 was not able to be called correctly. "
                  "Please check the '-b' flag.")
            raise
    else:
        raise EnvironmentError("Error: Bowtie2 was not found or was not "
                               "executable at location: {}.\nPlease specify "
                               "correctly with the '-b' "
                               "flag".format(arg_dict['bowtie']))

    return arg_dict


def run_bowtie2(args):

    # Bowtie2 for mapping reads
    if args['VERBOSE']:
        print('Using Bowtie2 to find initial seeds using {} threads:'.format(
            args['num_threads']))

    bowtie_base_cmd = "{} --very-sensitive --quiet \
        --threads {} --no-hd".format(args['bowtie'], args['num_threads'])

    # awk extracts read seq (col 10), ref sequence name of alignment (col 3)
    # CIGAR representation of alignment (6), 1-based offset into the forward
    # reference strand where leftmost character of the alignment occurs (4)
    awk_cmd = "awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'"
    bowtie_options = {'hv': '--norc',
                      'hc': '--nofw',
                      'lv': '--norc',
                      'lc': '--nofw'}

    if args['read_type'] == 'single':
        seq_options = "-U {}".format(args['FASTQ'])
    elif args['read_type'] == 'paired':
        seq_options = "-U {},{}".format(args['LEFT'], args['RIGHT'])

    # chain types- heavy variable (hv), heavy constant (hc); likewise for light
    for chain_type in ['hv', 'hc', 'lv', 'lc']:
        output_path = '{}/{}.{}'.format(
            args['tmpdir'], args['name'], chain_type)
        if args['VERBOSE']:
            print('Mapping reads to {}...'.format(chain_type))

        cmd = "{} {} -x {}/{} {} | {} > {} ".format(bowtie_base_cmd,
                                                    bowtie_options[chain_type],
                                                    args['db_path'],
                                                    chain_type,
                                                    seq_options,
                                                    awk_cmd,
                                                    output_path)

        if args['VERBOSE']:
            print("Calling bowtie with command:\n{}".format(cmd))

        subprocess.check_call(cmd, shell=True)


def main():

    parsed_args = parse_args()

    args = config_args(parsed_args)

    run_bowtie2(args)

    # find anchors and collapse reads into dictionary
    anchors_str = {}
    anchors_dict = {}
    max_read_length = 0
    chains = ['hv', 'hc', 'lv', 'lc']
    for chain_type in chains:
        if args['VERBOSE']:
            print("Searching for an anchor in {} ...".format(chain_type))

        output_path = "{}/{}.{}".format(
            args['tmpdir'], args['name'], chain_type)
        anchors_str[chain_type], anchors_dict[chain_type], max_rl = find_anchor(output_path)
        max_read_length = max(max_read_length, max_rl)

        # exit if no anchor found and the allow_partial argument was not used
        if not anchors_str[chain_type] and not args['allow_partial']:
            print('Error: reads did not map to {}'.format(chain_type))
            print('Program terminated.')
            exit(1)

        if args['VERBOSE']:
            if chain_type.endswith('v'):
                print("anchor: {}".format(anchors_str[chain_type]))
            else:
                print("anchor: {}".format(revcom(anchors_str[chain_type])))

    # stitch ends of each chain together in separate processes
    if args['type'] == 'BCR':
        assembly_names = ['heavy', 'light']
    elif args['type'] == 'TCR':
        assembly_names = ['TRA', 'TRB']

    # results stored in a Queue
    q = Queue()
    p1 = Process(target=stitch,
                 args=(anchors_str['hc'], anchors_str['hv'],
                       anchors_dict['hc'], anchors_dict['hv'],
                       args['VERBOSE'], max_read_length, args['name'],
                       assembly_names[0], q))
    p2 = Process(target=stitch,
                 args=(anchors_str['lc'], anchors_str['lv'],
                       anchors_dict['lc'], anchors_dict['lv'],
                       args['VERBOSE'], max_read_length, args['name'],
                       assembly_names[1], q))
    p1.start()
    p2.start()

    # get the result of the two processes
    fasta_res = {}
    for i in [0, 1]:
        fasta_res.update(q.get())

    # wait for processes to finish
    p1.join()
    p2.join()

    combined_fasta_str = ''.join([fasta_res[x] for x in assembly_names])

    outfile = '{}/{}.fasta'.format(args['output_location'], args['name'])
    write_fasta_str(combined_fasta_str, outfile)

    # remove intermediate bowtie outputs
    for chain_type in chains:
        output_path = "{}/{}.{}".format(
            args['tmpdir'], args['name'], chain_type)
        os.remove(output_path)


if __name__ == '__main__':
    start_time = time.time()
    main()

    print(("Done. Total run time: %s seconds" % (time.time() - start_time)))
