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


def splitread(s, rl):
    table = []
    rl = min(rl,len(s))
    for x in range(rl - 1, 14, -1):
    	table.append(s[-x:])
    return table


def extend_right(s, d, verb, rl):
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
            	match_dict = {k: v for (k, v) in cache_dict.iteritems() if re.match(hc_start_block, k)}
            except AttributeError:
            	match_dict = {k: v for (k, v) in cache_dict.items() if re.match(hc_start_block, k)}

            if len(match_dict) == 0:
                continue
            n = float(sum(match_dict.values()) + .5)
            v = list(match_dict.values())
            p = [(float(x) / n) for x in v]
            entropy = sum([(-1 * x * math.log(x, 2)) for x in p])
            k = list(match_dict.keys())
            if (entropy <= min(edges)) and (max(v) > 1):
                best_k = k[v.index(max(v))]
                best_c = -(rl - len(hc_start_block))
                edges.append(entropy)
            else:
                continue
        if len(edges) == 1:
            return s
        if best_k.strip('N') in s:
            return s
        s += best_k[best_c:]
        s = s.strip('N')
        if verb: print(s)
    return s


def extend_rc_left(s, d, verb, rl):
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
            	match_dict = {k: v for (k, v) in cache_dict.iteritems() if re.match(hc_start_block, k)}
            except AttributeError:
            	match_dict = {k: v for (k, v) in cache_dict.items() if re.match(hc_start_block, k)}
            
            if len(match_dict) == 0:
                continue
            n = float(sum(match_dict.values()) + .5)
            v = list(match_dict.values())
            p = [(float(x) / n) for x in v]
            entropy = sum([(-1 * x * math.log(x, 2)) for x in p])
            k = list(match_dict.keys())
            if (entropy <= min(edges)) and (max(v) > 1):
                best_k = k[v.index(max(v))]
                best_c = -(rl - len(hc_start_block))
                edges.append(entropy)
            else:
                continue
        if len(edges) == 1:
            return s
        if best_k.strip('N') in s:
            return s
        s += best_k[best_c:]
        s = s.strip('N')
        if verb: print((revcom(s)))
    return s



def heavy(hc_end, hc_start, c_hc, v_hc, verb, read_length, cellid, output_location):
    if verb: print("Stitching heavy chain sequence (5' <--- 3') ...")
    tmp = revcom(extend_rc_left(hc_end, c_hc, verb, read_length))
    if verb: print("Stitching heavy chain sequence (5' ---> 3') ...")
    const_token = extend_right(tmp, v_hc, verb, read_length)

    if (hc_start in const_token) or (revcom(hc_start) in const_token):
        if verb: print("Path found!")
        result = ">cell_id=" + cellid + ";heavy_chain;BASIC\n" + const_token + "\n"
    else:
        if verb: print("Re-stitching heavy chain sequence (5' ---> 3') ...")
        tmp = extend_right(hc_start, v_hc, verb, read_length)
        if verb: print("Re-stitching heavy chain sequence (5' <--- 3') ...")
        var_token = revcom(extend_rc_left(revcom(tmp), c_hc, verb, read_length))
        if (revcom(hc_end) in var_token) or (hc_end in var_token):
            if verb: print("Path found!")
            result = ">cell_id=" + cellid + ";heavy_chain;BASIC\n" + var_token + "\n"
        else:
            if verb: print("Path not found!")
            result = ">cell_id=" + cellid + ";heavy_chain[variable_region_contig];BASIC\n" + var_token + "\n" + ">cell_id=" + cellid + ";heavy_chain[constant_region_contig];BASIC\n" + const_token + "\n"
    cmd = output_location + "/" + cellid + ".fasta"
    with open(cmd, "a") as text_file:
        text_file.write(result)


def light(lc_end, lc_start, c_lc, v_lc, verb, read_length, cellid, output_location):
    if verb: print("Stitching light chain sequence (5' <--- 3') ...")
    tmp = revcom(extend_rc_left(lc_end, c_lc, verb, read_length))
    if verb: print("Stitching light chain sequence (5' ---> 3') ...")
    const_token = extend_right(tmp, v_lc, verb, read_length)

    if (lc_start in const_token) or (revcom(lc_start) in const_token):
        if verb: print("Path found!")
        result = ">cell_id=" + cellid + ";light_chain;BASIC\n" + const_token + "\n"
    else:
        if verb: print("Re-stitching light chain sequence (5' ---> 3') ...")
        tmp = extend_right(lc_start, v_lc, verb, read_length)
        if verb: print("Re-stitching light chain sequence (5' <--- 3') ...")
        var_token = revcom(extend_rc_left(revcom(tmp), c_lc, verb, read_length))
        if (revcom(lc_end) in var_token) or (lc_end in var_token):
            if verb: print("Path found!")
            result = ">cell_id=" + cellid + ";light_chain;BASIC\n" + var_token + "\n"
        else:
            if verb: print("Path not found!")
            result = ">cell_id=" + cellid + ";light_chain[variable_region_contig];BASIC\n" + var_token + "\n" + ">cell_id=" + cellid + ";light_chain[constant_region_contig];BASIC\n" + const_token + "\n"
    cmd = output_location + "/" + cellid + ".fasta"
    with open(cmd, "a") as text_file:
        text_file.write(result)


parser = argparse.ArgumentParser()

parser.add_argument('-i', action='store', dest='type',
                    default='BCR',
                    help='BCR, TCR, or HLA (default: BCR)')

parser.add_argument('-p', action='store', dest='constant_value',
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
                    default='hg19',
                    help='hg19 or mm10 (default: hg19)')

parser.add_argument('-b', action='store', dest='bowtie',
                    default='./',
                    help='Absolute path to directory that contains the bowtie2 executable')

parser.add_argument('-o', action='store', dest='output_location',
                    default='./',
                    help='Output dir (default: none -- current working directory)')

parser.add_argument('-v', action='store_true', dest='VERBOSE',
                    default=False,
                    help='Turns on verbosity (more details)')

parser.add_argument('--version', action='version', version='%(prog)s 1.0.1')

single = 0
paired = 0

results = parser.parse_args()
start_time = time.time()

output_location = "./" + str(re.sub('\W+', '', results.output_location))
output_file = str(re.sub(r'\W+', '', results.name))
tmp_name = results.name
database = os.path.dirname(os.path.realpath(sys.argv[0])) + "/db/"

if results.VERBOSE: print(('Run ID: ' + tmp_name ))

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

valid_genomes = ['hg19', 'mm10']
if not (results.genome in valid_genomes):
    print('Unsupported genome')
    print('Program terminated.')
    exit(0)
else:
    db_path = database + results.genome + '_ighc'
    check_read = os.access(db_path, os.R_OK)
    if check_read:
        db_path = database + results.genome + '_*'
        if results.VERBOSE: print(('Using ' +  db_path + ' genome files ...'))
    else:
        db_path = database + results.genome + '_*'
        print(('Genome index files not found: ' + db_path))
        print('Check database path (-d) and that index files are readable')
        print('Program terminated.')
        exit(0)

if results.type != 'BCR':
    print('Currently only BCR assembly is supported')
    print('Program terminated.')
    exit(0)
else:
    if results.VERBOSE: print(('Assembling ' + results.type + ' sequences ...'))

if results.VERBOSE: print(('Using ' + results.constant_value + ' threads ...'))

if single == 1:
    if results.VERBOSE == True: print('Using Bowtie2 to find initial seeds:')
    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --norc --no-hd -x " + database + results.genome + "_ighv " + \
              "--threads " + results.constant_value + " -U " + results.FASTQ + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".ighv"
        if results.VERBOSE == True: print("Mapping reads to heavy chain variable region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive --quiet --nofw --no-hd -x " + database + results.genome + "_ighc " + \
              "--threads " + results.constant_value + " -U " + results.FASTQ + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".ighc"
        if results.VERBOSE == True: print("Mapping reads to heavy chain constant region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --norc  --no-hd -x " + database + results.genome + "_iglv " + \
              "--threads " + results.constant_value + " -U " + results.FASTQ + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}' >" + output_location + "/" + tmp_name + ".iglv"
        if results.VERBOSE == True: print("Mapping reads to light chain variable region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --nofw  --no-hd -x " + database + results.genome + "_iglc " + \
              "--threads " + results.constant_value + " -U " + results.FASTQ + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".iglc"
        if results.VERBOSE == True: print("Mapping reads to light chain constant region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

elif paired == 1:
    if results.VERBOSE == True: print('Using Bowtie2 to find initial seeds:')
    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --norc --no-hd -x " + database + results.genome + "_ighv " + \
              "--threads " + results.constant_value + " -U " + results.LEFT + "," + results.RIGHT + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".ighv"
        if results.VERBOSE == True: print("Mapping reads to heavy chain variable region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive --quiet --nofw --no-hd -x " + database + results.genome + "_ighc " + \
              "--threads " + results.constant_value + " -U " + results.LEFT + "," + results.RIGHT + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".ighc"
        if results.VERBOSE == True: print("Mapping reads to heavy chain constant region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --norc  --no-hd -x " + database + results.genome + "_iglv " + \
              "--threads " + results.constant_value + " -U " + results.LEFT + "," + results.RIGHT + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}'  >" + output_location + "/" + tmp_name + ".iglv"
        if results.VERBOSE == True: print("Mapping reads to light chain variable region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

    try:
        cmd = results.bowtie + "/bowtie2 --very-sensitive  --quiet --nofw  --no-hd -x " + database + results.genome + "_iglc " + \
              "--threads " + results.constant_value + " -U " + results.LEFT + "," + results.RIGHT + " | awk -F'\t'  '{print $10\"\t\"$3\"\t\"$6\"\t\"$4}' >" + output_location + "/" + tmp_name + ".iglc"
        if results.VERBOSE == True: print("Mapping reads to light chain constant region ...")
        proc = subprocess.check_call(cmd, shell=True)
    except:
        print('Error Bowtie2 not installed or not found in path')
        print('Program terminated.')
        exit(0)

ptrn = output_location + "/" + tmp_name + ".ighv"
with open(ptrn, 'r') as f:
    line = f.readline()
ax, bx, cx, dx = line.split('\t')
read_length = len(ax)
ptrn = "\t*\t*\t"
d = defaultdict(int)

# anchor start of heavy chain
hv = defaultdict(int)
if results.VERBOSE: print("Searching for an anchor in heavy chain variable region ...")
cmd = output_location + "/" + tmp_name + ".ighv"
with open(cmd) as f:
    for line in f:
        ax, bx, cx, dx = line.split('\t')
        if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
            continue
        hv[ax] += 1
        hv[revcom(ax)] += 1
        if ptrn not in line:
            d[ax] += 1
try:
    hc_start = max(d, key=d.get)
    hc_start = hc_start.strip('N')
except (ValueError, TypeError):
    print(('Error: ' + tmp_name + ' did not map to heavy chain variable region.'))
    print('Program terminated.')
    exit(0)
if results.VERBOSE: print(("anchor: " + hc_start))
d.clear()

# anchor end of heavy chain
hc = defaultdict(int)
if results.VERBOSE: print("Searching for an anchor in heavy chain constant region ...")
cmd = output_location + "/" + tmp_name + ".ighc"
with open(cmd) as f:
    for line in f:
        ax, bx, cx, dx = line.split('\t')
        if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
            continue
        hc[ax] += 1
        hc[revcom(ax)] += 1
        if ptrn not in line:
            d[ax] += 1
try:
    hc_end = max(d, key=d.get)
    hc_end = hc_end.strip('N')
except (ValueError, TypeError):
    print(('Error: ' + tmp_name + ' did not map to heavy chain constant region.'))
    print('Program terminated.')
    exit(0)
if results.VERBOSE: print(("anchor: " + revcom(hc_end)))
d.clear()

# anchor start of light chain
lv = defaultdict(int)
if results.VERBOSE: print("Searching for an anchor in light chain variable region ...")
cmd = output_location + "/" + tmp_name + ".iglv"
with open(cmd) as f:
    for line in f:
        ax, bx, cx, dx = line.split('\t')
        if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
            continue
        lv[ax] += 1
        lv[revcom(ax)] += 1
        if ptrn not in line:
            d[ax] += 1
try:
    lc_start = max(d, key=d.get)
    lc_start = lc_start.strip('N')
except (ValueError, TypeError):
    print(('Error: ' + tmp_name + ' did not map to light chain variable region.'))
    print('Program terminated.')
    exit(0)
if results.VERBOSE: print(("anchor: " + lc_start))
d.clear()

# anchor end of light chain
lc = defaultdict(int)
if results.VERBOSE: print("Searching for an anchor in light chain constant region ...")
cmd = output_location + "/" + tmp_name + ".iglc"
with open(cmd) as f:
    for line in f:
        ax, bx, cx, dx = line.split('\t')
        if max(float(ax.count('A'))/read_length, float(ax.count('C'))/read_length, float(ax.count('G'))/read_length, float(ax.count('T'))/read_length) >= .5:
            continue
        lc[ax] += 1
        lc[revcom(ax)] += 1
        if ptrn not in line:
            d[ax] += 1
try:
    lc_end = max(d, key=d.get)
    lc_end = lc_end.strip('N')
except (ValueError, TypeError):
    print(('Error: ' + tmp_name + ' did not map to light chain constant region.'))
    print('Program terminated.')
    exit(0)
if results.VERBOSE: print(("anchor: " + revcom(lc_end)))
d.clear()
### End of single end sequencing processing

if __name__ == '__main__':
    p1 = Process(target=heavy,
                 args=(hc_end, hc_start, hc, hv, results.VERBOSE, read_length, results.name, output_location))
    p2 = Process(target=light,
                 args=(lc_end, lc_start, lc, lv, results.VERBOSE, read_length, results.name, output_location))
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    
    ptrn = output_location + "/" + tmp_name + ".ighv"
    os.remove(ptrn)
    ptrn = output_location + "/" + tmp_name + ".ighc"
    os.remove(ptrn)
    ptrn = output_location + "/" + tmp_name + ".iglv"
    os.remove(ptrn)
    ptrn = output_location + "/" + tmp_name + ".iglc"
    os.remove(ptrn)

print(("Done. Total run time: %s seconds" % (time.time() - start_time)))
