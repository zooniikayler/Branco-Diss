import os

from itertools import zip_longest


def find_size(file):
    return os.stat(file).st_size


def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)


def batch(file, prefix):
    '''Split big boy files into parts'''
    if file:
        file_size = find_size(file)
    else:
        raise Exception("File not found")
    print("Original size... " + str(round(file_size/1073741824, 2)) + "GB")
        #file_size outputs bytes; 5368709120 bytes is 5GB
    if file_size > 5368709120:
        cut_num = 4*round((file_size/5368709120)/4)
        print("Batches... " + str(cut_num))
        lines_num = sum(1 for line in open(file))
        reads_num = lines_num/4
        print("Number of reads... " + str(reads_num))
        n = lines_num/cut_num
        print("Reads per batch... " + str(n/4))
        with open(file) as f:
            for i, g in enumerate(grouper(n, f, fillvalue=''), 1):
                with open('{}_{}.fq'.format(prefix,i), 'w') as fout:
                    fout.writelines(g)
    else:
        return print("File appropriate size")

batch("from_scratch/sample_R1.fq", "sample-R1")

from numpy import random

def genSeq(seqLen, numSeq, probs=True):
    '''Depends on random module'''
    bases = 'ACGTX'
    seqs = []
    for i in range(numSeq):
        seq = ''
        for j in range(seqLen):
            index = random.randint(0,4)
            seq += bases[index]
        seqs += [seq]
    return seqs

list_seq = genSeq(150, 10)

def genFastq(seqs, fname, seqlen, outtype, copynum):
    file_obj = open(fname, "w")
    read_num = 1
    qual = 'F' * seqlen
    if outtype == "fastq":
        for x in seqs:
            file_obj.write(("@read" + str(read_num) + "\n" + x + "\n" + "+" + "\n" + qual + "\n")*copynum)
            read_num = read_num + 1
        file_obj.close()
        print("File done")
    if outtype == "fasta":
        for x in seqs:
            file_obj.write(">Chrm" + str(read_num) + "\n" + x + "\n")
            read_num = read_num + 1
        file_obj.close()


genFastq(list_seq, "ref_MM_1.fa", 150, "fasta", 1)
genFastq(list_seq, "reads_MM_1.fq", 150, "fastq", 2)


