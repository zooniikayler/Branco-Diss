#Current script for apocrita

import mappy as mp
from output2 import SAMBAMWriter
import process_reads2
from multiprocessing import Process
from Bio import SeqIO


def get_output_var(hit):
    """
    Method to convert hit to dictonary to make mutable
    """
    out = {}
    out["ctg"] = hit.ctg
    out["r_st"] = hit.r_st
    out["r_en"] = hit.r_en
    out["q_st"] = hit.q_st
    out["q_en"] = hit.q_en
    out["mapq"] = hit.mapq
    out["cigar"] = hit.cigar_str
    out["strand"] = hit.strand
    return out

def filter_fastq(TE, R1, R2, out_fastq):
    '''Filter reads with a single read aligning to a given sequence
    '''
    reference = mp.Aligner(TE, preset="sr")  # load or build index
    if not reference:
        raise Exception("ERROR: failed to load/build index")

    out_fastq = open(out_fastq, "w")

    iterator1 = SeqIO.parse(R1, "fastq")
    iterator2 = SeqIO.parse(R2, "fastq")

    for r1 in iterator1:
        r2 = iterator2.__next__()

        r1_maps = list(reference.map(r1.seq))
        r2_maps = list(reference.map(r2.seq))

        if (len(r1_maps) >= 1) and not (len(r2_maps) >= 1):
            SeqIO.write(r2, out_fastq, 'fastq')
        elif not (len(r1_maps) >= 1) and (len(r2_maps) >= 1):
            SeqIO.write(r1, out_fastq, 'fastq')
    out_fastq.close()

def map_te_reads(read, refseq):
    '''Only take alignments with 1 hit'''
    alignment_objects = list(refseq.map(read.seq))
    if len(alignment_objects) == 1:
        out = [get_output_var(x) for x in alignment_objects]
    else:
        out = None
    return out

def run_polyte(reffile, r1name, fname, output_type, cut_site, min_len):
    '''Align filtered fastq to genome
    '''
    reference = mp.Aligner(reffile, preset="sr")
    print("Load in reference...") # load or build index
    if not reference:
        raise Exception("ERROR: failed to load/build index")
    print("Done")
    output_sam = SAMBAMWriter(fname, reference, output_type)
    print("Running alignment...")
    reads1 = mp.fastx_read(r1name)
    while True:
        try:
            read1 = process_reads2.Read(reads1.__next__())
            read1.split_read(cut_site, min_len)
            read1.qual_trim(10,10)
            if read1.seq:
                res = map_te_reads(read1, reference)
                if res:
                    output_sam.process_te_output(res, read1)
        except StopIteration:
            break



if __name__ == '__main__':
    p = Process(target=run_polyte, args=("from_scratch/test_genome.fa", "from_scratch/sample_R2.fq", "from_scratch/paired_test.sam", "SAM","GATC",20))
    q = Process(target=run_polyte, args=("from_scratch/test_genome.fa", "from_scratch/sample_R2.fq", "from_scratch/paired_test2.sam", "SAM","GATC",20))
    p.start()
    q.start()


