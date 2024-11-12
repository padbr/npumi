#!/usr/bin/env python

'''
This needs to be run on the output of `01_extract_amplicons.py`
or otherwise, the input must have the form:
<barcode><linker><UMI><primer><SEQUENCE><primer><UMI><linker><barcode>
  24 nt   15 nt  15nt  N nt    N nt      N nt   15nt  15nt    24nt

It will write the amplified biologically relevant portion of the
sequences to an output file in fastq format. It will annotate the
fastq headers with information about which UMI combination, barcode
and primers were found.
'''

import gzip
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# `primers' should be distributed with this script.
# It should be possible to manually edit the primers file.
from primers import linker, barcodes, f_primers, r_primers

degen_table = {
    'A':['A'],
    'T':['T'],
    'G':['G'],
    'C':['C'],
    'R' : ['A', 'G'],
    'Y' : ['C', 'T'],
    'M' : ['A', 'C'],
    'K' : ['G', 'T'],
    'S' : ['C', 'G'],
    'W' : ['A', 'T'],
    'H' : ['A', 'C', 'T'],
    'B' : ['C', 'G', 'T'],
    'V' : ['A', 'C', 'G'],
    'D' : ['A', 'G', 'T'],
    'N' : ['A', 'C', 'G', 'T']
}

comp_nucls = {
    'A':'T',
    'T':'A',
    'G':'C',
    'C':'G',
    'R' : 'Y',
    'Y' : 'R',
    'M' : 'M',
    'K' : 'K',
    'S' : 'S',
    'W' : 'W',
    'H' : 'D',
    'B' : 'V',
    'V' : 'B',
    'D' : 'H',
    'N' : 'N'
}

primer_masters = {
    'UMI_338Fa' : '338F',
    'UMI_338Fb' : '338F',
    'UMI_27Fa' : '27F',
    'UMI_27Fb' : '27F',
    'UMI_1391R' : '1391R',
    'UMI_1540R' : '1540R'
}

# The barcodes seem to be sequenced/basecalled with low accuracy.
# Therefore, barcode matching is done only for the 3'-half of
# the barcodes.
barcode_halves = { k : v[12:] for k,v in barcodes.items()}
all_primers = { k:v for k,v in f_primers.items() }
all_primers.update(r_primers)

def reverse_complement(seqstr):
    nucls = [comp_nucls[nucl] for nucl in seqstr]
    nucls.reverse()
    revc = ''.join(nucls)
    return revc

def match_nucls(n1, n2):
    '''Matches two nucleotides. n1: query (only ATGC) n2 (IUPAC degenerates allowed)'''
    if n1 in degen_table[n2]:
        return 1
    else:
        return 0

def match_score(subseq, primer):
    assert len(primer) == len(subseq)
    l1 = list(subseq)
    l2 = list(primer)
    vals = [match_nucls(n[0], n[1]) for n in zip(l1, l2)]
    return sum(vals)

def parse_parts(seqstr):
    bc1 = seqstr[:24]
    umi1 = seqstr[39:54]
    primer1s = [seqstr[54:71], seqstr[54:74], seqstr[54:75]]
    bc2 = seqstr[-24:]
    umi2 = seqstr[-54:-24]
    primer2s = [seqstr[-71:-54], seqstr[-74:-54], seqstr[-75:-54]]
    return bc1, umi1, primer1s, bc2, umi2, primer2s

def best_match(seqstrs, dbdict):
    bm_score = 0
    bm_str = ''
    bm_ord = ''
    for k,v in dbdict.items():
        for seqstr in seqstrs:
            try:
                mscore = match_score(seqstr, v)
            except AssertionError:
                continue
            if mscore > bm_score:
                bm_score = mscore
                bm_str = k
                bm_ord = '+'
    revcseqstrs = [reverse_complement(seqstr) for seqstr in seqstrs]
    for k,v in dbdict.items():
        for revcseqstr in revcseqstrs:
            try:
                mscore = match_score(revcseqstr, v)
            except AssertionError:
                continue
            if mscore > bm_score:
                bm_score = mscore
                bm_str = k
                bm_ord = '-'
    return bm_score, bm_ord, bm_str

def f_or_r_primer(pstr):
    if pstr in list(f_primers.keys()):
        return 'forward'
    elif pstr in list(r_primers.keys()):
        return 'reverse'
    else:
        raise ValueError(f"{pstr} not in primer list?")

def assess_rec(rec):
    seqstr = str(rec.seq).upper()
    bc1, umi1, primer1s, bc2, umi2, primer2s = parse_parts(seqstr)
    bc1full_score, bc1full_ord, bc1full_str = best_match([bc1], barcodes)
    bc1half_score, bc1half_ord, bc1half_str = best_match([bc1[12:]], barcode_halves)
    bc2full_score, bc2full_ord, bc2full_str = best_match([bc2], barcodes)
    bc2half_score, bc2half_ord, bc2half_str = best_match([bc2[:12]], barcode_halves)
    primer1_score, primer1_ord, primer1_str = best_match(primer1s, all_primers)
    primer2_score, primer2_ord, primer2_str = best_match(primer2s, all_primers)
    output = [bc1full_str, bc1full_score, bc1full_ord,
    bc2full_str, bc2full_score, bc2full_ord,
    bc1half_str, bc1half_score, bc1half_ord,
    bc2half_str, bc2half_score, bc2half_ord,
    primer1_str, primer1_score, primer1_ord,
    primer2_str, primer2_score, primer2_ord]
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Extracts biological sequence from amplicons, detailing barcodes and UMIs'
    )
    parser.add_argument('-i', '--infile', required=True,
        help='Input file: fastq format. Optionally gzip compressed')
    parser.add_argument('-o', '--outfile', required=True,
        help='Output file. Fastq format')
    parser.add_argument('-f', '--fasta', action="store_true",
        help='Input is in fastA format. So will the output. Otherwise fastQ will be used.')
    
    args = parser.parse_args()
    INFILE = args.infile
    OUTFILE = args.outfile
    
    if INFILE.endswith('.gz'):
        fh = gzip.open(INFILE, 'rt')
    else:
        fh = open(INFILE, 'rt')
    
    SEQFMT = 'fastq'
    if args.fasta:
        SEQFMT = 'fasta'
    
    outrecs = []
    for rec in SeqIO.parse(fh, 'fastq'):
        bc1full_str, bc1full_score, bc1full_ord, \
        bc2full_str, bc2full_score, bc2full_ord, \
        bc1half_str, bc1half_score, bc1half_ord, \
        bc2half_str, bc2half_score, bc2half_ord, \
        primer1_str, primer1_score, primer1_ord, \
        primer2_str, primer2_score, primer2_ord = assess_rec(rec)
        if bc1half_str != bc2half_str:
            continue
        if bc1half_ord == bc2half_ord:
            continue
        if primer1_ord == primer2_ord:
            continue
        if f_or_r_primer(primer1_str) == f_or_r_primer(primer2_str):
            continue
        umi1 = str(rec[39:54].seq)
        umi2 = str(rec[-54:-39].seq)
        if f_or_r_primer(primer1_str) == 'forward':
            start_coord = 54 + len(f_primers[primer1_str])
            end_coord = len(rec) - (54 + len(r_primers[primer2_str]))
            fprim = primer_masters[primer1_str]
            rprim = primer_masters[primer2_str]
            primer_group = fprim + '+' + rprim
            seqstr = str(rec[start_coord:end_coord].seq)
            umi_str = umi1 + '+' + umi2
            phred_scores = rec.letter_annotations['phred_quality'][start_coord:end_coord]
            desc = ' '.join([bc1half_str, umi_str, primer_group])
            newrec = SeqRecord(Seq(seqstr), id=str(rec.id), name=str(rec.id), description=desc)
            if SEQFMT == 'fastq':
                newrec.letter_annotations['phred_quality'] = phred_scores
            outrecs.append(newrec)
        elif f_or_r_primer(primer1_str) == 'reverse':
            start_coord = 54 + len(r_primers[primer1_str])
            end_coord = len(rec) - (54 + len(f_primers[primer2_str]))
            fprim = primer_masters[primer2_str]
            rprim = primer_masters[primer1_str]
            primer_group = fprim + '+' + rprim
            seqstr = str(rec[start_coord:end_coord].seq)
            seqstr = reverse_complement(seqstr)
            umi_str = reverse_complement(umi2) + '+' + reverse_complement(umi1)
            desc = ' '.join([bc1half_str, umi_str, primer_group])
            newrec = SeqRecord(Seq(seqstr), id=str(rec.id), name=str(rec.id), description=desc)
            if SEQFMT == 'fastq':
                phred_scores = rec.letter_annotations['phred_quality'][start_coord:end_coord]
                phred_scores.reverse()
                newrec.letter_annotations['phred_quality'] = phred_scores
            outrecs.append(newrec)
        else:
            print(f"Something went wrong with this record\n{str(rec.id)}\n{str(rec.seq)}")
            quit()
        
    
    _ = SeqIO.write(outrecs, OUTFILE, SEQFMT)
        




