#!/usr/bin/env python

import sys
import gzip
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

'''
This script reads a fastQ or a fastA file.
It is designed to work with a very specific 2-step PCR system.
In step 1, each primer will have the format:
5'<linker><UMI><target-binding-site>3'
The target-binding site is the actual primer that should bind to the
target of interest (i.e. your amplicon).
The UMI (univeral molecular identifier) is a random sequence of
fixed length in the PCR1 primer. If using only 2 cycles of the PCR,
any molecule of DNA from your sample seen in sequencing should be
identifyable using the UMI.
The linker is basically a priming site for PCR2. In this design, the
linker is identical for both the forward and reverse primers.
In PCR2 the PCR1 products are amplified using primers that target
the linker. The PCR2 primers would usually have a barcode and some other
buffering nucleotides.

In step 2 of the 2-step PCR, the linkers are used as targets for
barcoded primers which have the format:
<adapter><24nt barcode><linker>
The length of 24 nucleotides is hardcoded in this script.

This script looks for the linker in both orientations. It then
searches for where two consecutive instances of the linker occur in
the correct orientations to be an amplicon. A size filter may also
be applied in this script. The output of the script will have the
sequence from barcode to barcode with the format:
<24nt barcode><linker><UMI><primer><SEQ><primer><UMI><linker><24nt barcode>
Here, the barcodes should be the same for a "good product". The UMIs should
differ. SEQ refers to the biologically informative bit. The output
format will be in fastQ by default.
'''

# Standard IUPAC degenerate nucleotides
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

def find_amplicons(seqstr, primer, maxmismatch=2):
    primer = primer.upper()
    revmer = reverse_complement(primer)
    seqstr = seqstr.upper()
    f_matches = []
    r_matches = []
    cutoff = len(primer) - maxmismatch
    for i in range(len(seqstr) - len(primer)):
        subseq = seqstr[i: i+len(primer)]
        if match_score(subseq, primer) >= cutoff:
            f_matches.append(i)
        if match_score(subseq, revmer) >= cutoff:
            r_matches.append(i)
    all_matches = f_matches + r_matches
    all_matches = list(set(all_matches))
    all_matches.sort()
    coords = []
    for i in range(len(all_matches) -1):
        a = all_matches[i]
        b = all_matches[i+1]
        if a in f_matches and b in r_matches:
            coords.append([a, b])
    return coords

def child_rec_fasta(parent_rec, new_id, start_coord, stop_coord):
    seqstr = str(parent_rec[start_coord:stop_coord].seq)
    rec = SeqRecord(Seq(seqstr), id=new_id, name=new_id, description='')
    return rec

def child_rec_fastq(parent_rec, new_id, start_coord, stop_coord):
    rec = child_rec_fasta(parent_rec, new_id, start_coord, stop_coord)
    phreds = [elem for elem in parent_rec.letter_annotations['phred_quality'][start_coord:stop_coord]]
    rec.letter_annotations['phred_quality'] = phreds
    return rec

def child_rec(parent_rec, new_id, start_coord, stop_coord, seq_format):
    if seq_format == 'fasta':
        rec = child_rec_fasta(parent_rec, new_id, start_coord, stop_coord)
        return rec
    elif seq_format == 'fastq':
        rec = child_rec_fastq(parent_rec, new_id, start_coord, stop_coord)
        return rec
    else:
        raise ValueError("Sequence format (seq_format) must be one of 'fasta', or 'fastq'.")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Parses UMI sequences and coordinates of amplicons bound by PCR2 binding site"""
    )
    parser.add_argument('-i', '--infile', required=True,
        help="Input file: fastQ and fastA supported. May optionally be gzip compressed")
    parser.add_argument('-l', '--linker', required=True,
        help="Linker sequence - i.e. PCR2 adapter sequence")
    parser.add_argument('-u', '--umilen', required=True, type=int,
        help="UMI length: Must be an integer")
    parser.add_argument('-b', '--bclen', required=True, type=int,
        help="Barcode length: Must be an integer")
    parser.add_argument('-m', '--minlen', required=False, type=int, default=0,
        help="Minimum length of amplicon: Defaults to 0. Probably good to set it to a reasonable minimum. Length excludes linker and UMI")
    parser.add_argument('-M', '--maxlen', required=False, type=int, default=np.inf,
        help="Maximum length of amplicon: Defaults to infinity. Probably good to set it to a reasonable maximum. Length excludes linker and UMI")
    parser.add_argument('-f', '--fasta', action="store_true",
        help="Set this if fastA format is in use. Otherwise fastQ is assumed")
    

    args = parser.parse_args()
    INFILE = args.infile
    linker = args.linker
    MINLEN = args.minlen
    MAXLEN = args.maxlen
    UMILEN = args.umilen
    BCLEN = args.bclen
    SEQFMT = 'fastq'
    if args.fasta:
        SEQFMT = 'fasta'

    if INFILE.endswith('.gz'):
        fh = gzip.open(INFILE, 'rt')
    else:
        fh = open(INFILE, 'rt')

    progress = 0
    for rec in SeqIO.parse(fh, SEQFMT):
        progress += 1
        recid = str(rec.id)
        seqstr = str(rec.seq)
        recid = str(rec.id)
        coords = find_amplicons(seqstr, linker, maxmismatch=3)
        for i,coord in enumerate(coords):
            c1, c2 = coord[0], coord[1]
            amplen = (c2 - UMILEN) - (c1 + len(linker) + UMILEN)
            if amplen < MINLEN:
                continue
            if amplen > MAXLEN:
                continue
            start_coord = c1 - BCLEN
            stop_coord = c2 + len(linker) + BCLEN
            if start_coord < 0:
                continue
            if stop_coord >= len(rec):
                continue
            new_id = recid + '_' + str(i+1)
            newrec = child_rec(rec, new_id, start_coord, stop_coord, SEQFMT)
            _ = SeqIO.write([newrec], sys.stdout, SEQFMT)
        if progress % 10000 == 0:
            print(f"Processed {progress} records", file=sys.stderr)
