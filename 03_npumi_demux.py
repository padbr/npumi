#!/usr/bin/env python

import os
import gzip
import argparse
from Bio import SeqIO

'''
This is intended to be run on the output of '02_annotate_amplicons.py'.
The header of the input files must have the format:
Sequence_identifier Barcode_identifier UMI_combination Primers_used
These must be separated by single spaces.
'''

def get_fh(infile):
    if infile.endswith('.gz'):
        return gzip.open(infile, 'rt')
    else:
        return open(infile, 'rt')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Reads all files in current directory with a particular ending and sorts the reads by barcodes annotated in the headers.'
    )
    parser.add_argument('-e', '--ending', required=True,
        help='The ending of all files to be scanned in the current directory. FastA/Q supported. Optionally gzip compressed.')
    parser.add_argument('-n', '--numbarcodes', type=int, default=192,
        help='The number of potential barcodes. Integer only. Default: 192')
    parser.add_argument('-f', '--fasta', action="store_true",
        help='Use this flag to specify fastA input/output. Otherwise fastQ will be used.')
    parser.add_argument('-l', '--lowmem', action="store_true",
        help='Use this flag if your system is running out of memory. It will be much slower with this option.')
    
    args = parser.parse_args()
    ENDING = args.ending
    NUMBCS = int(args.numbarcodes)
    SEQFMT = 'fastq'
    if args.fasta:
        SEQFMT = 'fasta'

    bc_indices = {} # Barcode ID : index in recss
    recss = []
    recss_indices = []

    infiles = [File for File in os.listdir(os.getcwd()) if str(File).endswith(ENDING)]
    if len(infiles) == 0:
        print(f"Error: There are no files ending with '{ENDING}' in the current directory.")
        quit()
    
    if args.lowmem:
        print("Running in low memory mode.\nThis will be slow!\n")
        print("Parsing barcode identifiers.")
        barcodes = []
        for infile in infiles:
            fh = get_fh(infile)
            recs = list(SeqIO.parse(fh, SEQFMT))
            for rec in recs:
                bc_id = str(rec.description).split(' ')[1]
                if not bc_id in barcodes:
                    barcodes.append(bc_id)
            fh.close()
        
        for barcode in barcodes:
            print(f"Addressing {barcode}")
            outrecs = []
            for infile in infiles:
                fh = get_fh(infile)
                recs = list(SeqIO.parse(fh, SEQFMT))
                for rec in recs:
                    bc_id = str(rec.description).split(' ')[1]
                    if bc_id == barcode:
                        outrecs.append(rec)
                fh.close()
            ofn = barcode + '.' + SEQFMT
            _ = SeqIO.write(outrecs, ofn, SEQFMT)

    
    else:
        print("Running in fast mode.")
        print("If too much memory is being consumed, stop and rerun with '-l' switch to reduce memory use.")
        for infile in infiles:
            fh = get_fh(infile)
            recs = list(SeqIO.parse(fh, SEQFMT))
            for rec in recs:
                bc_id = str(rec.description).split(' ')[1]
                if not bc_id in list(bc_indices.keys()):
                    recss.append([])
                    i = len(recss) - 1
                    bc_indices[bc_id] = i
                i = bc_indices[bc_id]
                recss[i].append(rec)
            fh.close()
        
        for bc_id, i in bc_indices.items():
            ofn = bc_id + '.' + SEQFMT
            recs = recss[i]
            _ = SeqIO.write(recs, ofn, SEQFMT)
