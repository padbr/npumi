#!/usr/bin/env python

import sys
import gzip
import random
import argparse
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_umi(rec):
    '''Parse a UMI string for a sequence record header'''
    return str(rec.description).split(' ')[2]

def get_prims(rec):
    '''Parse the primer-pair string from a sequence record header'''
    return str(rec.description).split(' ')[3]

def random_xfix(n=6):
    '''Returns 6 random letters to use as a basename for temporary files'''
    return ''.join([random.choice([chr(97 + i) for i in range(26)]) for _ in range(n)])

def del_temp_files(Flist):
    for F in Flist:
        cmd = ['rm', F]
        cmd = ' '.join(cmd)
        try:
            _ = subprocess.check_call(cmd, shell=True)
        except:
            pass

def get_consensus(recs):
    xfix = random_xfix()
    xfix = 'temp.' + xfix
    rawreads = xfix + '.raw.fastq'
    minimap_aln = xfix + '.mp2.sam'
    racon_out = xfix + '.racon.fasta'
    clustal_out = xfix + '.clustalo.fasta'
    _ = SeqIO.write(recs, rawreads, 'fastq')
    # all vs all minimap2
    cmd = ['minimap2', '-x', 'map-ont', '-t', '6', '-a', '-o', minimap_aln, rawreads, rawreads]
    cmd = ' '.join(cmd)
    _ = subprocess.check_call(cmd, shell=True)
    # racon
    cmd = ['racon', '-f', '-t', '6', rawreads, minimap_aln, rawreads]
    cmd = ' '.join(cmd)
    try:
        output = subprocess.check_output(cmd, shell=True)
    except subprocess.CalledProcessError:
        # clean up
        del_temp_files([rawreads, minimap_aln, racon_out, clustal_out])
        return None
    output = output.decode('utf-8')
    output = output.strip()
    with open(racon_out, 'wt') as ofh:
        _ = ofh.write(output)
        ofh.close()
    racon_recs = list(SeqIO.parse(racon_out, 'fasta'))
    racon_seqs = list(set([str(rec.seq) for rec in racon_recs]))
    del_temp_files([rawreads, minimap_aln, racon_out, clustal_out])
    if len(racon_seqs) == 1:
        return racon_seqs[0]
    else:
        return None


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Sorts an annotated fastq file into clusters of sequences based on barcodes and UMIs and attempts to generate consensus ASVs.'
    )
    parser.add_argument('-i', '--infile', required=True,
        help='Input file. Must be in fastQ format. Optionally gzip compressed.')
    parser.add_argument('-o', '--outfile', required=True,
        help='Output file name. Will be in fastA format.')
    parser.add_argument('-m', '--min_cluster_seqs', default=6, type=int,
        help='Minimum number of sequences required to attempt clustering to consensus. Default=6')
    parser.add_argument('-p', '--primer_proportion', default=0.9, type=float,
        help='At least this proportion of the primer pairs must be the same to avoid filtering on mismatching primers. Default=0.9')
    
    args = parser.parse_args()
    INFILE = args.infile
    OUTFILE = args.outfile
    MIN_UMI_CLUST_LEN = int(args.min_cluster_seqs)
    SEQFMT = 'fastq'
    PROPORTION_PRIMER = float(args.primer_proportion)
    if PROPORTION_PRIMER <= 0.5:
        print("Error:\nPlease reset the '--primer_proportion' argument to greater than 0.5",
              file=sys.stderr)
        quit()
        '''
        At least this proportion of the primers for a UMI must be the same.
        Otherwise the UMI is discarded.
        Must be greater than 0.5 to be deterministic.
        '''
    
    # Parse infile and get lists of UMI strings and primer strings
    if INFILE.endswith('.gz'):
        in_fh = gzip.open(INFILE, 'rt')
    else:
        in_fh = open(INFILE, 'rt')
    recs = list(SeqIO.parse(in_fh, SEQFMT))
    umis = [get_umi(rec) for rec in recs]
    prims = [get_prims(rec) for rec in recs]
    
    # Count the occurrences of each unique UMI string and sort in descending order
    unq_umis = list(set(umis))
    counts = [umis.count(umi) for umi in unq_umis]
    umi_counts = [[umi, count] for umi, count in zip(unq_umis, counts)]
    umi_counts.sort(key=lambda x:x[1], reverse=True)

    if len(umi_counts) == 0: # This could happen if the input is an empty file
        print(f"No UMI clusters found for {OUTFILE}", file=sys.stderr)
        quit()

    # If UMI count is below a limit there aren't enough of them to get a consensus sequence
    while umi_counts[-1][1] < MIN_UMI_CLUST_LEN:
        _ = umi_counts.pop()
        if len(umi_counts) == 0:
            print(f"No UMI clusters for {OUTFILE}", file=sys.stderr)
            quit()
    
    '''
    Sort through the UMIs and filter out chimeras.
    If PCR was perfect, each UMI and combination of two UMIs on either end
    should have an astronomically low chance of being seen more than once.
    Therefore, UMIs seen more than once will be assumed to be indicative of
    chimeras. The uniquely most abundant instance of a UMI is being accepted,
    with the assumption that it is the parent of chimeras and not a chimera
    itself.
    '''
    good_umis = []
    bad_fumis, bad_rumis = [], []
    good_frumis = []
    while umi_counts:
        query = umi_counts.pop(0)
        qfumi, qrumi = query[0].split('+')
        qcount = query[1]
        for subject in umi_counts:
            if subject[1] < qcount:
                break # less abundant umis will be excluded when it is seen that there is already a more abundant one
            sfumi, srumi = subject[0].split('+')
            if qfumi == sfumi or qrumi == srumi:
                bad_fumis.append(qfumi)
                bad_rumis.append(qrumi)
        if qfumi in bad_fumis or qrumi in bad_rumis:
            bad_fumis.append(qfumi)
            bad_rumis.append(qrumi)
        elif qfumi in good_frumis or qrumi in good_frumis:
            bad_fumis.append(qfumi)
            bad_rumis.append(qrumi)
        else:
            good_umis.append(query)
            good_frumis.append(qfumi)
            good_frumis.append(qrumi)
        bad_fumis = list(set(bad_fumis))
        bad_rumis = list(set(bad_rumis))
        good_frumis = list(set(good_frumis))
    
    print(f"Found {len(good_umis)} good UMIs", file=sys.stderr)
    for good_umi in good_umis:
        print(good_umi, file=sys.stderr)
    
    # Gather together records for the good UMIs and filter out ones with (too many) mismatched primers
    umi_recss = []
    for umi_count in good_umis:
        umi = umi_count[0]
        umi_recs = [rec for rec in recs if get_umi(rec) == umi]
        assert len(umi_recs) == umi_count[1], f"Problem with this UMI:\n{umi}"
        prims = [get_prims(rec) for rec in umi_recs]
        unq_prims = list(set(prims))
        prim_counts = [prims.count(prim) for prim in unq_prims]
        if max(prim_counts) / sum(prim_counts) < PROPORTION_PRIMER:
            print(f"Dropping UMI for clashing primers: {umi}", file=sys.stderr)
            continue
        i = prim_counts.index(max(prim_counts))
        prim = prims[i]
        umi_recs = [rec for rec in umi_recs if get_prims(rec) == prim]
        if len(umi_recs) < MIN_UMI_CLUST_LEN:
            print(f"Dropping {umi} due to low abundance after filtering out bad primers.",
                  file=sys.stderr)
            continue
        umi_recss.append(umi_recs)
    
    print(f"Ther are now {len(umi_recss)} good UMI clusters", file=sys.stderr)
    
    consensus_recs = []

    for i,umi_recs in enumerate(umi_recss):
        consensus = get_consensus(umi_recs)
        if consensus == None:
            continue
        desc = str(umi_recs[0].description)
        recid = 'asv' + str(i+1)
        cons_rec = SeqRecord(Seq(consensus), id=recid, description=desc)
        consensus_recs.append(cons_rec)
    
    if len(consensus_recs) >= 1:
        i = SeqIO.write(consensus_recs, OUTFILE, 'fasta')
        print(f"Wrote {i} sequences to {OUTFILE}", file=sys.stderr)
    else:
        print(f"No output for {OUTFILE}", file=sys.stderr)

   
    
        
