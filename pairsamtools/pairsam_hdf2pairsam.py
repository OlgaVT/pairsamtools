#!/usr/bin/env python
# -*- coding: utf-8 -*-
from collections import OrderedDict
import subprocess
import fileinput
import itertools
import click
import pipes
import sys
import os
import io
import copy
import h5py
import pickle

from . import _fileio, _pairsam_format, _headerops, cli, common_io_options
from .pairsam_stats import PairCounter


UTIL_NAME = 'pairsam_hdf2pairsam'

EXTRA_COLUMNS = [
    'mapq',
    'pos5',
    'pos3',
    'cigar',
    'read_len',
    'matched_bp',
    'algn_ref_span',
    'algn_read_span',
    'dist_to_5',
    'dist_to_3',
    'rfrag',
    'rfrag_dist',
    'rfrag_dist_up',
    'rfrag_dist_down'
]

@cli.command()
@click.argument(
    'hdf_path',
    type=str,
    required=False)
@click.option(
    "-o", "--output", 
    type=str, 
    default="", 
    help='output file. '
        ' If the path ends with .gz or .lz4, the output is pbgzip-/lz4-compressed.'
         'By default, the output is printed into stdout. ')
@click.option(
    "-f", "--frags",
    type=str,
    required=False,
    help='a tab-separated BED file with the positions of restriction fragments '
         '(chrom, start, end). Can be generated using cooler digest.')

@common_io_options

def hdf2pairsam(hdf_path, output, **kwargs):
    '''parse .hdf5 and make .pairsam.

    SAM_PATH : input .sam file. If the path ends with .bam, the input is 
    decompressed from bam. By default, the input is read from stdin.
    '''
    parse_hdf(hdf_path, output, **kwargs)


def parse_hdf(hdf_path, output, **kwargs):


    infile = h5py.File(hdf_path)

    outstream = (_fileio.auto_open(output, mode='w',
                                   nproc=kwargs.get('nproc_out'),
                                   command=kwargs.get('cmd_out', None)) 
                 if output else sys.stdout)

    write_pairsam(infile, outstream, **kwargs)

    if outstream != sys.stdout:
        outstream.close()

def write_pairsam(infile, out_file, **kwargs):

    infile_d = {"chrms1": infile['chrms1'].value,
                "chrms2": infile['chrms2'].value,
                "pos1": infile['cuts1'].value,
                "pos2": infile['cuts2'].value,
                "strand1": infile['strands1'].value,
                "strand2": infile['strands2'].value,
                }

    idx2label = pickle.loads(infile['misc'].value)['genome']['idx2label']

    #reading rfrags
    if len(kwargs['frags'])>0:
        frags = kwargs['frags']

        import numpy as np
        from numpy.lib.recfunctions import append_fields  # for rfrags indexing

        rfrags = np.genfromtxt(
            frags, delimiter='\t', comments='#', dtype=None,
            names=['chrom', 'start', 'end', 'idx'])

        rfrags.sort(order=['chrom', 'start', 'end'])

        rfrags = append_fields(rfrags, 'idx', np.arange(len(rfrags)))
        rfrags['end'] += 1

        chrom_borders = np.r_[0,
                              1 + np.where(rfrags['chrom'][:-1] != rfrags['chrom'][1:])[0],
                              rfrags.shape[0]]
        rfrags = {rfrags['chrom'][i]: rfrags[['end', 'idx']][i:j]
                  for i, j in zip(chrom_borders[:-1], chrom_borders[1:])}

        print('Rfrags read')


    out_file.write("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rfrag1 rfrag2")
    out_file.write("\n")

    for i in range(len(infile_d["chrms1"])):

        if (infile_d['chrms1'][i]<0) or (infile_d['chrms2'][i]<0):
            continue

        chr1 = "chr"+idx2label[infile_d['chrms1'][i]]
        chr2 = "chr"+idx2label[infile_d['chrms2'][i]]

        rfrag1, _, _ = \
            find_rfrag(rfrags, chr1, infile_d['pos1'][i] + (10 if infile_d['strand1'][i] else -10))
        rfrag2, _, _ = \
            find_rfrag(rfrags, chr2, infile_d['pos2'][i] + (10 if infile_d['strand2'][i] else -10))

        if rfrag1<rfrag2:
            towrite = [i, chr1, infile_d['pos1'][i], chr2, infile_d['pos2'][i],
                      "+" if infile_d['strand1'][i] else "-",  "+" if infile_d['strand2'][i] else "-", 'JJ', rfrag1, rfrag2]
        else:
            towrite = [i, chr2, infile_d['pos2'][i], chr1, infile_d['pos1'][i],
                      "+" if infile_d['strand2'][i] else "-",  "+" if infile_d['strand1'][i] else "-", 'JJ', rfrag2, rfrag1]

        for x in towrite:
            out_file.write(str(x))
            out_file.write(_pairsam_format.PAIRSAM_SEP)

        out_file.write("\n")

    out_file.write('\n')


if __name__ == '__main__':
    parse_hdf()

def find_rfrag(rfrags, chrom, pos):
    if chrom.encode('ascii') in rfrags.keys():
        rsites_chrom = rfrags[chrom.encode('ascii')]['end']
        rsites_idx = rfrags[chrom.encode('ascii')]['idx']
        idx = min(max(0, rsites_chrom.searchsorted(pos, 'right') - 1), len(rsites_chrom) - 2)
        return rsites_idx[idx], rsites_chrom[idx], rsites_chrom[idx + 1]
    else:
        return -1, -1, -1