"""
Single cell data processing utilities. TODO: add annotation of functions. 
"""

import numpy as np
import pandas as pd
import h5py

from mirnylib import genome
import hiclib
from hiclib import fragmentHiC

import glob

import logging
logging.basicConfig(level=logging.INFO)


import time
from datetime import timedelta

import subprocess
import os

def call_and_check_errors(command):
    
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            shell=True, executable='/bin/bash')
    (stdout, stderr) = proc.communicate()
    logging.info("Check stdout: {}".format(stdout))
    if stderr:
        logging.info("Stderr is not empty. Might be an error in call_and_check_errors for the command: {}".format(command))
        logging.info("Check stderr: {}".format(stderr))
        return stderr   # Error, very bad!
    else:
        return 0        # No error, great!

def run_command(command, force=False):

    logging.info(command)

    possible_outfile = command.split('>')

    if len(possible_outfile)>1:
        possible_outfile = possible_outfile[-1]
        if os.path.isfile(possible_outfile):
            if force:
                logging.info("Outfile {} exists. It will be overwritten!".format(possible_outfile))
            else:
                raise Exception("Outfile {} exists. Please, delete it, or use force=True to overwrite it.".format(possible_outfile))

    cmd_bgn_time = time.time()
    is_err = call_and_check_errors(command)
    cmd_end_time = time.time()
    
    return is_err

def read_pairsams(filenames, experiment_ids, cell_name='', filter_type='JJ', chunk_size=1000):
    
    assert len(filenames)==len(experiment_ids)
    
    with open(filenames[0], 'r') as tmp_file:
        while True:
            line = tmp_file.readline()
            if not line.startswith('#'):
                break
            line_prev = line
        colnames = line_prev.split()[1:]
    
    df = []
    for exp, filename in zip(experiment_ids, filenames):
        iter_csv = pd.read_csv(filename, iterator=True, 
                               chunksize=chunk_size, comment='#', 
                               header=None, sep='\t', index_col=None)
        df_tmp = pd.concat([chunk[chunk.iloc[:,7] == filter_type] for chunk in iter_csv])
        df_tmp.columns = colnames
        df_tmp['exp'] = exp
        df_tmp['cell'] = cell_name
        df.append(df_tmp.copy())
        
    df = pd.concat(df).reset_index(drop=True)
    
    return df

def filter_pair_df(df, 
                   filter_mirrors=True, 
                   filter_exact_duplicates=True, 
                   filter_rfrags_type='rfrag_ends_count', 
                   max_rfrag_counts=2, rfrags_total=332981):
    
    df['rfrag_code_pair'] = df.rfrag1.values*rfrags_total + df.rfrag2.values 
    df['rfrag1_directed'] = df.rfrag1.values*10 + (df.strand1.values=='+')
    df['rfrag2_directed'] = df.rfrag2.values*10 + (df.strand2.values=='+')
    df['rfrag_code_pair_directed'] = df.rfrag1_directed.values*rfrags_total*10 + df.rfrag2_directed.values
    
    stats = {}
    
    N = len(df)
    logging.info("Initial size: {}".format(N))
    # Filter mirrors:
    if filter_mirrors:
        df = df.query('rfrag1!=rfrag2')
        logging.info("Mirrors filtered: {}".format(N-len(df)))
        stats['01_mirrors'] = N-len(df)
        N = len(df)
    
    #Filter duplicates
    if filter_exact_duplicates:
        df = df.drop_duplicates('rfrag_code_pair_directed')
        logging.info("Duplicates filtered: {}".format(N-len(df)))
        stats['02_duplicates'] = N-len(df)
        N = len(df)
    
    #Filter rfrags with too much contacts (>rfrag_end_count)
    if filter_rfrags_type=='rfrag_count':
        target_vector = np.concatenate([df.rfrag1.values.astype(int), df.rfrag2.values.astype(int)])
        v = np.bincount(target_vector)
        idx = np.where(v<=max_rfrag_counts)[0]
        v1 = np.in1d( df.rfrag1.values, idx )
        v2 = np.in1d( df.rfrag2.values, idx )
        df = df[v1&v2]

        logging.info("Filtering rfrags {}: {}".format(filter_rfrags_type, N-len(df)))
        stats['03_rfrag_filtered_{}'.format(filter_rfrags_type)] = N-len(df)
        N = len(df)
        
    elif filter_rfrags_type=='rfrag_ends_count':
        target_vector = np.concatenate([df.rfrag1_directed.values.astype(int), df.rfrag2_directed.values.astype(int)])
        v = np.bincount(target_vector)
        idx = np.where(v<=max_rfrag_counts)[0]
        v1 = np.in1d( df.rfrag1_directed.values, idx )
        v2 = np.in1d( df.rfrag2_directed.values, idx )
        df = df[v1&v2]

        logging.info("Filtering rfrags {}: {}".format(filter_rfrags_type, N-len(df)))
        stats['03_rfrag_filtered_{}'.format(filter_rfrags_type)] = N-len(df)
        N = len(df)
        
    elif filter_rfrags_type==None:
        pass
    else:
        logging.warn('filter_rfrags_type is not specified, skipping step: filter rfrags with too much contacts')

    logging.info("Resulting size: {}".format(N))
    stats['04_unique_contacts'] = N
    return df, stats
        
def create_pairix(df_input, output):
    """ #columns: readID chr1 pos1 chr2 pos2 strand1 strand2 """    
    logging.info("Creating pairix: {}".format(output))
    df_all = df_input.copy()
    df_all = df_all.query("(chrom1!='chrM')&(chrom2!='chrM')")
    df_all['cuts1'] = df_all.apply( lambda r: r.pos1-10 if r['strand1']=='+' else r.pos1+10, axis=1 )
    df_all['cuts2'] = df_all.apply( lambda r: r.pos2-10 if r['strand2']=='+' else r.pos2+10, axis=1 )
    
    df_all[['readID', 'chrom1', 'cuts1', 'chrom2', 'cuts2', 'strand1', 'strand2']].to_csv(output, index=False, header=False, sep='\t')
    return df_all

def create_cooler(df, pairix_file, cool_mask, chr_sizes, resolutions_list=[20, 100]):
    
    create_pairix(df, pairix_file)
    
    for res in resolutions_list:
        command1 = "cooler csort -c1 2 -c2 4 -p1 3 -p2 5 {pairix} {chr_sizes}".format(pairix=pairix_file, chr_sizes=chr_sizes)
        command2 = "cooler cload pairix -p 4 {chr_sizes}:{res} {pairix}.blksrt.gz {output}".format(res=res*1000, 
                                                                                   chr_sizes=chr_sizes, 
                                                                                   output=cool_mask.format(res),
                                                                                   pairix=pairix_file
                                                                                  )
        run_command(command1)
        run_command(command2)
        
def cooler2txt_chr(infile, outfile, fmt='mtx', chromosome='chrX'):

    import cooler
    import numpy as np
    
    if fmt=='mtx':
        c = cooler.Cooler(infile)
        mtx = c.matrix(balance=False, as_pixels=False).fetch(chromosome, chromosome)
        np.savetxt(outfile, mtx, fmt='%.0f')
    elif fmt=='sparse_bins':
        c = cooler.Cooler(infile)
        res = c.binsize
        mat_df = c.matrix(balance=False, as_pixels=True).fetch(chromosome, chromosome)
        
        # shiftig the bin ids so that the numeration is from 0
        chr_number = np.where(np.array(c.chromnames)==chromosome)[0][0]
        if chr_number==0:
            chr_start = 0
        else:
            chr_start = np.cumsum(c.chromsizes//res+1)[ chr_number-1 ]
        mat_df.bin1_id = mat_df.bin1_id-chr_start
        mat_df.bin2_id = mat_df.bin2_id-chr_start
        mat_df.to_csv(outfile, index=False)
        
    elif fmt=='sparse_coords':
        c = cooler.Cooler(infile)
        res = c.binsize
        mat_df = c.matrix(balance=False, as_pixels=True, join=True, ignore_index=False).fetch(chromosome, chromosome)
        mat_df.to_csv(outfile, index=False)
    else:
        raise Exception("Data saving format is unknown: {}".format(fmt))
        
def cooler2png_chr(infile, outfile, cmap='jet', chromosome='chrX'):
    import cooler
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors

    c = cooler.Cooler(infile)
    mtx = c.matrix(balance=False, as_pixels=False).fetch(chromosome, chromosome)
    
    fig = plt.figure(figsize=(20,20))
    if isinstance(cmap, str):
        plt.imshow(mtx, cmap='jet')
    elif isinstance(cmap, dict):
        #d = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
        cmap = {k:colors.to_rgb(cmap[k]) for k in cmap}
        mx = np.max(list(cmap))
        mtx[mtx>mx] = mx
        mtx_out = np.empty([mtx.shape[0], mtx.shape[1], 3])
        for k in cmap:
            mtx_out[mtx==k] = cmap[k]
        plt.imshow(mtx_out)
    else:
        raise Exception("color map unknown: {}".format(cmap))
    plt.xticks([])
    plt.yticks([])
    fig.tight_layout()
    plt.savefig(outfile, dpi=len(mtx)/20, bbox_inches='tight')
    plt.close(fig)