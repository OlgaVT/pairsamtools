"""
Single cell data processing utilities. TODO: add annotation of functions.

Works for Linux system only.

Non-trivial requirements for import:
h5py
mirnylib
hiclib
cooler

Non-trivial requirements for run:
java
juicer_tools .jar file
"""

import numpy as np
import pandas as pd
import h5py
import pickle

from mirnylib import genome
import hiclib
from hiclib import fragmentHiC

import glob
import cooler

import logging
logging.basicConfig(level=logging.INFO)


import time
from datetime import timedelta

import subprocess
import os


######### Basic utils to run linux commands #########

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


######### Single-cell Hi-C specific tools #########

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
                   max_rfrag_counts=4, rfrags_total=332981):
    
    df['rfrag_code_pair'] = df.rfrag1.values*rfrags_total + df.rfrag2.values 
    # Note that the pirs are always located R1-R2 (direction is also flipped in pairsamtools parse)
    # We need to save this information to retrieve rfrag ends
    df['rfrag1_directed'] = df.rfrag1.values*10 + (df.strand1.values=='+')
    df['rfrag2_directed'] = df.rfrag2.values*10 + (df.strand2.values=='-')
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

        
######### Files formats converters #########

def hdf2cool(infile, outfile, chrms_sizes, assembly='dm3', correct=True):
    """
    This function converts hiclib .hdf5 to .cool file.
    Note that attributes "heatmap" (whole-genome heatmap), "resolution" and "genomeIdxToLabel" are required in .hdf5 file.

    :param infile: input .hdf5 file
    :param outfile: output .cool file
    :param chrms_sizes: tab-separated file with chromosome lengths
    :param assembly: genome assembly (dm3 by default)
    :param correct: iteratively correct the heatmap? (True by default)

    :return: Python cooler object which was written to the file.
    """
    
    a = h5py.File(infile, 'r')
    heatmap = a['heatmap'].value

    chrms = list( map( lambda x: 'chr'+x if 'chr' not in x else x, pickle.loads(a['genomeIdxToLabel'].value).values() ))

    chromsizes = pd.read_csv(chrms_sizes, sep='\t', names=['name', 'length']).set_index('name').loc[chrms, 'length']
    binsize = pickle.loads(a['resolution'].value)
    bins    = cooler.binnify(chromsizes, binsize)

    iterator = cooler.io.ArrayLoader(bins, heatmap, binsize)
    cooler.io.create(outfile, bins, iterator, assembly=assembly)

    c = cooler.Cooler(outfile)
    if correct:
        bias, stats = cooler.ice.iterative_correction(c, store=c)
    
    return c
        
def cooler2txt_chr(infile, outfile, fmt='mtx', chromosome='chrX', writing_mode='w', separator=None):
    """
    Converts .cool to .txt matrix of three various formats (dense or sparse). For one chromosome.
    Writing_mode might be 'a' to append to the file, only for sparse data formats.

    :param infile: input .cool file
    :param outfile: output .txt file
    :param fmt: output format ('mtx' -- dense matrix, 'sparse_bins' -- sparse bins numbers, 'sparse_coords' -- sparse
        genomic coordinates). 'mtx' by default
    :param chromosome: name of chromosome to write into the file ('chr1', 'chr2', 'chrX' by default)
    :param writing_mode: only for sparse outputs. Any writing mode that can be used for pandas.to_csv 'mode' parameter.
        If "a" appends to the end of specified file without header created.
    :param separator: field separator for output files
    :return: None
    """

    if separator is None:
        separator = '\t'

    if fmt=='mtx': # dense matrix
        c = cooler.Cooler(infile)
        mtx = c.matrix(balance=False, as_pixels=False).fetch(chromosome, chromosome)
        np.savetxt(outfile, mtx, fmt='%.0f', delimiter=separator)

    elif fmt=='sparse_bins': # sparse matrix
        c = cooler.Cooler(infile)
        res = c.binsize
        mat_df = c.matrix(balance=False, as_pixels=True).fetch(chromosome, chromosome)
        
        # shifting the bin ids so that the numeration is from 0
        chr_number = np.where(np.array(c.chromnames)==chromosome)[0][0]
        if chr_number==0:
            chr_start = 0
        else:
            chr_start = np.cumsum(c.chromsizes//res+1)[ chr_number-1 ]
        mat_df.bin1_id = mat_df.bin1_id-chr_start
        mat_df.bin2_id = mat_df.bin2_id-chr_start
        mat_df.to_csv(outfile, index=False, mode=writing_mode, sep=separator, header=True if writing_mode=='w' else False)
        
    elif fmt=='sparse_coords':
        c = cooler.Cooler(infile)
        mat_df = c.matrix(balance=False, as_pixels=True, join=True, ignore_index=False).fetch(chromosome, chromosome)
        mat_df.to_csv(outfile, index=False, mode=writing_mode, sep=separator, header=True if writing_mode=='w' else False)
    else:
        logging.error("Data saving format is unknown: {}".format(fmt))


def cooler2png_chr(infile, outfile, cmap='jet', chromosome='chrX', remove_diagonal=0, balance=False, scale='linear'):
    """
    Converts cooler to png pixel image. The size of image is adjusted so that one cell of matrix corresponds to one pixel (TODO: check).

    :param infile: input .cool file
    :param outfile: output .png file
    :param cmap: color map for the heatmap. Can be either string with cmap name or ~matplotlib.colors.Colormap or
        a dictionary with all possible values in heatmap as key, 'jet' by default
    :param chromosome: chromosome to plot, 'chrX' by default
    :param remove_diagonal: number of diagonals to remove for plotting (to avoid too bright and noisy pixels), 0 by default
    :param balance: cooler balanced or not balanced heatmap should be loaded, False by default
    :param scale: scale of plot, 'log' or 'linear' (default)
    :return: None
    """
    import cooler
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib import colors

    c = cooler.Cooler(infile)
    mtx = c.matrix(balance=balance, as_pixels=False).fetch(chromosome, chromosome)
    
    if remove_diagonal:
      np.fill_diagonal(mtx, 0)

    if scale=='log':
      mtx = np.log2(mtx)

    fig = plt.figure(figsize=(20,20))

    if isinstance(cmap, str) or callable(cmap):
        plt.imshow(mtx, cmap=cmap)
    elif isinstance(cmap, dict):
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

def cooler2hic(cool, outfile_hic,
               genome = 'dm3',
               resolutions = [10000,20000,100000],
               remove_intermediary_files = False,
               juicer_path = './juicer_tools.1.8.9_jcuda.0.8.jar'):
    """
    Converts .cool to Lieberman's .hic file.

    :param cool: input .cool file
    :param outfile_hic: output .hic file
    :param genome: genome annotation name (one of hg18, hg19, hg38, dMel, mm9, mm10, anasPlat1, bTaurus3, canFam3,
        equCab2, galGal4, Pf3D7, sacCer3, sCerS288c, susScr3, or TAIR10) or a tab-delimited file with chromosomes sizes, default 'dm3'
    :param resolutions: list of resolutions that should be present in .hic file, default [10000,20000,100000]
    :param remove_intermediary_files: whether to remove intermediary .txt files, default False
    :param juicer_path: path to juicer .jar file, default './juicer_tools.1.8.9_jcuda.0.8.jar'
    :return: None
    """
    
    outfile_txt = outfile_hic + '.txt'
    outfile_tmp = outfile_hic + '.tmp'
    
    c = cooler.Cooler(cool)
    chromosomes = c.chromnames

    with open(outfile_tmp, 'w'):
        pass

    for chrom in chromosomes:
        cooler2txt_chr(cool, outfile_tmp,
                       fmt='sparse_coords',
                       chromosome=chrom,
                       writing_mode='a',
                       separator='\t')

    command1 = "awk '{{print 0, $1, $2, 0, 0, $4, $5, 1, $7}}' {} > {}".format(outfile_tmp, outfile_txt)
    command2 = "gzip -f {}".format(outfile_txt)
    command3 = "java -Xmx2g -jar {} pre -r {} {}.gz {} {}".format(juicer_path, ','.join(list(map(str, resolutions))), outfile_txt, outfile_hic, genome)

    run_command(command1)
    run_command(command2)
    run_command(command3)
    
    if remove_intermediary_files:
        os.remove(outfile_txt+'.gz')
        os.remove(outfile_tmp)


######### Cooler manipulations #########

def merge_single_cells(output, input_mask=None, files=None):
    """
    This functions creates a cooler file with merged single-cell data.
    Basically it's an explicit python wrapper over cooler merge.

    :param output: Filename of output cooler
    :param input_mask: Mask for all the files that are to be merged
    :param files: List of files to include into merged file
    :return: None
    """
    if input_mask:
        files = glob.glob(input_mask)

    command = "cooler merge {} {}".format(output, " ".join(files))
    run_command(command)


######### Scaling utilities #########
######### This part was adopted from other sources #########

from mirnylib.numutils import logbinsnew

def getMatrixScaling(inMatrix, inMask=[], measureType='sum', scaleType='log', logFactor=1.3):

    inMatrix = np.array(inMatrix, dtype=np.double)
    N = len(inMatrix)

    if len(inMask) > 0:
        mask2d = inMask
        inMatrix *= mask2d
    else:
        marginals = np.nansum(inMatrix, axis=0)
        mask = marginals > 0
        mask2d = mask[:, None] * mask[None, :]

    if scaleType == 'log':
        bins = logbinsnew(1, N, logFactor)
    else:
        bins = np.arange(0, N)

    mids = (0.5 * (bins[:-1] + bins[1:]))
    Pc = []
    for st, end in zip(bins[:-1], bins[1:]):
        curmean = 0
        maskmean = 0
        for i in range(st, end):
            if measureType == 'sum':
                curmean += np.nansum(np.diagonal(inMatrix, i))
                maskmean += np.nansum(np.diagonal(mask2d, i))
            else:
                curmean += np.nanmean(np.diagonal(inMatrix, i))
                maskmean += np.nanmean(np.diagonal(mask2d, i))

        Pc.append(curmean / maskmean)
    mids = np.r_[mids, N]
    Pc = np.r_[Pc, np.sqrt((Pc[-1] / Pc[-2])) * Pc[-1]]
    return Pc, mids

def getCoolerScaling(c, cc=None, scaleType='log', logFactor=1.15, chrom='all', measureType='sum'):

    Pc_list = []
    mids_list = []

    if chrom == 'all':
        chrs = c.chromnames
    else:
        chrs = chrom

    for chrom in chrs:
        try:
            inMatrix = c.matrix(balance=False).fetch('{0}'.format(chrom))

            if cc is None:
                Pc, mids = getMatrixScaling(inMatrix, measureType=measureType, logFactor=logFactor, scaleType=scaleType)
            else:
                marginals = np.nansum(cc.matrix(balance=False).fetch('{0}'.format(chrom)), axis=0)
                mask = marginals > 0
                mask2d = mask[:, None] * mask[None, :]
                Pc, mids = getMatrixScaling(inMatrix, inMask=mask2d, measureType=measureType, logFactor=logFactor, scaleType=scaleType)

            Pc_list.append(Pc)
            mids_list.append(mids)
        except Exception as e:
            print('Could not process chromosome:{0} Error: \n {1}'.format(chrom, e))

    # get average value genome-wide
    biggest_val = np.max([len(x) for x in Pc_list])

    Pc = np.zeros((len(Pc_list), biggest_val)) * np.nan
    for si, s in enumerate(Pc_list):
        Pc[si, 0:len(s)] = s
    Pc = np.nanmean(Pc, axis=0)

    mids = mids_list[0]
    for m in mids_list:
        if len(m) > len(mids):
            mids = m
    return Pc, mids, Pc_list

def get_scalings_df(files, labels, merged_file=None, chromosomes='all', scaleType='log', logFactor=1.15, measureType='sum'):
    """
    Calculates scalings for all the files provided and creates output dataframe.
    :param files:
    :param labels:
    :param merged_file:
    :param chromosomes:
    :param scaleType:
    :param logFactor:
    :param measureType:
    :return: dataframe with columns ['Pc', 'mids', 'label']
    """
    assert len(files)==len(labels)

    df = {x: [] for x in ['Pc', 'mids', 'label']}

    if not merged_file is None:
        merged_cooler  = cooler.Cooler(merged_file)
    else:
        merged_cooler = None

    for filename, label in zip(files, labels):
        logging.info("Reading cooler: {}".format(filename))

        c = cooler.Cooler(filename)

        Pc, mids, Pc_list = getCoolerScaling(c, merged_cooler, chrom=chromosomes, logFactor=logFactor, scaleType=scaleType, measureType=measureType)

        df['Pc'] += list(Pc)
        df['mids'] += list(mids)
        df['label'] += [label]*len(Pc)

    df = pd.DataFrame(df)

    return df


######### TADs utilities #########

import lavaburst

def produce_segmentation(mtx, gamma, good_bins='default', method='modularity', max_intertad_size=3, max_tad_size=10000):
    """
    Produces single segmentation (TADs or CDs calling) of mtx with one gamma with the algorithm provided.
    :param mtx: input numpy matrix
    :param gamma: parameter for segmentation calling
    :param good_bins: bool vector with length of len(mtx) with False corresponding to masked bins, 'default' is that
        good bins are all columns/rows with sum > 0
    :param method: 'modularity' (default) or 'armatus'
    :param max_intertad_size: max size of segmentation unit that is considered as interTAD
    :return:  2D numpy array where segments[:,0] are segment starts and segments[:,1] are segments end, each row corresponding to one segment
    """
    if np.any(np.isnan(mtx)):
        logging.warning("NaNs in dataset, pease remove them first.")
    
    if np.diagonal(mtx).sum()>0:
        logging.warning("Note that diagonal is not removed. you might want to delete it to avoid noisy and not stable results. ")
    
    if method=='modularity':
        score = lavaburst.scoring.modularity_score
    elif method=='armatus':
        score = lavaburst.scoring.armatus_score
    else:
        pass
        
    if good_bins=='default':
        good_bins = mtx.astype(bool).sum(axis=0) > 0

    S = score(mtx, gamma=gamma, binmask=good_bins)
    model = lavaburst.model.SegModel(S)

    segments = model.optimal_segmentation()
    
    v = segments[:,1]-segments[:,0]
    mask = (v>max_intertad_size) & (np.isfinite(v)) & (v<max_tad_size)

    segments = segments[mask]
    
    return segments

def find_optimal_gamma(mtx, optimum_nbins, 
                       epsilon=0.2, 
                       gammas=[1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 12.0, 14.0, 15.0, 20.0], 
                       good_bins='default', 
                       method='modularity', 
                       max_intertad_size=3, 
                       max_tad_size=10000,
                       optimization_function=np.median, 
                       min_tads_number=10):
    """
    Finds optimal gamma that produces the closest segments size to what is expected.
    :param mtx: input numpy matrix
    :param optimum_nbins: expected size of TADs/CCs/segments
    :param epsilon: the precision with which expected size should be approximated, note that this precision is not
        guaranteed but algorithmically it improves computation time for the wort case
    :param gammas: list of gammas to probe the segmentation calling
    :param good_bins: bool vector with length of len(mtx) with False corresponding to masked bins, 'default' is that
        good bins are all columns/rows with sum > 0
    :param method: 'modularity' (default) or 'armatus'
    :param max_intertad_size: max size of segmentation unit that is considered as interTAD
    :return: the best approximated mean segments size and optimal gamma for this segmentation
    """

    means = []
    
    for n, gamma in enumerate(gammas):
        segments = produce_segmentation(mtx, gamma, good_bins=good_bins, method=method, max_intertad_size=max_intertad_size, max_tad_size=max_tad_size)
        
        v = segments[:,1]-segments[:,0]
        if len(v)<min_tads_number:
            means.append(0)
            continue

        means.append(optimization_function(v))
        
        if np.abs(np.array(optimization_function(v)-optimum_nbins))<epsilon:
            break
            
    if len(means)==0:
      logging.error("TADs calling run with bad parameters set: empty means array!")

    idx = np.argmin(np.abs(np.array(means)-optimum_nbins))
    opt_mean, opt_gamma = means[idx], gammas[idx]
        
    return opt_mean, opt_gamma

def find_optimal_gamma_2opt(mtx, optimum_nbins,  
                       gammas=[1.5, 2.0, 3.0, 4.0, 6.0, 10.0, 12.0, 14.0, 15.0, 20.0], 
                       good_bins='default', 
                       method='modularity', 
                       max_intertad_size=3, 
                       max_tad_size=10000,
                       optimization_function=np.median):
    """
    Finds optimal gamma that produces the closest segments size to what is expected, but with two optimization steps.
    First, optimization_function for all segmentations is minimized. Second, over all best cases the segmentation with 
    largest number of TADs is selected. 
    
    The reasoning behind this implementation is that we want not only to obtain segmentation with the best medium or 
    mean size of TADs, but also get the best coverage of chromosome with TADs. 
    
    :param mtx: input numpy matrix
    :param optimum_nbins: expected size of TADs/CCs/segments
    :param gammas: list of gammas to probe the segmentation calling
    :param good_bins: bool vector with length of len(mtx) with False corresponding to masked bins, 'default' is that
        good bins are all columns/rows with sum > 0
    :param method: 'modularity' (default) or 'armatus'
    :param max_intertad_size: max size of segmentation unit that is considered as interTAD
    :return: the best approximated mean segments size and optimal gamma for this segmentation
    """

    means = []
    ntads = []

    for n, gamma in enumerate(gammas):
        segments = produce_segmentation(mtx, gamma, 
                                        good_bins='default', method=method, 
                                        max_intertad_size=max_intertad_size, max_tad_size=max_tad_size)

        v = segments[:,1]-segments[:,0]

        means.append(optimization_function(v))
        ntads.append(len(v))
        
    diff = np.abs(np.array(means)-optimum_nbins)
    idxs = np.where( np.abs(diff-min(diff))<0.0001)[0]
    idx = idxs[ np.argmax(np.array(ntads)[idxs]) ]
    opt_mean, opt_gamma, opt_ntads = means[idx], gammas[idx], ntads[idx]
        
    return opt_mean, opt_gamma, opt_ntads

def segmentations_to_bed(by_chr_dct, outfile, resolution):
    """
    Writes dictionary with segmentations to bed file.

    :param by_chr_dct: dictionary with segments, keys are chromosomes names, values are numpy 2D arrays
    :param outfile:
    :param resolution:
    :return: None
    """
    with open(outfile, 'w') as outf:
        for k in by_chr_dct.keys():
            for i in by_chr_dct[k]:
                outf.write("\t".join([k, str(i[0]*resolution),str(int(i[1]-0.5)*resolution)]) + "\n")
                
def segmentations_to_2Djuice(by_chr_dct, outfile, resolution):
    """
    Writes dictionary with segmentations to 2Dannotation Juicebox file.
    :param by_chr_dct: dictionary with segments, keys are chromosomes names, values are numpy 2D arrays
    :param outfile:
    :param resolution:
    :return: None

# .2Dannot file example:
# chr1   x1         x2         chr2   y1         y2         color     comment
# chrX   85000000   89000000   chrX   85000000   89000000   0,255,0   My green region
# chrX   90000000   99100000   chrX   90000000   99100000   0,0,255   My blue region
    """

    with open(outfile, 'w') as outf:
        outf.write("\t".join("chr1 x1 x2 chr2 y1 y2 color comment".split())+"\n")
        for k in by_chr_dct.keys():
            for i in by_chr_dct[k]:
                line = [k, str(i[0]*resolution),str((i[1])*resolution), k, str(i[0]*resolution),str((i[1])*resolution), "0,0,255", "TAD"]
                outf.write("\t".join(line) + "\n")


######### Insulation Score (IS) utilities #########

from cooltools import insulation
from functools import reduce

def calculate_insulation_for_ranges(cool, square_size_range, diagonals_range, resolution):
    """
    Calculates insulation for range of values of parameters for insulation score. Basically, this is a wrapper over cooltools insulation function
    :param cool: input .cool file
    :param square_size_range:
    :param diagonals_range:
    :param resolution:
    :return:
    """
    
    dfs = [ insulation.find_insulating_boundaries(cool, window_bp = x, ignore_diags=y)\
        .drop(['bad_bin_masked', 'boundary_strength_{}'.format(x)], axis=1)\
        .rename(columns = {'chrom':'chrom', 'start':'start', 'end':'end', 
                           'log2_insulation_score_{}'.format(x) : 'log2_IS_sq:{}_ind:{}'.format(x, y*resolution)} ) for x in square_size_range for y in diagonals_range ]

    df = reduce(lambda left,right: pd.merge(left,right,on=['chrom', 'start', 'end']), dfs)
    
    return df


######### Permutations utilities #########

def permute_matrix(mtx):
    pass

def permute_segmentation(segmentations):
    pass

######### Saddle plot utilities #########
