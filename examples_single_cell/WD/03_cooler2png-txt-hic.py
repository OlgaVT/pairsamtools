""" 
The script reads cooler files in target directory and converts them to txt, png and hicfiles. 

Note: It requires juicer java and juicer jar file to create .hic files for juicebox, for example: juicer_tools.1.8.9_jcuda.0.8.jar. Please, check the presence of this file and path to it. 
If you donit need .hic files, please, comment the corresponding line.

TODO: Rewrite so that it requires a single file input.

Example run:
python 03_cooler2png-txt-hic.py '../DATA/COOL/A6*.*.*'

"""

from sys import argv
import matplotlib as mpl
mpl.use('agg')

import glob
import os
from basic_utils import *

DIR_IMG = '../DATA/IMG'
DIR_TXT = '../DATA/TXT'
DIR_HIC = '../DATA/HIC'

if not os.path.isdir(DIR_IMG):
    os.mkdir(DIR_IMG)
if not os.path.isdir(DIR_TXT):
    os.mkdir(DIR_TXT)
if not os.path.isdir(DIR_HIC):
    os.mkdir(DIR_HIC)
    
#infile = argv[1]
#outfile = argv[2]

infiles_mask = argv[1] # '../DATA/COOL/A6*.*.*'
coolers = glob.glob(infiles_mask)

cmap = {0:'#241F20', 1:'#5C006D', 2:'#A30080', 3:'#FB5500', 4:'#FFB300', 5:'#F3F773'}
chromosomes = ['chr4', 'chrX', 'chr2L', 'chr2R', 'chr3L', 'chr3R']

for cool in coolers:
    logging.info("Processing cooler: {}".format(cool))
    
    basename = cool.split('/')[-1].split('.')
    cell = basename[0]
    res = basename[1]
    if len(basename)>3:
        exp = basename[3]
    else:
        exp = cell
        
    # if exp==cell or 'merged' in cool or 'Dros' in cool:
    #     for chrom in chromosomes:
    #        outfile = os.path.join(DIR_TXT, '{}.{}.{}.{}.txt'.format('mtx', cell, res, chrom))
    #        cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
    #        outfile = os.path.join(DIR_TXT, '{}.{}.{}.{}.txt'.format('sparse_bins', cell, res, chrom))
    #        cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
    #
    #
    # if res=='10' or res=='100' or res=='20' or res=='10000' or res=='20000' or res=='100000':
    #     for chrom in chromosomes:
    #         outfile = os.path.join(DIR_IMG, '{}.{}.{}.{}.png'.format(cell, exp, res, chrom))
    #         if 'Dros' in cool or 'merged' in cool:
    #           cmap='jet'
    #           balance=False
    #           if 'Dros' in cool:
    #             balance=True
    #             scale='log'
    #           remove_diagonal=True
    #           scale='linear'
    #         else:
    #           balance=False
    #           remove_diagonal=False
    #           scale='linear'
    #         cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom, remove_diagonal=remove_diagonal, balance=balance, scale=scale)


    cooler2hic(cool, os.path.join(DIR_HIC, '{}.{}.{}.hic'.format(cell, exp, res) ), 
               genome='dm3', resolutions=[10000,20000,100000, 1000000],
               remove_intermediary_files=True, juicer_path="~/soft/juicer/juicer_tools.1.8.9_jcuda.0.8.jar")
