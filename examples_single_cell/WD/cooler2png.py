""" 
The script reads cooler files in target directory and converts them to txt and png files. 
"""

#from sys import argv
import matplotlib as mpl
mpl.use('agg')

import glob
import os
from sc_lib import *

if not os.path.isdir('../DATA/IMG'):
    os.mkdir('../DATA/IMG')
if not os.path.isdir('../DATA/TXT'):
    os.mkdir('../DATA/TXT')

#infile = argv[1]
#outfile = argv[2]

coolers = glob.glob('../DATA/COOL/*')

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
        
    if exp==cell:
         for chrom in chromosomes:
            outfile = '../DATA/TXT/{}.{}.{}.{}.txt'.format('mtx', cell, res, chrom)
            cooler2txt_chr(cool, outfile, fmt='mtx', chromosome=chrom)
            outfile = '../DATA/TXT/{}.{}.{}.{}.txt'.format('sparse_bins', cell, res, chrom)
            cooler2txt_chr(cool, outfile, fmt='sparse_bins', chromosome=chrom)
            outfile = '../DATA/TXT/{}.{}.{}.{}.txt'.format('sparse_coords', cell, res, chrom)
            cooler2txt_chr(cool, outfile, fmt='sparse_coords', chromosome=chrom)
        
    if res=='10' or res=='100':
        for chrom in chromosomes:
            outfile = '../DATA/IMG/{}.{}.{}.{}.png'.format(cell, exp, res, chrom)
            cooler2png_chr(cool, outfile, cmap=cmap, chromosome=chrom)    