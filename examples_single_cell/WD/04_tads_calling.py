"""
Python script that calculates the TADs  with parameters optimizing the expected TADs size.

Usage:
python 04_tads_calling.py <cooler input file> <prefix for output with TADs> <method name: modularity, armatus> <retrieve matrix: balanced, non-balanced> <expected TADs size> <max size of interTAD>

Example usage:
python 04_tads_calling.py ../DATA/COOL/A6.10.cool ../DATA/TAD/tmp modularity non-balanced 12 3

Output files:
<prefix>.bed bed file with end coordinate -- 0.5*resolution to produce visible results in juicebox with this file
<prefix>.2Dannot -- 2Dannotation file for juicebox
<prefix>.opt_gamma.pickle -- pickle file with dictionary with optimal gammas
<prefix>.opt_segmentation.pickle-- pickle file with dictionary with optimal segmentation

"""

# libraries import
from basic_utils import *
import numpy as np
import pickle

from sys import argv

# Parameters setting
input_cooler      = argv[1]
output_prefix     = argv[2]
calling_algorithm = argv[3] # armatus of modularity
reading_mode      = argv[4] # balanced or not
tad_size          = int(argv[5]) # median expected TADs size in bins
intertad_size     = int(argv[6]) # max length for segmentation unit to be considered as interTAD
max_tad_size      = 10000 #int(argv[7]) # min length for segmentation unit to be considered as interTAD (too large TAD)

output_bed        = output_prefix + '.bed'
output_2Dannot    = output_prefix + '.2Dannot'
output_opt_gammas = output_prefix + '.opt_gamma.pickle'
output_opt_segmentation = output_prefix + '.opt_segmentation.pickle'

if reading_mode=='balanced':
    balance = True
else:
    balance = False
    
# Reading input file
c = cooler.Cooler(input_cooler)

segmentations = {}
opt_gammas    = {}

# Iteration over chromosomes
for chrom in c.chromnames:

    if chrom=='chrM':
        continue

    mtx = c.matrix(balance=balance).fetch('{0}'.format(chrom)).astype(float)
    mtx[np.isnan(mtx)] = 0
    np.fill_diagonal(mtx, 0)

    # Algorithm-specific parameters setup and calculations
    if calling_algorithm == 'modularity':
        step1 = 1
        step2 = 0.001
        mx1 = 100
        mx2 = 1
    elif calling_algorithm == 'armatus':
        step1 = 0.1
        step2 = 0.001
        mx1 = 5
        mx2 = 0.2
        
        # log of matrix, scaling to the positive numbers only (required by lavaburst scoring)
        mn = np.percentile(mtx[mtx>0], 1)
        mx = np.percentile(mtx[mtx>0], 99)

        mtx[mtx<=mn] = mn
        mtx[mtx>=mx] = mx
        
        mtx = np.log(mtx)
        mtx = mtx-np.min(mtx)

    else:
        logging.error("Algorithm {} not known!".format(calling_algorithm))

    # Iterative search of optimal gammas
    gammas = np.arange(step1, mx1, step1)
    #opt_mean, opt_gamma  = find_optimal_gamma(mtx, tad_size, gammas=gammas, method=calling_algorithm, max_intertad_size=intertad_size, epsilon=0.1, max_tad_size=max_tad_size, min_tads_number=0.5*len(mtx)/tad_size)
    opt_mean, opt_gamma, opt_ntads  = find_optimal_gamma_2opt(mtx, tad_size, gammas=gammas, 
                                              method=calling_algorithm, max_intertad_size=intertad_size, 
                                              max_tad_size=max_tad_size)
    
    gammas = np.arange(max(opt_gamma-mx2,step2), opt_gamma+mx2, step2)
    #opt_mean, opt_gamma  = find_optimal_gamma(mtx, tad_size, gammas=gammas, method=calling_algorithm, max_intertad_size=intertad_size, epsilon=0.1, max_tad_size=max_tad_size, min_tads_number=0.5*len(mtx)/tad_size)
    opt_mean, opt_gamma, opt_ntads  = find_optimal_gamma_2opt(mtx, tad_size, gammas=gammas, 
                                              method=calling_algorithm, max_intertad_size=intertad_size, 
                                              max_tad_size=max_tad_size)

    segments = produce_segmentation(mtx, opt_gamma, method=calling_algorithm, 
                                    max_intertad_size=intertad_size, max_tad_size=max_tad_size)
    
    logging.info( "{}:\t{}\t{}\t{}\t{}\t{}".format(chrom, len(segments[:,1]), np.max(segments[:,1]-segments[:,0]), np.median(segments[:,1]-segments[:,0]), np.mean(segments[:,1]-segments[:,0]), opt_gamma) )

    segmentations[chrom] = segments.copy()
    opt_gammas[chrom]    = opt_gamma


pickle.dump(opt_gammas,    open(output_opt_gammas,       'wb'))
pickle.dump(segmentations, open(output_opt_segmentation, 'wb'))

segmentations_to_bed(segmentations,     output_bed,     c.binsize)
segmentations_to_2Djuice(segmentations, output_2Dannot, c.binsize)