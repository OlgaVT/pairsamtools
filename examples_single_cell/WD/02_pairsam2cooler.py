"""
Converts pset of pairsam files for one cell into cooler files. 
TODO: optimize experiments handling.
Example run in bash: 
 python pairsam2cooler.py "../DATA/PAIR/${pref}*.pairsam" "$pref" "../DATA/PAIRIX/${pref}.pairix" "../DATA/COOL/${pref}.{}.cool" "../DATA/STATS/${pref}.filter_stats"
"""

from sys import argv
import glob
from basic_utils import *

mask = argv[1]
cell = argv[2]
out_pairix = argv[3] # "../DATA/PAIRIX/A9_S38.pairsam"
out_cool_mask = argv[4] # "../DATA/COOL/A9_S38.{}.cool"
out_stats = argv[5]

print(mask, cell, out_pairix, out_cool_mask)

filelist = glob.glob(mask)
exp_list = [x.split('/')[-1].split('.')[0] for x in filelist]

print(filelist)
print(exp_list)

df = read_pairsams(filelist, exp_list, cell)

df_filtered, stats = filter_pair_df(df)
stats['cell'] = cell
stats['exp'] = cell
with open(out_stats, 'w') as outf:
    for k in sorted(stats.keys()):
        outf.write('{}\t{}\n'.format(k, stats[k]))
        
resolutions = [100, 20, 10, 1]

create_cooler(df_filtered, out_pairix, out_cool_mask, 
              "../DATA/GENOME/dm3.reduced.chrom.sizes", resolutions_list=resolutions)

for exp in exp_list:
    df_tmp = df.query("exp=='{}'".format(exp))
    df_tmp_filtered, stats = filter_pair_df(df_tmp)
    
    stats['cell'] = cell
    stats['exp'] = exp
    with open(out_stats+"."+exp, 'w') as outf:
        for k in sorted(stats.keys()):
            outf.write('{}\t{}\n'.format(k, stats[k]))
            
    create_cooler(df_tmp_filtered, out_pairix, out_cool_mask+"."+exp, 
              "../DATA/GENOME/dm3.reduced.chrom.sizes", resolutions_list=resolutions)
