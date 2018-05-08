#!/bin/bash
#PBS -l walltime=100:00:00,mem=4gb,nodes=1:ppn=4
#PBS -t 0-1
#PBS -d.

### Params for SGE:
##$ -S /bin/bash
##$ -V
##$ -cwd
##$ -t 1-7
##$ -pe smp 4
##$ -l mem_free=4G
## 0-124

#set -o errexit
#set -o nounset
#set -o pipefail

# Some specific to platform actions:
unset PYTHONPATH
source activate distiller-editable 

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

GENOME="dm3" 
INDIR="../DATA/FASTQ/"
BAMDIR="../DATA/BAM/"
PAIRDIR="../DATA/PAIR/"
COOLDIR="../DATA/COOL/"
STATSDIR="../DATA/STATS/"
GENOMEDIR="../DATA/GENOME/"

mkdir -p $BAMDIR
mkdir -p $PAIRDIR
mkdir -p $STATSDIR
mkdir -p ../DATA/PAIRIX
mkdir -p ../DATA/COOL

#MYDIR=(`ls $INDIR/[A,B,s,D]*`)
MYDIR=(`ls $INDIR/D*`)

for ((i=${PBS_ARRAYID}; i < ${PBS_ARRAYID}+1; i++))
do
  
  j1=$(($i*2))
  
  dir="${MYDIR[$j1]}"
  pref=${dir##*/}
  pref=${pref%_*_*}

  if [[ $GENOME =~ 'dm3' ]] 
  then 
    file1="${pref}_R1_001.fastq.gz"
    file2="${pref}_R2_001.fastq.gz"
  else
    file1="${pref}_1.fastq.gz"
    file2="${pref}_2.fastq.gz"
  fi
  
  echo $pref $file1 $file2
  #bwa mem -t 4 -v 3 -SP ${GENOMEDIR}/${GENOME}.fa.gz ${INDIR}/$file2 ${INDIR}/$file1 | samtools view -bS > ${BAMDIR}/${pref}.bam

  ##pairsamtools parse ${BAMDIR}/${pref}.bam -c $GENOMEDIR/${GENOME}.reduced.chrom.sizes --output ${PAIRDIR}/${pref}.pairsam -f $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --walks-policy ligation_junctions --drop-sam --add-columns rfrag --output-stats ${STATSDIR}/${pref}.txt

  if [[ $pref =~ 'Dros' ]]
  then
    #pairsamtools parse ${BAMDIR}/${pref}.bam -c $GENOMEDIR/${GENOME}.reduced.chrom.sizes --output ${PAIRDIR}/${pref}.pairsam.full -f $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --walks-policy mask --drop-sam --output-stats ${STATSDIR}/${pref}.txt.full

    pairsamtools sort --nproc 1 -o ${PAIRDIR}/${pref}.pairsam.full.sorted ${PAIRDIR}/${pref}.pairsam.full

    ##pairsamtools restrict ${PAIRDIR}/${pref}.pairsam.full --frags $GENOMEDIR/rfrags_${GENOME}_DpnII.txt --output ${PAIRDIR}/${pref}.pairsam.full.restricted
  fi

done

