#!/bin/bash

input=$1
k=$2
outdir=$3/
#anchor=$3

TMP=${outdir}tmp
#COUNT_DIR=${outdir}kmc_count
#ID_DIR=${outdir}kmc_id
#BITMAP_DIR=${outdir}kmc_bitvec

mkdir -p ${TMP}

#BIN_D=/home/skovaka/code/pan_kmers
#KMC_DIR="/grid/schatz/home/jenike/skovaka/KMC"
KMC_DIR="/home/skovaka/code/pan_kmers/KMC"

source ${KMC_DIR}/py_kmc_api/set_path.sh

#input is the file with the list of files 
#input="genomes.tsv"
#k=21

i=0
j=0
samp_count=0
db_count=1
while read name fasta
do

  if (($i >= 32)); then
        i=0
        j=$((j + 1))
        db_count=$((db_count+1))
  fi

  let "id = 2**i"
  echo "Counting ${name}"  

  ##Run kmc for each file
  /usr/bin/time -v ${KMC_DIR}/bin/kmc -k${k} -t12 -m64 -ci1 -cs10000000 -fm ${fasta} ${outdir}${name}.cnt ${TMP} > ${outdir}${name}.cnt.err 2>&1

  #Modify the kmc database 
  /usr/bin/time -v ${KMC_DIR}/bin/kmc_tools -t12 transform ${outdir}${name}.cnt set_counts ${id} ${outdir}${name}.id > ${outdir}${name}.id.err 2>&1

  i=$((i + 1))
  samp_count=$((samp_count + 1))

done < $input

echo "Merging"

dbs=""

for i in `seq 0 $j`; do

    let "h = ($i+1)*32"

    if (($db_count == 1 || $i < $j)); then
        t=32
    else
        let "t = $samp_count-32"
    fi
    echo $n

    echo "INPUT:" > ${outdir}opdef.${i}.txt
    head -n $h ${input} | tail -n $t | awk -v out=$outdir '{print $1" = "out$1".id"}' >> ${outdir}opdef.${i}.txt
    echo "OUTPUT:" >> ${outdir}opdef.${i}.txt
    head -n $h ${input} | tail -n $t | awk -v i=$i -v out=$outdir 'NR==1 {s=$1} NR>1 {s=s" + "$1} END {print out"bitvecs."i" = "s}' >> ${outdir}opdef.${i}.txt
    echo "-ocsum" >> ${outdir}opdef.${i}.txt

    #Add all of the kmc databases together
    ${KMC_DIR}/bin/kmc_tools complex ${outdir}opdef.${i}.txt

    dbs="$dbs${outdir}bitvecs.$i "

done

echo $dbs

#Run the "get_counts" kmc function 
panagram anchor ${input} ${dbs} -o ${outdir}anchor
bgzip -r ${outdir}anchor.pank -I ${outdir}anchor.panx
bgzip -r ${outdir}anchor.100nt.pank -I ${outdir}anchor.100nt.panx
