#!/bin/bash     
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-100 #parallel cores
#SBATCH -t 12:00:00 # time
#SBATCH -o ../outputs_tradeoffs/outputs_12_02/outfile # send stdout to outfile
#SBATCH -e ../outputs_tradeoffs/outputs_12_02/errfile # send stderr to errfile

if [[ -v SLURM_ARRAY_TASK_ID ]];
then
   echo $SLURM_ARRAY_TASK_ID; simID=$SLURM_ARRAY_TASK_ID
else
   simID=1
fi

<<oldcode
mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2;
mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2/mutdata;

for r in {0..5}; do
    for mk in {0..74};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mk}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script10_exp_2.slim > ../outputs_tradeoffs/outputs_11_16_pair_2/dump.txt;
    done;
done

mkdir -p ../outputs_tradeoffs/outputs_11_25_pair_2;
mkdir -p ../outputs_tradeoffs/outputs_11_25_pair_2/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_11_25_pair_2/mutdata;

for r in {0..9}; do
    for mk in {0..14};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mk}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script10_exp_2new.slim > ../outputs_tradeoffs/outputs_11_25_pair_2/dump.txt;
    done;
done


mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2_extra;
mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2_extra/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_11_16_pair_2_extra/mutdata;

for r in {0..5}; do
    for mk in {75..89};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mk}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script10_exp_2.slim > ../outputs_tradeoffs/outputs_11_16_pair_2_extra/dump.txt;
    done;
done
mks=(18 43 68 17 42 67 16 41 66 23 48 73 22 47 72 21 46 71)
mks=(64 10 17 37 44)


#mkdir -p ../outputs_tradeoffs/outputs_12_02;
#mkdir -p ../outputs_tradeoffs/outputs_12_02/logfiles;
#mkdir -p ../outputs_tradeoffs/outputs_12_02/mutdata;

mks=(65 70 80 85)
for mk in ${!mks[@]};do
    echo "simID: $simID"
    echo "rep: 1"
    echo "mutation kernel: ${mks[$mk]}"
    /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=1 -d mes=$simID -d mk=${mks[$mk]} script10_exp_2_extinct.slim > ../outputs_tradeoffs/outputs_12_02/dump.txt;
done
oldcode

/cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=1 -d mes=$simID  -d mk=14 script10_exp_2_extinct.slim > ../outputs_tradeoffs/outputs_12_02/dump.txt;

