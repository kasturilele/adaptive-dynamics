#!/bin/bash     
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-100 #parallel cores
#SBATCH -t 96:00:00 # time
#SBATCH -o ../outputs_tradeoffs/outputs_02_03_3/outfile # send stdout to outfile
#SBATCH -e ../outputs_tradeoffs/outputs_02_03_3/errfile # send stderr to errfile

if [[ -v SLURM_ARRAY_TASK_ID ]];
then
   echo $SLURM_ARRAY_TASK_ID; simID=$SLURM_ARRAY_TASK_ID
else
   simID=1
fi

#mkdir -p ../outputs_tradeoffs/outputs_11_16_3;
#mkdir -p ../outputs_tradeoffs/outputs_11_16_3/logfiles;
#mkdir -p ../outputs_tradeoffs/outputs_11_16_3/mutdata;

<< oldcode
for r in {0..4}; do
    for mk in {0..74};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mk}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script9_exp_3.slim > ../outputs_tradeoffs/outputs_11_16_3/dump.txt;
    done;
done



for r in {0..4}; do
    for mk in {75..89};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mk}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script9_exp_3.slim > ../outputs_tradeoffs/outputs_11_16_3/dump.txt;
        #/cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mk} script9_exp_4.slim > ../outputs_tradeoffs/outputs_11_16_4/dump.txt;
    done;
done



mkdir -p ../outputs_tradeoffs/outputs_02_03_3;
mkdir -p ../outputs_tradeoffs/outputs_02_03_3/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_02_03_3/mutdata;

mks=(10 11 12 13 14 35 36 37 38 39 60 61 62 63 64 15 16 17 18 19 40 41 42 43 44 65 66 67 68 69 20 21 22 23 24 45 46 47 48 49 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89)
for r in {0..6}; do
    for mk in ${!mks[@]};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mks[$mk]}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mks[$mk]} script9_exp_3new.slim > ../outputs_tradeoffs/outputs_02_03_3/dump.txt;
    done;
done

oldcode

mks=(90 91 92)
for r in {0..6}; do
    for mk in ${!mks[@]};do
        echo "simID: $simID"
        echo "rep: ${r}"
        echo "mutation kernel: ${mks[$mk]}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d rep=${r} -d mes=$simID -d mk=${mks[$mk]} script9_exp_3new.slim > ../outputs_tradeoffs/outputs_02_03_3/dump.txt;
    done;
done
