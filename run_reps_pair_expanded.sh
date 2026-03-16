#!/bin/bash
#SBATCH -J init_params_4
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-100 #parallel cores
#SBATCH -t 96:00:00 # time
#SBATCH -o ../outputs_tradeoffs/outputs_01_02/outfile # send stdout to outfile
#SBATCH -e ../outputs_tradeoffs/outputs_01_02/errfile # send stderr to errfile

if [[ -v SLURM_ARRAY_TASK_ID ]];
then
   echo $SLURM_ARRAY_TASK_ID; simID=$SLURM_ARRAY_TASK_ID
else
   simID=1
fi

<< oldcode
mkdir -p ../outputs_tradeoffs/outputs_11_28;
mkdir -p ../outputs_tradeoffs/outputs_11_28/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_11_28/mutdata;

for mk1 in {0..26}; do
    for mk2 in {0..26};do
        echo "simID: $simID"
        echo "mutation kernel 1: ${mk1}"
        echo "mutation kernel 2: ${mk2}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d mes=$simID -d mk_a=${mk1} -d mk_b=${mk2} script10_exp_2_expanded.slim > ../outputs_tradeoffs/outputs_11_28/dump.txt;
    done;
done

#mkdir -p ../outputs_tradeoffs/outputs_12_11;
#mkdir -p ../outputs_tradeoffs/outputs_12_11/logfiles;
#mkdir -p ../outputs_tradeoffs/outputs_12_11/mutdata;
#mkdir -p ../outputs_tradeoffs/outputs_12_11/muts_15;
#mkdir -p ../outputs_tradeoffs/outputs_12_11/muts_16;

echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2"  > ../outputs_tradeoffs/outputs_12_11/muts_15/mutation_params_$simID.csv;
#echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2"  > ../outputs_tradeoffs/outputs_12_11/muts_16/mutation_params_$simID.csv;

for rep in {0..99}; do
    echo "rep: ${rep}"
    echo "simID: $simID"
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=15 -d mk_b=15 script10_2_expanded.slim > ../outputs_tradeoffs/outputs_12_11/dump.txt;
    #/cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=16 -d mk_b=16 script10_2_expanded.slim > ../outputs_tradeoffs/outputs_12_11/dump.txt;
done

mkdir -p ../outputs_tradeoffs/outputs_12_14;
mkdir -p ../outputs_tradeoffs/outputs_12_14/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_12_14/mutdata;

echo "rep,mes,mk1,mk2,mk3,mk4,mut,mutr_sp1,mutr_sp2"  > ../outputs_tradeoffs/outputs_12_14/mutation_params_$simID.csv;

for rep in {0..199}; do
    echo "rep: ${rep}"
    echo "simID: $simID"
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=15 -d mk_b=15 script10_2_expanded_edit.slim > ../outputs_tradeoffs/outputs_12_14/dump.txt;
done

mkdir -p ../outputs_tradeoffs/outputs_12_18;
mkdir -p ../outputs_tradeoffs/outputs_12_18/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_12_18/mutdata;
mkdir -p ../outputs_tradeoffs/outputs_12_18/mutdata;
mkdir -p ../outputs_tradeoffs/outputs_12_18/muts_15;
mkdir -p ../outputs_tradeoffs/outputs_12_18/muts_16;

echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2"  > ../outputs_tradeoffs/outputs_12_18/muts_15/mutation_params_$simID.csv;
echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2"  > ../outputs_tradeoffs/outputs_12_18/muts_16/mutation_params_$simID.csv;

for rep in {0..99}; do
    echo "rep: ${rep}"
    echo "simID: $simID"
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=15 -d mk_b=15 script10_2_expanded.slim > ../outputs_tradeoffs/outputs_12_11/dump.txt;
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=16 -d mk_b=16 script10_2_expanded.slim > ../outputs_tradeoffs/outputs_12_11/dump.txt;
done

mkdir -p ../outputs_tradeoffs/outputs_12_22;
mkdir -p ../outputs_tradeoffs/outputs_12_22/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_12_22/mutdata;

echo "rep,mes,mk1,mk2,mk3,mk4,mut,mutr_sp1,mutr_sp2"  > ../outputs_tradeoffs/outputs_12_22/mutation_params_$simID.csv;

for rep in {0..99}; do
    echo "rep: ${rep}"
    echo "simID: $simID"
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=15 -d mk_b=15 script10_2_expanded_edit.slim > ../outputs_tradeoffs/outputs_12_22/dump.txt;
done

mkdir -p ../outputs_tradeoffs/outputs_12_23_3;
mkdir -p ../outputs_tradeoffs/outputs_12_23_3/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_12_23_3/mutdata;

for mk1 in {0..26}; do
    for mk2 in {0..26};do
        echo "simID: $simID"
        echo "mutation kernel 1: ${mk1}"
        echo "mutation kernel 2: ${mk2}"
        /cluster/tufts/uricchiolab/software/build/slim -s $simID -d mes=$simID -d mk_a=${mk1} -d mk_b=${mk2} script10_exp_2_expanded3.slim > ../outputs_tradeoffs/outputs_12_23_3/dump.txt;
    done;
done

oldcode
mkdir -p ../outputs_tradeoffs/outputs_01_02;
mkdir -p ../outputs_tradeoffs/outputs_01_02/logfiles;
mkdir -p ../outputs_tradeoffs/outputs_01_02/mutdata;

for rep in {50..99}; do
    echo "rep: ${rep}"
    echo "simID: $simID"
    /cluster/tufts/uricchiolab/software/build/slim -d mes=$simID -d rep=${rep} -d mk_a=1 -d mk_b=2 script10_3.slim > ../outputs_tradeoffs/outputs_01_02/dump.txt;
done

