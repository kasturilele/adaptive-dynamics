#!/bin/bash
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-1
#SBATCH -t 20:00:00 # time
#SBATCH -o ../scripts/outfile # send stdout to outfile
#SBATCH -e ../scripts/errfile # send stderr to errfile

module load anaconda/2024.10

echo "mutation_kernel,rep,mes,tick,num_individuals_species1,r1,a11"  > ../outputs_tradeoffs/log_combined_single3_extra.csv;
#echo "rep,mes,mutation_kernel,time_origin,time_fixed,effect_size,mutation_type,blank,blank"  > ../outputs_tradeoffs/fixed_combined_single3.csv;
#echo "mutation_kernel,rep,mes,tick,num_individuals_species1,r1,a11"  > ../outputs_tradeoffs/log_combined_single4.csv;
#echo "rep,mes,mutation_kernel,time_origin,time_fixed,effect_size,mutation_type,blank,blank"  > ../outputs_tradeoffs/fixed_combined_single4.csv;

#mks=(10 11 12 13 14 35 36 37 38 39 60 61 62 63 64 15 16 17 18 19 40 41 42 43 44 65 66 67 68 69 20 21 22 23 24 45 46 47 48 49 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89)
mks=(35 40 45 60 65 70 80 90 91 92)

for rep in {0..6}; do 
    for mes in {1..100}; do
        for mk in ${!mks[@]};do
        cat ../outputs_tradeoffs/outputs_02_03_3/logfiles/output_${rep}_${mes}_${mks[$mk]}.csv | awk -v rep="$rep" -v mes="$mes" -v mk="${mks[$mk]}" '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_single3_extra.csv;
        #cat ../outputs_tradeoffs/outputs_11_16_3/mutdata/fixed_${rep}_${mes}_${mks[$mk]}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/fixed_combined_single3.csv;
        #cat ../outputs_tradeoffs/outputs_11_16_4/logfiles/output_${rep}_${mes}_${mks[$mk]}.csv | awk -v rep="$rep" -v mes="$mes" -v mk="${mks[$mk]}" '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_single4.csv;
        #cat ../outputs_tradeoffs/outputs_11_16_4/mutdata/fixed_${rep}_${mes}_${mks[$mk]}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/fixed_combined_single4.csv;

        done;
    done;
done;