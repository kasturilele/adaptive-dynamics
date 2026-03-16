#!/bin/bash
#SBATCH -N 1 # nodes requested
#SBATCH -n 1 # tasks requested
#SBATCH -c 4 # cores requested
#SBATCH --mem=8000 # memory in Mb
#SBATCH --array=1-1
#SBATCH -t 20:00:00 # time
#SBATCH -o ../scripts/outfile2 # send stdout to outfile
#SBATCH -e ../scripts/errfile2 # send stderr to errfile

module load anaconda/2024.10

#echo "rep,mes,mutation_kernel,species,time_origin,time_fixed,effect_size,mutation_type,blank,blank"  > ../outputs_tradeoffs/fixed_combined_cox_ext2.csv;

<<oldcode

echo "rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_exp_15_ext2.csv;
echo "rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_exp_16_ext2.csv;
echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2" > ../outputs_tradeoffs/mutkern_15_ext2.csv;
echo "rep,mes,mk1,mk2,mk3,mk4,mutr,muta,mut_sp1,mut_sp2" > ../outputs_tradeoffs/mutkern_16_ext2.csv;


for rep in {0..5}; do 
    for mes in {1..100}; do
        for mk in {0..74}; do
        cat ../outputs_tradeoffs/outputs_11_16_pair_2/logfiles/10_output_${rep}_${mes}_${mk}.csv | awk -v rep="$rep" -v mes="$mes" -v mk="$mk" '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_pair2new.csv;
        cat ../outputs_tradeoffs/outputs_11_16_pair_2/mutdata/10_fixed_${rep}_${mes}_${mk}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/fixed_combined_pair2new.csv;

        done;
    done;
done;

for rep in {0..5}; do 
    for mes in {1..100}; do
        for mk in {75..89}; do
        cat ../outputs_tradeoffs/outputs_11_16_pair_2_extra/logfiles/10_output_${rep}_${mes}_${mk}.csv | awk -v rep="$rep" -v mes="$mes" -v mk="$mk" '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_pair2new.csv;
        cat ../outputs_tradeoffs/outputs_11_16_pair_2_extra/mutdata/10_fixed_${rep}_${mes}_${mk}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/fixed_combined_pair2new.csv;

        done;
    done;
done;
mks=(18 43 68 17 42 67 16 41 66 23 48 73 22 47 72 21 46 71)

for mes in {1..100}; do
    for mk in ${!mks[@]};do
    cat ../outputs_tradeoffs/outputs_12_02/logfiles/10_output_1_${mes}_${mks[$mk]}.csv | awk -v rep=1 -v mes="$mes" -v mk=${mks[$mk]} '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_extinct.csv;

    done;
done;
mks=(64 10 17 37 44)
mks=(64 10 17 37 44 14)
mks=(65 70 80 85)

for mes in {1..100}; do
    for mk in ${!mks[@]};do
    cat ../outputs_tradeoffs/outputs_12_02/logfiles/10_output_1_${mes}_${mks[$mk]}.csv | awk -v rep=1 -v mes="$mes" -v mk=${mks[$mk]} '{if (NR > 1) print mk, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_combined_cox_ext2.csv;
    cat ../outputs_tradeoffs/outputs_12_02/mutdata/10_fixed_1_${mes}_${mks[$mk]}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/fixed_combined_cox_ext2.csv;

    done;
done;

for mk1 in {0..26}; do
    for mk2 in {0..26};do
        for mes in {1..100}; do
            cat ../outputs_tradeoffs/outputs_11_28/logfiles/10_output_1_${mes}_${mk1}_${mk2}.csv | awk -v rep=1 -v mes="$mes" -v mk1="$mk1" -v mk2="$mk2" '{if (NR > 1) print mk1, ",", mk2, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_expand.csv;
        done;
    done;
done

for mes in {1..100}; do
    cat ../outputs_tradeoffs/outputs_11_28/logfiles/10_output_1_${mes}_24_13.csv | awk -v rep=1 -v mes="$mes" -v mk1=24 -v mk2=13 '{if (NR > 1) print mk1, ",", mk2, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_exp_temp.csv;
done;

for mes in {1..200}; do 
    for rep in {0..99}; do
        cat ../outputs_tradeoffs/outputs_12_11/logfiles/10_output_${rep}_${mes}_15_15.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_exp_15_ext2.csv;
        cat ../outputs_tradeoffs/outputs_12_11/logfiles/10_output_${rep}_${mes}_16_16.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_exp_16_ext2.csv;
    done;
    cat ../outputs_tradeoffs/outputs_12_11/muts_15/mutation_params_${mes}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/mutkern_15_ext2.csv;
    cat ../outputs_tradeoffs/outputs_12_11/muts_16/mutation_params_${mes}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/mutkern_16_ext2.csv;
done;

for mes in {201..500}; do 
    for rep in {0..99}; do
        cat ../outputs_tradeoffs/outputs_12_18/logfiles/10_output_${rep}_${mes}_15_15.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_exp_15_ext2.csv;
        cat ../outputs_tradeoffs/outputs_12_18/logfiles/10_output_${rep}_${mes}_16_16.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_exp_16_ext2.csv;
    done;
    cat ../outputs_tradeoffs/outputs_12_18/muts_15/mutation_params_${mes}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/mutkern_15_ext2.csv;
    cat ../outputs_tradeoffs/outputs_12_18/muts_16/mutation_params_${mes}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/mutkern_16_ext2.csv;
    
done;

for mes in {1..200}; do 
    for rep in {0..99}; do
        cat ../outputs_tradeoffs/outputs_12_22/logfiles/10_output_${rep}_${mes}_15_15.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_me_new.csv;
    done;
    cat ../outputs_tradeoffs/outputs_12_22/mutation_params_${mes}.csv | awk '{if (NR > 1) print}' >> ../outputs_tradeoffs/mutkern_me_new.csv;
done;

echo "mk1,mk2,rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_init0.csv;
echo "mk1,mk2,rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_init2.csv;
echo "mk1,mk2,rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_init3.csv;

for mk1 in {0..26}; do
    for mk2 in {0..26};do
        for mes in {1..100}; do
            cat ../outputs_tradeoffs/outputs_12_23_0/logfiles/10_output_0_${mes}_${mk1}_${mk2}.csv | awk -v rep=1 -v mes="$mes" -v mk1="$mk1" -v mk2="$mk2" '{if (NR > 1) print mk1, ",", mk2, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_init0.csv;
            cat ../outputs_tradeoffs/outputs_12_23_2/logfiles/10_output_2_${mes}_${mk1}_${mk2}.csv | awk -v rep=1 -v mes="$mes" -v mk1="$mk1" -v mk2="$mk2" '{if (NR > 1) print mk1, ",", mk2, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_init2.csv;
            cat ../outputs_tradeoffs/outputs_12_23_3/logfiles/10_output_3_${mes}_${mk1}_${mk2}.csv | awk -v rep=1 -v mes="$mes" -v mk1="$mk1" -v mk2="$mk2" '{if (NR > 1) print mk1, ",", mk2, ",", rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_init3.csv;
        done;
    done;
done

oldcode

echo "rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_init_00.csv;
echo "rep,mes,tick,Strain_1,Strain_2,num_individuals_species1,num_individuals_species2,r1,a11,r2,a22"  > ../outputs_tradeoffs/log_comb_init_12.csv;

for mes in {1..100}; do 
    for rep in {0..99}; do
        cat ../outputs_tradeoffs/outputs_01_02/logfiles/10_output_${rep}_${mes}_0_0.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_init_00.csv;
        cat ../outputs_tradeoffs/outputs_01_02/logfiles/10_output_${rep}_${mes}_1_2.csv | awk -v rep="$rep" -v mes="$mes" '{if (NR > 1) print rep,",", mes,",", $0}' >> ../outputs_tradeoffs/log_comb_init_12.csv;
    done;
done;

