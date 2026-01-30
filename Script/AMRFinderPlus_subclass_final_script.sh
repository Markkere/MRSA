#!/bin/bash

dir="$1" #Fasta file dir
output_dir="$2" #desired output dir
seq_input="$3" #Illumina or Nanopore
species="$4" #Desired species name according to AMRFinderPlus available species
seq_input_lower=$(echo $seq_input | tr '[:upper:]' '[:lower:]')
##Take 3 input antibiotics
ant_1="$5"
ant_2="$6"
ant_3="$7"

#Check for 7 argument inputs
if [ "$#" -lt 7 ]; then
 echo "7 arguments are not sastisfied"
 exit 1
fi

if [ -z "$4" ]; then
 echo "Please enter your desired species"
 exit 1
fi

if [ -z "$5" ]; then
 echo "Please enter your desired antibiotic"
 exit 1
fi

#Convert to lower case
ant_1_lower=$(echo $ant_1 | tr '[:upper:]' '[:lower:]')
ant_2_lower=$(echo $ant_2 | tr '[:upper:]' '[:lower:]')

#for case where user want to look at only 1 antibiotic
if [[ "$ant_2_lower" == 'na' ]]; then
   ant_2="NA_1"
fi
ant_3_lower=$(echo $ant_3 | tr '[:upper:]' '[:lower:]')
if [[ "$ant_3_lower" == 'na' ]]; then
   ant_3="NA_2"	
fi

for file in "$dir"/*.fna;do
 name=$(basename "$file")
 ending=".csv"
 file_name="$name$ending"
 amrfinder -n "$file" -O "$species" --plus > "$output_dir/$file_name"
done

result_dir="$2"
current_dir=$(pwd)
#Final name
final_name="AMRFinderPlus_${seq_input_lower}_assembly_result.tsv" 

##templates
tsv=".tsv"
exp_ant1_gene_t="exp_ant1_gene"
exp_ant2_gene_t="exp_ant2_gene"
exp_ant3_gene_t="exp_ant3_gene"
exp_ant1_seq_t="exp_ant1_seq"
exp_ant2_seq_t="exp_ant2_seq"
exp_ant3_seq_t="exp_ant3_seq"
exp_ant1_iden_t="exp_ant1_iden"
exp_ant2_iden_t="exp_ant2_iden"
exp_ant3_iden_t="exp_ant3_iden"
exp_ant1_closest_seq_t="exp_ant1_clst_seq"
exp_ant2_closest_seq_t="exp_ant2_clst_seq"
exp_ant3_closest_seq_t="exp_ant3_clst_seq"
exp_ant1_subclass_t="exp_ant1_subclass"
exp_ant2_subclass_t="exp_ant2_subclass"
exp_ant3_subclass_t="exp_ant3_subclass"
pred='_pred'
gene='_gene_symbol'
seq='_seqeunce_name'
subclass='_subclass'
clst_seq='_closest_seq'
percent_iden='_percent_identity_to_reference'
exp_stress_gene="_stress_gene"
exp_virulence_gene="_virulence_gene"
exp_stress_seq="_stress_seq"
exp_virulence_seq="_virulence_seq"
echo -e 'genome_id\t'"$ant_1$pred"'\t'"$ant_1$gene"'\t'"$ant_1$seq"'\t'"$ant_1$subclass"'\t'"$ant_1$clst_seq"'\t'"$ant_1$percent_iden"'\t'"$ant_2$pred"'\t'"$ant_2$gene"'\t'"$ant_2$seq"'\t'"$ant_2$subclass"'\t'"$ant_2$clst_seq"'\t'"$ant_2$percent_iden"'\t'"$ant_3$pred"'\t'"$ant_3$gene"'\t'"$ant_3$seq"'\t'"$ant_3$subclass"'\t'"$ant_3$clst_seq"'\t'"$ant_3$percent_iden"'\t'"$species$exp_stress_gene"'\t'"$species$exp_virulence_gene"'\t'"$species$exp_stress_seq"'\t'"$species$exp_virulence_seq"> "$current_dir/header.tsv"

#Loop through all result files  
for file in "$result_dir"/*.csv; do
 echo "$file"
 genome_id="$(basename "$file")" #extract genome id
 exp_ant1_gene="$(awk -F '\t' -v ab="$ant_1_lower" 'tolower($12) ~ ab' "$file" | cut -f6 | sort | uniq | paste -sd ",")" #Look for desired antibiotic in the subclass column and extrating just important columns
 exp_ant2_gene="$(awk -F '\t' -v ab="$ant_2_lower" 'tolower($12) ~ ab' "$file" | cut -f6 | sort | uniq | paste -sd ",")"
 exp_ant3_gene="$(awk -F '\t' -v ab="$ant_3_lower" 'tolower($12) ~ ab' "$file" | cut -f6 | sort | uniq | paste -sd ",")"
 exp_ant1_seq="$(awk -F '\t' -v ab="$ant_1_lower" 'tolower($12) ~ ab' "$file" | cut -f7 | sort | uniq | paste -sd ",")"
 exp_ant2_seq="$(awk -F '\t' -v ab="$ant_2_lower" 'tolower($12) ~ ab' "$file" | cut -f7 | sort | uniq | paste -sd ",")"
 exp_ant3_seq="$(awk -F '\t' -v ab="$ant_3_lower" 'tolower($12) ~ ab' "$file" | cut -f7 | sort | uniq | paste -sd ",")"
 ant1_iden="$(awk -F '\t' -v ab="$ant_1_lower" 'tolower($12) ~ ab' "$file" | cut -f17 | sort | uniq | paste -sd ",")"
 ant2_iden="$(awk -F '\t' -v ab="$ant_2_lower" 'tolower($12) ~ ab' "$file" | cut -f17 | sort | uniq | paste -sd ",")"
 ant3_iden="$(awk -F '\t' -v ab="$ant_3_lower" 'tolower($12) ~ ab' "$file" | cut -f17 | sort | uniq | paste -sd ",")"
 ant1_closest_seq="$(awk -F '\t' -v ab="$ant_1_lower" 'tolower($12) ~ ab' "$file" | cut -f20 | sort | uniq | paste -sd ",")"
 ant2_closest_seq="$(awk -F '\t' -v ab="$ant_2_lower" 'tolower($12) ~ ab' "$file" | cut -f20 | sort | uniq | paste -sd ",")"
 ant3_closest_seq="$(awk -F '\t' -v ab="$ant_3_lower" 'tolower($12) ~ ab' "$file" | cut -f20 | sort | uniq | paste -sd ",")"
 ant1_subclass="$(awk -F '\t' -v ab="$ant_1_lower" 'tolower($12) ~ ab' "$file" | cut -f12 | sort | uniq | paste -sd ",")"
 ant2_subclass="$(awk -F '\t' -v ab="$ant_2_lower" 'tolower($12) ~ ab' "$file" | cut -f12 | sort | uniq | paste -sd ",")"
 ant3_subclass="$(awk -F '\t' -v ab="$ant_3_lower" 'tolower($12) ~ ab' "$file" | cut -f12 | sort | uniq | paste -sd ",")"
 stress_gene="$(grep "STRESS" "$file" | cut -f6 | sort | uniq | paste -sd ",")"
 virulence_gene="$(grep "VIRULENCE" "$file" | cut -f6 | sort | uniq | paste -sd ",")"
 stress_seq="$(grep "STRESS" "$file" | cut -f7 | sort | uniq | paste -sd ",")"
 virulence_seq="$(grep "VIRULENCE" "$file" | cut -f7 | sort | uniq | paste -sd ",")"
 echo "$genome_id" >> "genome_id.tsv"
 if [ -n "$ant1_subclass" ] && [ -n "$ant2_subclass" ] && [ -n "$ant3_subclass" ]; then #assign accordingly to each case
   echo 'R' >> "$current_dir/$ant_1$tsv"
   echo 'R' >> "$current_dir/$ant_2$tsv"
   echo 'R' >> "$current_dir/$ant_3$tsv"
   echo "$exp_ant1_gene" >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo "$exp_ant2_gene" >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo "$exp_ant3_gene" >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo "$exp_ant1_seq" >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo "$exp_ant2_seq" >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo "$exp_ant3_seq" >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo "$ant1_iden" >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo "$ant2_iden" >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo "$ant3_iden" >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo "$ant1_closest_seq" >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo "$ant2_closest_seq" >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo "$ant3_closest_seq" >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo "$ant1_subclass" >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo "$ant2_subclass" >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo "$ant3_subclass" >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant1_subclass" ] && [ -n "$ant2_subclass" ]; then
   echo 'R' >> "$current_dir/$ant_1$tsv"
   echo 'R' >> "$current_dir/$ant_2$tsv"
   echo 'S' >> "$current_dir/$ant_3$tsv"
   echo "$exp_ant1_gene" >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo "$exp_ant2_gene" >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo "$exp_ant1_seq" >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo "$exp_ant2_seq" >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo "$ant1_iden" >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo "$ant2_iden" >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo "$ant1_closest_seq" >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo "$ant2_closest_seq" >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo "$ant1_subclass" >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo "$ant2_subclass" >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant1_subclass" ] && [ -n "$ant3_subclass" ]; then
   echo 'R' >> "$current_dir/$ant_1$tsv"
   echo 'S' >> "$current_dir/$ant_2$tsv"
   echo 'R' >> "$current_dir/$ant_3$tsv"
   echo "$exp_ant1_gene" >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo "$exp_ant3_gene" >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo "$exp_ant1_seq" >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_seq_t$tsv" 
   echo "$exp_ant3_seq" >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo "$ant1_iden" >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo "$ant3_iden" >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo "$ant1_closest_seq" >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo "$ant3_closest_seq" >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo "$ant1_subclass" >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo "$ant3_subclass" >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant2_subclass" ] && [ -n "$ant3_subclass" ]; then
   echo 'S' >> "$current_dir/$ant_1$tsv"
   echo 'R' >> "$current_dir/$ant_2$tsv"
   echo 'R' >> "$current_dir/$ant_3$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo "$exp_ant2_gene" >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo "$exp_ant3_gene" >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo "$exp_ant2_seq" >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo "$exp_ant3_seq" >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo "$ant2_iden" >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo "$ant3_iden" >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo "$ant2_closest_seq" >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo "$ant3_closest_seq" >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo "$ant2_subclass" >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo "$ant3_subclass" >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant1_subclass" ];then
   echo 'R' >> "$current_dir/$ant_1$tsv"
   echo 'S' >> "$current_dir/$ant_2$tsv"
   echo 'S' >> "$current_dir/$ant_3$tsv"
   echo "$exp_ant1_gene" >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo "$exp_ant1_seq" >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo "$ant1_iden" >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo "$ant1_closest_seq" >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo "$ant1_subclass" >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant2_subclass" ];then
   echo 'S' >> "$current_dir/$ant_1$tsv"
   echo 'R' >> "$current_dir/$ant_2$tsv"
   echo 'S' >> "$current_dir/$ant_3$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo "$exp_ant2_gene" >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo "$exp_ant2_seq"  >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo "$ant2_iden" >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo "$ant2_closest_seq"  >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo "$ant2_subclass" >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 elif [ -n "$ant3_subclass" ];then
   echo 'S' >> "$current_dir/$ant_1$tsv"
   echo 'S' >> "$current_dir/$ant_2$tsv"
   echo 'R' >> "$current_dir/$ant_3$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo "$exp_ant3_gene" >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo "$exp_ant3_seq" >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo "$ant3_iden" >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo "$ant3_closest_seq"  >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo "$ant3_subclass" >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 else
   echo 'S' >> "$current_dir/$ant_1$tsv"
   echo 'S' >> "$current_dir/$ant_2$tsv"
   echo 'S' >> "$current_dir/$ant_3$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_gene_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_iden_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant2_closest_seq_t$tsv"
   echo 'NA'  >> "$current_dir/$exp_ant3_closest_seq_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant1_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant2_subclass_t$tsv"
   echo 'NA' >> "$current_dir/$exp_ant3_subclass_t$tsv"
   echo "$stress_gene" >> "$current_dir/$exp_stress_gene$tsv"
   echo "$virulence_gene" >> "$current_dir/$exp_virulence_gene$tsv"
   echo "$stress_seq" >> "$current_dir/$exp_stress_seq$tsv"
   echo "$virulence_seq" >> "$current_dir/$exp_virulence_seq$tsv"
 fi
done

#Paste every intermediate files
paste "$current_dir/genome_id.tsv" "$current_dir/$ant_1$tsv" "$current_dir/$exp_ant1_gene_t$tsv" "$current_dir/$exp_ant1_seq_t$tsv" "$current_dir/$exp_ant1_subclass_t$tsv" "$current_dir/$exp_ant1_closest_seq_t$tsv" "$current_dir/$exp_ant1_iden_t$tsv" "$current_dir/$ant_2$tsv" "$current_dir/$exp_ant2_gene_t$tsv" "$current_dir/$exp_ant2_seq_t$tsv" "$current_dir/$exp_ant2_subclass_t$tsv" "$current_dir/$exp_ant2_closest_seq_t$tsv" "$current_dir/$exp_ant2_iden_t$tsv" "$current_dir/$ant_3$tsv" "$current_dir/$exp_ant3_gene_t$tsv" "$current_dir/$exp_ant3_seq_t$tsv" "$current_dir/$exp_ant3_subclass_t$tsv" "$current_dir/$exp_ant3_closest_seq_t$tsv" "$current_dir/$exp_ant3_iden_t$tsv" "$current_dir/$exp_stress_gene$tsv" "$current_dir/$exp_stress_seq$tsv" "$current_dir/$exp_virulence_gene$tsv" "$current_dir/$exp_virulence_seq$tsv" > "$current_dir/Final_noheader.tsv"
#Concat header with all results
cat "$current_dir/header.tsv" "$current_dir/Final_noheader.tsv" > "$current_dir/$final_name" 

#cut other antibiotics column when they choose to input NA
if [[ "$ant_2_lower" == "na" && "$ant_3_lower" == "na" ]]; then
 new_name="AMRFinderPlus_${seq_input_lower}_assembly_result_new.tsv"
 cut -d$'\t' -f 1-7 "$current_dir/$final_name" > "$current_dir/$new_name"
 rm "$current_dir/$final_name"
elif [[ "$ant_3_lower" == "na" ]]; then
 new_name="AMRFinderPlus_${seq_input_lower}_assembly_result_new.tsv"
 cut -d$'\t' -f 1-13 "$current_dir/$final_name" > "$current_dir/$new_name"
 rm "$current_dir/$final_name"
fi
#Remove all intermediate files
rm "$current_dir/genome_id.tsv" "$current_dir/$ant_1$tsv" "$current_dir/$exp_ant1_gene_t$tsv" "$current_dir/$exp_ant1_seq_t$tsv" "$current_dir/$exp_ant1_subclass_t$tsv" "$current_dir/$exp_ant1_closest_seq_t$tsv" "$current_dir/$exp_ant1_iden_t$tsv" "$current_dir/$ant_2$tsv" "$current_dir/$exp_ant2_gene_t$tsv" "$current_dir/$exp_ant2_seq_t$tsv" "$current_dir/$exp_ant2_subclass_t$tsv" "$current_dir/$exp_ant2_closest_seq_t$tsv" "$current_dir/$exp_ant2_iden_t$tsv" "$current_dir/$ant_3$tsv" "$current_dir/$exp_ant3_gene_t$tsv" "$current_dir/$exp_ant3_seq_t$tsv" "$current_dir/$exp_ant3_subclass_t$tsv" "$current_dir/$exp_ant3_closest_seq_t$tsv" "$current_dir/$exp_ant3_iden_t$tsv" "$current_dir/header.tsv" "$current_dir/Final_noheader.tsv" "$current_dir/$exp_stress_gene$tsv" "$current_dir/$exp_stress_seq$tsv" "$current_dir/$exp_virulence_gene$tsv" "$current_dir/$exp_virulence_seq$tsv"
