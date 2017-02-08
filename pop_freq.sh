pop=( AFR AMR EAS EUR SAS )

for i in "${pop[@]}"
do
	grep -w ${i} integrated_call_samples_v3.20130502.ALL.panel | cut -f 1 | awk '{ print $1, $1 }' > pop_freq/$i.sample
done


# Get frequencies across all populations
parallel --progress 'plink1.90 \
    --bfile chr{} \
    --freq \
    --out pop_freq/ALL_chr{} > /dev/null 2>&1' ::: {1..23}

rm -f pop_freq/snplist_0.01.txt 
for i in {1..23}
do
echo ${i}
awk '{ if($5 > 0.01) print $2 }' pop_freq/ALL_chr${i}.frq >> pop_freq/snplist_0.01.txt
done 


parallel --progress 'plink1.90 \
    --bfile chr{} \
    --keep pop_freq/AFR.sample \
    --extract pop_freq/snplist_0.01.txt \
    --freq \
    --out pop_freq/AFR_chr{} > /dev/null 2>&1' ::: {1..23}

parallel --progress 'plink1.90 \
    --bfile chr{} \
    --keep pop_freq/AMR.sample \
    --extract pop_freq/snplist_0.01.txt \
    --freq \
    --out pop_freq/AMR_chr{} > /dev/null 2>&1' ::: {1..23}

parallel --progress 'plink1.90 \
    --bfile chr{} \
    --keep pop_freq/EAS.sample \
    --extract pop_freq/snplist_0.01.txt \
    --freq \
    --out pop_freq/EAS_chr{} > /dev/null 2>&1' ::: {1..23}

parallel --progress 'plink1.90 \
    --bfile chr{} \
    --keep pop_freq/EUR.sample \
    --extract pop_freq/snplist_0.01.txt \
    --freq \
    --out pop_freq/EUR_chr{} > /dev/null 2>&1' ::: {1..23}

parallel --progress 'plink1.90 \
    --bfile chr{} \
    --keep pop_freq/SAS.sample \
    --extract pop_freq/snplist_0.01.txt \
    --freq \
    --out pop_freq/SAS_chr{} > /dev/null 2>&1' ::: {1..23}

for i in "${pop[@]}"
do
echo ${i}
cat pop_freq/${i}_chr1.frq > pop_freq/${i}.frq
for j in {2..23}
do
echo ${j}
sed 1d pop_freq/${i}_chr${j}.frq >> pop_freq/${i}.frq
done
gzip -f pop_freq/${i}.frq
done




parallel --progress 'zgrep -v "^#" ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | cut -f 1,2,3,8 | grep "AA=" | cut -d ";" -f 1,11 | cut -d "|" -f 1 > aa{}.txt' ::: {1..22}

parallel --progress 'tr ";" "\t" < aa{}.txt | sed "s@AA=@@g" | cut -f 1,2,3,5 | gzip -c > aa{}.txt.gz' ::: {1..22}

for i in {1..22}
do
	cat aa${i}.txt.gz
done > aa.txt.gz


zcat aa.txt.gz | awk '{ if($4 == "A" || $4 == "G" || $4 == "C" || $4 == "T") print "chr"$1":"$2":SNP", $4 }' | gzip -c > aa_clean.txt.gz


