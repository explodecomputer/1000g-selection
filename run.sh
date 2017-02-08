#!/bin/bash


for i in {1..22}
do
	wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
done

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

grep -w "EUR" integrated_call_samples_v3.20130502.ALL.panel | awk '{ print $1, $1 }' > eur.sample

parallel 'plink1.90 \
	--vcf ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
	--make-bed \
	--out chr{}' ::: {1..22}

plink1.90 \
	--vcf ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz \
	--make-bed \
	--out chr23

doawk() {
	cp chr$1.bim chr$1.bim.orig;
	awk '{
		if (($5 == "A" || $5 == "T" || $5 == "C" || $5=="G") &&  ($6 == "A" || $6 == "T" || $6 == "C" || $6=="G")) 
			print $1, "chr"$1":"$4":SNP", $3, $4, $5, $6;
		else 
			print $1, "chr"$1":"$4":INDEL", $3, $4, $5, $6;
		}' chr$1.bim.orig > chr$1.bim
}
export -f doawk
parallel --progress doawk {} ::: {1..23}

parallel --progress 'plink1.90 \
	--bfile chr{} \
	--keep eur.sample \
	--freq \
	--out chr{} > /dev/null 2>&1' ::: {1..23}

for i in {1..22}
do
	sed 1d chr${i}.frq | awk '{ print $2, $3, $4, $5 }'
done > 1kg_phase3_eur_aut.frq
gzip 1kg_phase3_eur_aut.frq
zcat 1kg_phase3_eur_aut.frq.gz | awk '{ if($4 > 0) print $0 }' | gzip -c > 1kg_phase3_eur_aut_polymorphic.frq.gz
zcat 1kg_phase3_eur.frq.gz | awk '{ if($4 > 0) print $0 }' | gzip -c > 1kg_phase3_eur_polymorphic.frq.gz


