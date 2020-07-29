#!/usr/bin/env nextflow

params.gsurl37 = "gs://gatk-legacy-bundles/b37"
params.exome_file = "/data/bdigby/WES/assets/20130108.exome.targets.bed"

process dl_fasta{

	publishDir "$baseDir/DNAref", mode: 'copy'

	output:
	tuple file('*noChr.fasta'), file('*noChr.fasta.fai') into (fasta_bwa, fasta_seqza, fasta_msi, fasta_dict, fasta_2bit, fasta_exome_biall, fasta_wgs_biall)

	script:
	"""
        ##http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
        gsutil cp ${params.gsurl37}/human_g1k_v37.fasta.gz ./human_g1k_v37.fasta.gz
        gunzip -c human_g1k_v37.fasta.gz | sed 's/>chr/>/g' > human_g1k_v37.noChr.fasta
        samtools faidx human_g1k_v37.noChr.fasta
        """
}


process dict_pr {

	publishDir "$baseDir/DNAref", mode: 'copy'
	
	input:
	tuple file(fa), file(fai) from fasta_dict

	output:
	file('*.dict') into dict_win
 	tuple file(fa), file(fai), file('*.dict') into (fasta_dict_exome, fasta_dict_wgs, fasta_dict_gensiz, fasta_dict_gridss)

	script:
	"""
	DICTO=\$(echo $fa | sed 's/fasta/dict/')
        picard CreateSequenceDictionary \
        R=$fa \
        O=\$DICTO
	"""
}


process dbsnp_dl{

	output:
	file('*.vcf') into vcf_tabix
	file('KG_phase1.snps.high_confidence.*.vcf') into ascatloci

	script:	
	"""
	gsutil cp ${params.gsurl37}/1000G_phase1.snps.high_confidence.b37.vcf.gz ./KG_phase1.snps.high_confidence.b37.vcf.gz
	gsutil cp ${params.gsurl37}/dbsnp_138.b37.vcf.gz ./dbsnp_138.b37.vcf.gz
    	gsutil cp ${params.gsurl37}/hapmap_3.3.b37.vcf.gz ./hapmap_3.3.b37.vcf.gz
    	gsutil cp ${params.gsurl37}/1000G_omni2.5.b37.vcf.gz ./KG_omni2.5.b37.vcf.gz
    	gsutil cp ${params.gsurl37}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz ./Mills_KG_gold.indels.b37.vcf.gz
    	gunzip -cd dbsnp_138.b37.vcf.gz | sed 's/chr//g' > dbsnp_138.b37.vcf
    	gunzip -cd hapmap_3.3.b37.vcf.gz | sed 's/chr//g' > hapmap_3.3.b37.sites.vcf
    	gunzip -cd KG_omni2.5.b37.vcf.gz | sed 's/chr//g' > KG_omni2.5.b37.vcf
    	gunzip -cd KG_phase1.snps.high_confidence.b37.vcf.gz | sed 's/chr//g' > KG_phase1.snps.high_confidence.b37.vcf
    	gunzip -cd Mills_KG_gold.indels.b37.vcf.gz | sed 's/chr//g' > Mills_KG_gold.indels.b37.vcf
	"""
}

process ascat_loci{
	
	publishDir "$baseDir/DNAref", mode: 'copy'

	input:
	file(vcf) from ascatloci

	output:
	file('*loci') into complete_ascat

	script:
	"""
  	LOCIFILE=\$(echo $vcf | sed 's/vcf/maf0.3.loci/')
  	cat $vcf | \
  	perl -ane '@s=split(/[=\\;]/,\$F[7]);if(\$s[3]>0.3){print "\$F[0]\\t\$F[1]\\n";}' > \$LOCIFILE
	"""
}


process dict_pr2 {
	
	publishDir "$baseDir/DNAref", mode: 'copy'

	input:
	file(win_dict) from dict_win

	output:
	file('*') into complete_dict

	script:
	"""
  	perl -ane 'if(\$F[0]=~m/SQ\$/){@sc=split(/:/,\$F[1]);@ss=split(/:/,\$F[2]); if(\$sc[1]!~m/[GLMT]/){ print "\$sc[1]\\t\$ss[1]\\n";}}' $win_dict > seq.dict.chr-size
	bedtools makewindows -g seq.dict.chr-size -w 35000000 | perl -ane 'if(\$F[1]==0){\$F[1]++;};print "\$F[0]:\$F[1]-\$F[2]\n";' > 35MB-window.bed
	"""
}


process exome_file {
	
	publishDir "$baseDir/DNAref", mode: 'copy'

	input:
	path(exome) from params.exome_file

	output:
	file('*.file.bed') into exome_bed

	script:
	"""
	if [[ $exome =~ bed\$ ]]; then
        ##remove any non-chr, coord lines in top of file
        CHR=\$(tail -n1 $exome | perl -ane 'print \$F[0];')
        if [[ \$CHR =~ "chr" ]]; then
           perl -ane 'if(\$F[0]=~m/^chr/){print \$_;}' $exome >  exome.file.bed
        else
           perl -ane 'if(\$F[0]=~m/^[0-9MXY]/){print \$_;}' $exome >  exome.file.bed
        fi
       
	else
        echo "BED file is not a BED file, please retry"
        exit 147
    	fi
	"""
}

process exome_bed_pr {

 	publishDir "$baseDir/DNAref", mode: "copy", pattern: "*[.interval_list,.bed]"

  	input:
  	tuple file(fa), file(fai), file(dict) from fasta_dict_exome
  	file(exomelift) from exome_bed

  	output:
  	file("*.bed.interval_list") into complete_exome
  	file("*.bed") into (exome_tabix, exome_biallgz)

  	script:
  	"""
  	##must test if all chr in fasta are in exome, else manta cries
  	##must test if all regions are greater than length zero or strelka cries
  	##must test if all seq.dict chrs are in bed and only they or BedToIntervalList cries
  	perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' $exomelift | grep -v chrM | sed 's/chr//g' > tmp.bed
   	grep @SQ $dict | cut -f2 | sed 's/SN://' | while read CHR; do
   	TESTCHR=\$(awk -v chrs=\$CHR '\$1 == chrs' tmp.bed | wc -l)
   	if [[ \$TESTCHR != 0 ]];then
    	 awk -v chrs=\$CHR '\$1 == chrs' tmp.bed
   	fi
  	done >> tmp.dict.bed
  
	##always make interval list so we are in line with fasta
  	picard BedToIntervalList I=tmp.dict.bed O=exome.interval_list SD=$dict
  	##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  	perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' exome.interval_list > exome.bed.interval_list
  	##output BED
  	grep -v "@" exome.bed.interval_list | cut -f 1,2,3,5 > exome.bed
  	"""
}


process tabix_files{
	
	publishDir "$baseDir/DNAref/exome/", mode: 'copy', pattern: 'exome'

	input:
	file(bed) from exome_tabix

	output:
	tuple file("${bed}.gz"), file("${bed}.gz.tbi") into complete_tabix

	script:
	"""
	bgzip $bed
	tabix $bed".gz"
	"""
}

process exome_biall{
	
	publishDir "$baseDir/DNAref/exome", mode: 'copy'

	input:
	file(exomebed) from exome_biallgz
	tuple file(fasta), file(fai) from fasta_exome_biall

	output:
  	tuple file('af-only-gnomad.*.noChr.vcf.gz'), file('af-only-gnomad.*.noChr.vcf.gz.tbi') into exome_biallelicgz
  	file('exome.biall.bed') into pcgrtoml_exome

  	script:
  	"""
  	cut -f 1,2,3 $exomebed > exome.biall.bed
  	
    	gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ./
   	bgzip af-only-gnomad.raw.sites.vcf
    	tabix af-only-gnomad.raw.sites.vcf.gz
    	gunzip -c af-only-gnomad.raw.sites.vcf.gz |
    	bcftools view -R exome.biall.bed af-only-gnomad.raw.sites.vcf.gz | bcftools sort -T '.' > af-only-gnomad.exomerh.hg19.noChr.vcf
    	perl ${workflow.projectDir}/bin/reheader_vcf_fai.pl af-only-gnomad.exomerh.hg19.noChr.vcf $fai > af-only-gnomad.exome.hg19.noChr.vcf
    	bgzip af-only-gnomad.exome.hg19.noChr.vcf
    	tabix af-only-gnomad.exome.hg19.noChr.vcf.gz
	"""
}


process index_feature_files {
	
	publishDir "$baseDir/DNAseq", mode: 'copy'

	input:
  	file(tbtbx) from vcf_tabix.flatten()

  	output:
  	file('*') into indexfeatured

  	script:
  	"""
  	bgzip $tbtbx
  	gatk IndexFeatureFile -F $tbtbx".gz"
  	"""
}

process bwa_index {

        publishDir "$baseDir/DNAref", mode: 'copy'

        input:
        tuple file(fa), file(fai) from fasta_bwa

        output:
        file('*') into complete_bwa

        script:
        """
        bwa index -a bwtsw $fa
        """
}
