#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '-------------------------------------------------------------------'
  log.info 'NEXTFLOW: MAKE HUMAN GRCH38 REFERENCE FOR DNASEQ NEXTFLOW PIPELINES'
  log.info '-------------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run main.nf'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --version   STRING    GRCh37 or GRCh38 (default)'
  log.info '    --outDir    STRING    output directory path; NB ${params.version} dir is created therein'
  log.info '    --exomebedurl     STRING      URL to exome bed file for intervals; NB assumes GRCh37 as all Illumina exomes are'
  log.info ''
  exit 1
}

/* 0.0: Global Variables
*/
if(params.outDir){
  params.refDir = "${params.outDir}"
}
if(!params.outDir){
  params.refDir = "${workflow.launchDir}/${params.version}"
}

//base URL for GRCh37, 38
params.gsurl37 = "gs://gatk-legacy-bundles/b37"
params.gsurl38 = "gs://genomics-public-data/resources/broad/hg38/v0"

//lowercase version
params.versionlc = "${params.version}".toLowerCase()

/* 1.0: Download GATK4 resource bundle fasta
*/
process fasta_dl {

  publishDir path: "$params.refDir", mode: "copy"
  validExitStatus 0,1,2
  errorStrategy 'retry'
  maxRetries 3

  output:
  set file('*noChr.fasta'), file('*noChr.fasta.fai') into (fasta_bwa, fasta_seqza, fasta_msi, fasta_dict, fasta_2bit)

  script:
  if( params.version == 'GRCh37' )
    """
    ##http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    gsutil cp ${params.gsurl37}/human_g1k_v37.fasta.gz ./human_g1k_v37.fasta.gz
    gunzip -c human_g1k_v37.fasta.gz | sed 's/>chr/>/g' > human_g1k_v37.noChr.fasta
    samtools faidx human_g1k_v37.noChr.fasta
    """

  else
    """
    ##http://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    ##moved to Verily as gs bucket more reliable
    gsutil cp gs://genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa ./
    cat GRCh38_Verily_v1.genome.fa | sed 's/>chr/>/g' > GRCh38_Verily_v1.genome.noChr.fasta
    samtools faidx GRCh38_Verily_v1.genome.noChr.fasta
    """
}

/* 1.1: Dictionary for fasta
*/
process dict_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_dict

  output:
  file('*.dict') into dict_win
  set file(fa), file(fai), file('*.dict') into (fasta_dict_exome, fasta_dict_wgs, fasta_dict_gensiz)

  """
  DICTO=\$(echo $fa | sed 's/fasta/dict/')
  picard CreateSequenceDictionary \
    R=$fa \
    O=\$DICTO
  """
}

/* 1.2: Download GATK4 resource bundle dbsnp, Mills
*/
process dbsnp_dl {

  validExitStatus 0,1,2
  errorStrategy 'retry'
  maxRetries 3

  output:
  file('*.vcf') into vcf_tabix
  file('KG_phase1.snps.high_confidence.*.vcf') into ascatloci

  when:
  !params.nodbsnp

  script:
  if( params.version == 'GRCh37' )
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
  else
    """
    gsutil cp gs://genomics-public-data/cwl-examples/gdc-dnaseq-cwl/input/dbsnp_144.hg38.vcf.gz ./dbsnp_144.hg38.vcf.gz
    gsutil cp ${params.gsurl38}/1000G_phase1.snps.high_confidence.hg38.vcf.gz ./KG_phase1.snps.high_confidence.hg38.vcf.gz
    gsutil cp ${params.gsurl38}/Homo_sapiens_assembly38.dbsnp138.vcf ./Homo_sapiens_assembly38.dbsnp138.vcf
    gsutil cp ${params.gsurl38}/hapmap_3.3.hg38.vcf.gz ./hapmap_3.3.hg38.vcf.gz
    gsutil cp ${params.gsurl38}/1000G_omni2.5.hg38.vcf.gz ./KG_omni2.5.hg38.vcf.gz
    gsutil cp ${params.gsurl38}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ./Mills_KG_gold.indels.hg38.vcf.gz

    gunzip -cd dbsnp_144.hg38.vcf.gz | sed 's/chr//g' > dbsnp_144.hg38.vcf
    gunzip -cd hapmap_3.3.hg38.vcf.gz | sed 's/chr//g' > hapmap_3.3.hg38.vcf
    gunzip -cd KG_omni2.5.hg38.vcf.gz | sed 's/chr//g' > KG_omni2.5.hg38.vcf
    gunzip -cd KG_phase1.snps.high_confidence.hg38.vcf.gz | sed 's/chr//g' > KG_phase1.snps.high_confidence.hg38.vcf
    gunzip -cd Mills_KG_gold.indels.hg38.vcf.gz | sed 's/chr//g' > Mills_KG_gold.indels.hg38.vcf
    """
}

/* 1.3: KG ASCAT loci
*/
process ascat_loci {

  publishDir path: "$params.refDir", mode: "copy"

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

/* 2.0: Fasta processing
*/
process fasta_pr {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_bwa

  output:
  file('*') into complete_bwa

  script:
  """
  ##https://gatkforums.broadinstitute.org/gatk/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk
  bwa index -a bwtsw $fa
  """
}

/* 2.1: Dict processing
*/
process dict_pr2 {

  publishDir path: "$params.refDir", mode: "copy"

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

/* 3.0: Exome bed file and liftOver
* of course why would Illumina supply exomes in GRCh38 in 2019?!?!?!?!
*/
process lift_over {

  publishDir path: "$params.refDir", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3

  output:
  set file("exome.bed"), file("README.exome.bed") into exome_bed

  when:
  params.exomebedurl

  script:
  """
  ##download URL
  echo "Exome bed used here is from:" > README.exome.bed
  echo ${params.exomebedurl} >> README.exome.bed

  wget ${params.exomebedurl}
  if [[ ${params.exomebedurl} =~ zip\$ ]]; then
    unzip -p *.zip > exome.url.bed
  elif [[ ${params.exomebedurl} =~ bed\$ ]]; then
    mv *bed exome.url.bed
  else
    echo "No ZIP or BED files resulting from ${params.exomebedurl}"
    echo "Please try another URL with ZIP or BED file resulting"
    exit 0
  fi

  if [[ ${params.version} != "GRCh37" ]]; then
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    liftOver exome.url.bed hg19ToHg38.over.chain.gz exome.bed unmapped
    echo "${params.exomebedurl} was liftOver'd from hg19 to hg38" >> README.exome.bed
  else
    mv exome.url.bed exome.bed
  fi
  """
}

/* 3.11: Parse bed for exome
*/
process exome_bed {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai), file(dict) from fasta_dict_exome
  set file(exomebed), file(readme) from exome_bed

  output:
  file('*') into complete_exome
  set file(fa), file(fai), file('exome.bed'), file(dict) into (exome_tabix, exome_fasta_biallgz)

  when:
  params.exomebedurl

  script:
  """
  ##must test if all chr in fasta are in exome, else manta cries
  ##must test if all regions are greater than length zero or strelka cries
  ##must test if all seq.dict chrs are in bed and only they or BedToIntervalList cries
  perl -ane 'if(\$F[1] == \$F[2]){\$F[2]++;} if(\$F[0] !~m/^chrM/){print join("\\t", @F[0..\$#F]) . "\\n";}' $exomebed | grep -v chrM | sed 's/chr//g' > tmp.bed

   grep @SQ $dict | cut -f2 | sed 's/SN://' | while read CHR; do
   TESTCHR=\$(awk -v chrs=\$CHR '\$1 == chrs' tmp.bed | wc -l)
   if [[ \$TESTCHR != 0 ]];then
    awk -v chrs=\$CHR '\$1 == chrs' tmp.bed
   fi
  done >> tmp.dict.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.dict.bed O=exome.bed.interval_list SD=$dict

  ##BedToIntervalList (reason unknown) makes 1bp interval to 0bp interval, replace with original
  perl -ane 'if(\$F[0]=~m/^@/){print \$_;next;} if(\$F[1] == \$F[2]){\$f=\$F[1]; \$f--; \$F[1]=\$f; print join("\\t", @F[0..\$#F]) . "\\n";} else{print \$_;}' exome.bed.interval_list > exome.bed.interval_list1
  mv exome.bed.interval_list1 exome.bed.interval_list

  ##output BED
  grep -v "@" exome.bed.interval_list | cut -f 1,2,3,5 > exome.bed
  """
}

/* 3.12: create bed for WGS
*/
process wgs_bed {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai), file(dict) from fasta_dict_wgs

  output:
  file('*') into complete_wgs
  set file(fa), file(fai), file('wgs.bed'), file(dict) into (wgs_tabix, wgs_fasta_biallgz)

  script:
  """
  ##WGS intervals = 1-LN for each chr
  grep @SQ $dict | cut -f 2,3 | perl -ane '\$chr=\$F[0];\$chr=~s/SN://;\$end=\$F[1];\$end=~s/LN://;print "\$chr\\t0\\t\$end\\n";' > tmp.wgs.dict.bed

  ##always make interval list so we are in line with fasta
  picard BedToIntervalList I=tmp.wgs.dict.bed O=wgs.bed.interval_list SD=$dict

  ##output BED
  grep -v "@" wgs.bed.interval_list | cut -f 1,2,3,5 > wgs.bed
  """
}

/* 3.2: Tabix those requiring tabixing
*/
wgs_tabix.concat(exome_tabix).set { bint_tabix }
process tabix_files {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai), file(bed), file(dict) from bint_tabix

  output:
  file('*') into complete_tabix

  script:
  """
  ##tabix
  bgzip $bed
  tabix $bed".gz"
  """
}

/* 3.31: Create Mutect2 af-only-gnomad file
*/
process exome_biall {

  publishDir path: "$params.refDir", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3
  label 'half_cpu_mem'

  input:
  set file(fa), file(fai), file(exomebed), file(dict) from exome_fasta_biallgz

  output:
  set file('af-only-gnomad.exome.*.noChr.vcf.gz'), file('af-only-gnomad.exome.*.noChr.vcf.gz.tbi') into exome_biallelicgz

  script:
  """
  cut -f 1,2,3 $exomebed > exome.biall.bed

  if [[ ${params.version} == "GRCh37" ]];then

    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ./
    bgzip af-only-gnomad.raw.sites.vcf
    tabix af-only-gnomad.raw.sites.vcf.gz
    bcftools view -R exome.biall.bed af-only-gnomad.raw.sites.vcf.gz | bcftools sort -T '.' | bgzip > af-only-gnomad.exome.hg19.noChr.vcf.gz
    tabix af-only-gnomad.exome.hg19.noChr.vcf.gz

  else

    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ./
    gunzip -c af-only-gnomad.hg38.vcf.gz | sed 's/chr//' | bgzip > af-only-gnomad.hg38.noChr.vcf.gz
    tabix af-only-gnomad.hg38.noChr.vcf.gz

    bcftools view -R exome.biall.bed af-only-gnomad.hg38.noChr.vcf.gz | bcftools sort -T '.' | bgzip > af-only-gnomad.exome.hg38.noChr.vcf.gz
    tabix af-only-gnomad.exome.hg38.noChr.vcf.gz
  fi
  """
}

/* 3.32: Create Mutect2 af-only-gnomad file
*/
process wgs_biall {

  publishDir path: "$params.refDir", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3
  label 'half_cpu_mem'

  input:
  set file(fa), file(fai), file(wgsbed), file(dict) from wgs_fasta_biallgz

  output:
  set file('af-only-gnomad.wgs.*.noChr.vcf.gz'), file('af-only-gnomad.wgs.*.noChr.vcf.gz.tbi') into wgs_biallelicgz

  script:
  """
  cut -f 1,2,3 $wgsbed > wgs.biall.bed

  if [[ ${params.version} == "GRCh37" ]];then

    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ./
    bgzip af-only-gnomad.raw.sites.vcf
    tabix af-only-gnomad.raw.sites.vcf.gz
    bcftools view -R wgs.biall.bed af-only-gnomad.raw.sites.vcf.gz | bcftools sort -T '.' | bgzip > af-only-gnomad.wgs.hg19.noChr.vcf.gz
    tabix af-only-gnomad.wgs.hg19.noChr.vcf.gz

  else

    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz ./
    gunzip -c af-only-gnomad.hg38.vcf.gz | sed 's/chr//' | bgzip > af-only-gnomad.hg38.noChr.vcf.gz
    tabix af-only-gnomad.hg38.noChr.vcf.gz

    bcftools view -R wgs.biall.bed af-only-gnomad.hg38.noChr.vcf.gz | bcftools sort -T '.' | bgzip > af-only-gnomad.wgs.hg38.noChr.vcf.gz
    tabix af-only-gnomad.wgs.hg38.noChr.vcf.gz
  fi
  """
}

/* 4.0 Index various VCFs
*/
process indexfeature_files {

  publishDir path: "$params.refDir", mode: "copy"

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

/* 5.0: Sequenza GC bias
*/
process seqnza {

  publishDir path: "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_seqza

  output:
  file('*') into sequenzaout

  script:
  """
  GENOMEGC50GZ=\$(echo $fa | sed -r 's/.fasta/.gc50Base.txt.gz/')
  sequenza−utils.py GC-windows −w 50 $fa | gzip > \$GENOMEGC50GZ
  """
}

/* 6.0: MSIsensor microsatellites
*/
process msisen {

  publishDir "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai) from fasta_msi

  output:
  file('*') into completedmsisensor

  script:
  """
  msisensor scan -d $fa -o msisensor_microsatellites.list
  """
}

/* 7.0: PCGR/CPSR data bundle
*/
process pcgr_data {

  publishDir "$params.refDir/pcgr", mode: "copy"
  errorStrategy 'retry'
  maxRetries 3

  output:
  file('data/') into completedpcgrdb
  file("data/${params.versionlc}/.vep/") into pcgrdbvep
  file("data/${params.versionlc}/RELEASE_NOTES") into pcgrreleasenotes

  when:
  !params.nopcgr

  script:
  if( params.version == "GRCh37" )
    """
    wget ${params.pcgrURL37}
    tar -xf *.tgz
    rm -rf *.tgz
    """
  else
    """
    wget ${params.pcgrURL38}
    tar -xf *.tgz
    rm -rf *.tgz
    """
}

/* 3.6: PCGR/CPSR VEP cache
*/
process vepdb {

  publishDir "$params.refDir/pcgr/data/${params.versionlc}", mode: "copy"

  input:
  file(releasenotes) from pcgrreleasenotes
  file(pcgrdbvepdir) from pcgrdbvep

  output:
  file('.vep/homo_sapiens') into complete_vepdb

  """
  #! /bin/bash
  ##build VEP cache using PCGR Singularity container 'vep_install' script
  ##however PCGR installs a version of vep cache, so test that matches required version, and only install if not

  ##variables for install and test
  VEP_INSTALL=\$(find /opt/miniconda/envs/somatic_exome_n-of-1/share/*/vep_install)
  VEP_VERSION=\$(cat $releasenotes | perl -ane 'if(\$F[0] eq "VEP"){@s=split(/\\./,\$F[5]); \$v=\$s[0]; \$v=~s/v//; print \$v;}')

  ls $pcgrdbvepdir/homo_sapiens/ | cut -d "_" -f 1 > test.match
  if [[ \$(grep \$VEP_VERSION test.match | wc -l) != 1 ]];then
    \$VEP_INSTALL \
      --AUTO cf \
      --CACHE_VERSION \$VEP_VERSION \
      --CACHEDIR "./" \
      --SPECIES "homo_sapiens" \
      --ASSEMBLY ${params.version} \
      --NO_UPDATE \
      --NO_HTSLIB \
      --NO_BIOPERL \
      --NO_TEST
  fi
  """
}

/* 3.7: GenomeSize.xml for Pisces
*/
process gensizxml {

  publishDir "$params.refDir", mode: "copy"

  input:
  set file(fa), file(fai), file(dict) from fasta_dict_gensiz

  output:
  file('*') into complete_gensiz

  script:
  """
  echo "<sequenceSizes genomeName=\"$dict\">" > GenomeSize.xml
  grep "@SQ" $dict | while read LINE; do
    CONTIGNAME=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[1]);print \$s[1];' | sed 's/chr//')
    TOTALBASES=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[2]);print \$s[1];')
    MD5SUM=\$(echo \$LINE | perl -ane '@s=split(/:/,\$F[3]);print \$s[1];')
    echo -e "\\t<chromosome fileName=\"$fa\" contigName=\"\$CONTIGNAME\" totalBases=\"\$TOTALBASES\" isCircular=\"false\" md5=\"\$MD5SUM\" ploidy=\"2\" knownBases=\"\$TOTALBASES\" type=\"Chromosome\" />" >> GenomeSize.xml
  done
  echo "</sequenceSizes>" >> GenomeSize.xml
  """
}
