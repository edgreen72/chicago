#!/usr/local/bin/tcsh

### Synatx works on edser (not edser2)

### Do a git clone of
### https://github.com/edgreen72/chicago
### into local filesystem and do make
setenv CHICAGO /soe/ed/src/chicago

### The folder where the raw, input fastq reads are
setenv DD /projects/redser2/raw/
setenv GENOME /projects/redser2/genomes/hg19/hg19.fa # The bwa indexed genome to align against
setenv FJT ${CHICAGO}/fq-jt
setenv IOS ${CHICAGO}/innies-outies-same-strand-lenhist.pl
setenv SST ${CHICAGO}/sesam2table.pl

### Junction sequence for truncating reads
setenv JUNCTION GATCGATC

### Make output directories
mkdir dat aln seq

### REPLACE FILE_NAME.R1, FILE_NAME.R2, and LIB-NAME
### with appropriate values

${FJT} -f ${DD}/FILE_NAME.R1.fastq.gz -t ${JUNCTION} -l 4 > seq/LIB-NAME.R1.t.fq
${FJT} -f ${DD}/FILE_NAME.R2.fastq.gz -t ${JUNCTION} -l 4 > seq/LIB-NAME.R2.t.fq
gzip seq/LIB-NAME.*.t.fq

foreach LIB ( "LIB-NAME" )
   echo ${LIB}
   bwa aln -t 12 ${GENOME} seq/${LIB}.R1.t.fq.gz > aln/${LIB}.R1.sai
   bwa aln -t 12 ${GENOME} seq/${LIB}.R2.t.fq.gz > aln/${LIB}.R2.sai

   bwa samse -f aln/${LIB}.R1.sam ${GENOME} aln/${LIB}.R1.sai seq/${LIB}.R1.t.fq.gz
   bwa samse -f aln/${LIB}.R2.sam ${GENOME} aln/${LIB}.R2.sai seq/${LIB}.R2.t.fq.gz

   samtools view -Sb aln/${LIB}.R1.sam | samtools sort - aln/${LIB}.R1.s
   samtools view -Sb aln/${LIB}.R2.sam | samtools sort - aln/${LIB}.R2.s

   preseq lc_extrap -B aln/${LIB}.R1.s.bam -e 105000000 -s 5000000 > dat/${LIB}.R1.preseq.dat
   preseq lc_extrap -B aln/${LIB}.R2.s.bam -e 105000000 -s 5000000 > dat/${LIB}.R2.preseq.dat

   ${STT} -f aln/${LIB}.R1.sam -r aln/${LIB}.R2.sam \
	  -s dat/${LIB}.same.tab -d dat/${LIB}.diff.tab
   ${IOS} -s dat/${LIB}.same.tab > dat/${LIB}.ios.txt

   rm aln/${LIB}.R1.sai
   rm aln/${LIB}.R1.sai   
   end

   
