# Analysis of Lu et al. 2019
## Andrew D. Johnston
### 03/23/20

For the revision of "A cellular stress response induced by the CRISPR/dCas9 activation system is not heritable through cell divisions", Lu et al 2019 data was analyzed to examine if a similar stress response was observed in their experiments. Specifically, they transfected Cas9 to knockout STAT3 in SKOV3 cells (ovarian cancer cell line). 

Not only did they compare the wild type SKOV3 cell line with KO-STAT3 via CRISPR cells BUT also examined transfection with Cas9 w/o gRNA (unclear if scrambled or absent) vs. the CRIPSR KO cell line. 

#### Table of Contents
1. [Downloading data](#download)
2. [Trimming and Aligning](#trim)



<a name="download"></a>
1. Downloading data

I need to download the raw fastq data from:
	
Lu T, Bankhead A 3rd, Ljungman M, Neamati N. Multi-omics profiling reveals key signaling pathways in ovarian cancer controlled by STAT3. Theranostics 2019;9(19):5478-5496. PMID: 31534498

GEO URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134375


"paired-end 50bp (SKOV3) sequencing" - so I should split files

SRA codes:
SRR9694244
.
.
SRR9694249


SRR9694262
.
.
SRR9694267

```bash
 # set directory
pwd
 #/home/greally-lab/Claudia_Andrew/Lu_et_al/fastq

 #download the files 
for i in {5..9}; 
do 
echo SRR969424$i
qsub -S /bin/bash -N down_$i -cwd -l h_vmem=10G -j y << EOF
module load sra-toolkit/2.9.6-1 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip SRR969424$i
EOF
done

for i in {2..7}; 
do 
echo SRR969426$i
qsub -S /bin/bash -N down_$i -cwd -l h_vmem=10G -j y << EOF
module load sra-toolkit/2.9.6-1 
fastq-dump --gzip --skip-technical --readids --dumpbase --split-files --clip SRR969426$i
EOF
done

 # Renaming the files
mv SRR9694244_1.fastq.gz SKOV3_WT_rep1_r1.fq.gz 
mv SRR9694244_2.fastq.gz SKOV3_WT_rep1_r2.fq.gz 
mv SRR9694245_1.fastq.gz SKOV3_WT_rep2_r1.fq.gz 
mv SRR9694245_2.fastq.gz SKOV3_WT_rep2_r2.fq.gz 
mv SRR9694246_1.fastq.gz SKOV3_WT_rep3_r1.fq.gz 
mv SRR9694246_2.fastq.gz SKOV3_WT_rep3_r2.fq.gz

mv SRR9694247_1.fastq.gz SKOV3_STAT3_rep1_r1.fq.gz
mv SRR9694247_2.fastq.gz SKOV3_STAT3_rep1_r2.fq.gz 
mv SRR9694248_1.fastq.gz SKOV3_STAT3_rep2_r1.fq.gz
mv SRR9694248_2.fastq.gz SKOV3_STAT3_rep2_r2.fq.gz
mv SRR9694249_1.fastq.gz SKOV3_STAT3_rep3_r1.fq.gz
mv SRR9694249_2.fastq.gz SKOV3_STAT3_rep3_r2.fq.gz

mv SRR9694262_1.fastq.gz SKOV3_Cas9_rep1_r1.fq.gz
mv SRR9694262_2.fastq.gz SKOV3_Cas9_rep1_r2.fq.gz
mv SRR9694263_1.fastq.gz SKOV3_Cas9_rep2_r1.fq.gz
mv SRR9694263_2.fastq.gz SKOV3_Cas9_rep2_r2.fq.gz
mv SRR9694264_1.fastq.gz SKOV3_Cas9_rep3_r1.fq.gz
mv SRR9694264_2.fastq.gz SKOV3_Cas9_rep3_r2.fq.gz

mv SRR9694265_1.fastq.gz SKOV3_STAT3KO_rep1_r1.fq.gz
mv SRR9694265_2.fastq.gz SKOV3_STAT3KO_rep1_r2.fq.gz
mv SRR9694266_1.fastq.gz SKOV3_STAT3KO_rep2_r1.fq.gz
mv SRR9694266_2.fastq.gz SKOV3_STAT3KO_rep2_r2.fq.gz
mv SRR9694267_1.fastq.gz SKOV3_STAT3KO_rep3_r1.fq.gz
mv SRR9694267_2.fastq.gz SKOV3_STAT3KO_rep3_r2.fq.gz
```

<a name="trim"></a>

2. Trimming and aligning the data

Then we use trim galore to trim and generate fastqc reports and then align to the hg38 assembly 


```bash
for f1 in *.fq.gz
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1,2)"
READNUM="$(echo ${f1} | cut -d '_' -f4 | cut -d '.' -f1 | tr -cd '[[:digit:]]')"
if [ $READNUM == 1 ]
then
SAMPLE1="$(echo ${f1})"
echo ${SAMPLE1}
fi
if [ $READNUM == 2 ]
then
echo ${f1}
qsub -S /bin/bash -N trim_${SAMPLE} -cwd -l h_vmem=10G -j y << EOF
module load trim_galore/0.4.1
module load cutadapt/2.8/python.3.7.3
module load FastQC/0.11.4/java.1.8.0_20
trim_galore --fastqc --paired --illumina $SAMPLE1 $f1 
EOF
fi
done

mkdir ../Mapped_STAR_79/

for f1 in *val*.fq.gz
do
SAMPLE="$(echo ${f1} | cut -d '_' -f1,2,3)"
READNUM="$(echo ${f1} | cut -d '_' -f6 | cut -d '.' -f1)"
if [ $READNUM == 1 ]
then
SAMPLE1="$(echo ${f1})"
echo ${SAMPLE1}
fi
echo ${f1} 
if [ $READNUM == 2 ]
then
qsub -S /bin/bash -N align_${SAMPLE} -cwd -l h_vmem=5.6G -j y -pe smp 10 << EOF
module load samtools
module load STAR
STAR --runThreadN 10 --genomeDir /home/greally-lab/indexes/Hg38_rel79_ERCC/STAR/ --readFilesIn ${SAMPLE1} ${f1} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --sjdbGTFfile /home/greally-lab/indexes/Hg38_rel79_ERCC/Homo_sapiens.GRCh38.79.ERCC.gtf --sjdbOverhang 99 --outFileNamePrefix ../Mapped_STAR_79/${SAMPLE}
EOF
fi
done

cd ../Mapped_STAR_79/
```
 
Please see `Analysis-03-23-20` for the downstream R analysis of differentially expressed genes. 