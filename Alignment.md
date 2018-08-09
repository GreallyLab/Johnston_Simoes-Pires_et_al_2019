# Alignment of RNAseq reads
### Andrew D. Johnston

Directional RNAseq libraries from samples from each experiment (6 in each) were single-end sequenced to 100bp on 1 lane of an Illumina HiSeq2500. The reads were trimmed using trim_galore using default settings (not shown here). The trimmed reads were aligned using STAR to an Ensembl release 79 reference genome plus ERCC spike-in controls.

The reference genome was generated as such:
```bash
mkdir Hg38_rel79_ERCC
cat Hg38_rel79/Homo_sapiens.GRCh38.dna.primary_assembly.fa ERCC92/ERCC92.fa > Hg38_rel79_ERCC/Homo_sapiens.GRCh38.rel79.cdna.all.ERCC.fa
cat hg38/ensembl_79/Homo_sapiens.GRCh38.79.gtf ERCC92/ERCC92.gtf > Hg38_rel79_ERCC/Homo_sapiens.GRCh38.79.ERCC.gtf

 # the genome build for STAR was made.
qsub -S /bin/bash -N Build_Star_rel79 -j y -cwd -q highmem.q -pe smp-highmem 20 -l h_vmem=10G << EOF
module load samtools
module load STAR
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/greally-lab/indexes/Hg38_rel79_ERCC --genomeFastaFiles /home/greally-lab/indexes/Hg38_rel79_ERCC/Homo_sapiens.GRCh38.rel79.cdna.all.ERCC.fa
EOF
```

Next, the trimmed reads were aligned to this genome. 

```bash
for f1 in *_trimmed.fq.gz;
do 
SAMPLE="$(echo ${f1} | cut -d '_' -f 5,6)"
echo $SAMPLE
echo $f1
qsub -S /bin/bash -N C_${SAMPLE}_STAR_79 -cwd -l h_vmem=5.6G -j y -pe smp 20 << EOF
module load samtools
module load STAR
STAR --runThreadN 20 --genomeDir /home/greally-lab/indexes/Hg38_rel79_ERCC/STAR/ --readFilesIn ${f1} --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --sjdbGTFfile /home/greally-lab/indexes/Hg38_rel79_ERCC/Homo_sapiens.GRCh38.79.ERCC.gtf --sjdbOverhang 99 --outFileNamePrefix Mapped_STAR_79_ERCC/${SAMPLE}
EOF
done 
```

The generated count tables were imported into R for downstream anaylsis. Please see <a href="DEG_analysis.Rmd">DEG_Analysis.Rmd</a>