prefix=22KHMTLT3_8_0420620629_PRO5_S72_L008_R2_001
hg38chrominfo=/broad/boxialab/shawn/projects/genome_sequence/human/hg38/hg38.clean.chrom.sizes

zcat ${prefix}.fastq.gz | fastx_clipper -a GATCGTCGGACTGTAGAACTCTGAAC -l 34 -d 0 -c -n -v -M 26 -o ${prefix}.fastq
umi_tools extract --stdin=${prefix}.fastq --bc-pattern=NNNNNN --log=${prefix}.log --stdout ${prefix}.fastq.gz --3prime
cutadapt -u 7 -o ${prefix}_trim.fastq.gz ${prefix}.fastq.gz
bowtie2 -p 16 --end-to-end --sensitive --no-unal -U ${prefix}_trim.fastq.gz -x /seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM 2>${prefix}_bowtie2.log | samtools view -@ 16 -Sb -q 20 -o ${prefix}_trim.bam
samtools sort -@ 8 ${prefix}_trim.bam -o ${prefix}_trim_sort.bam
samtools index ${prefix}_trim_sort.bam
umi_tools dedup -I ${prefix}_trim_sort.bam --output-stats=${prefix}_trim_sort_dedup.txt -S ${prefix}_trim_sort_dedup.bam  > ${prefix}_dedup.log 2>&1
bedtools bamtobed -i ${prefix}_trim_sort_dedup.bam > ${prefix}_trim_sort_dedup.bed
awk '$6 == "+"' ${prefix}_trim_sort_dedup.bed | genomeCoverageBed -i stdin -3 -bg -g $hg38chrominfo > ${prefix}_trim_sort_dedup_hg38_plus.bedgraph
awk '$6 == "-"' ${prefix}_trim_sort_dedup.bed | genomeCoverageBed -i stdin -3 -bg -g $hg38chrominfo > ${prefix}_trim_sort_dedup_hg38_m.bedgraph
awk '{$4=$4*-1; print}' ${prefix}_trim_sort_dedup_hg38_m.bedgraph > ${prefix}_trim_sort_dedup_hg38_minus.bedgraph
bedGraphToBigWig ${prefix}_trim_sort_dedup_hg38_plus.bedgraph $hg38chrominfo ${prefix}_trim_sort_dedup_hg38_plus.bw
bedGraphToBigWig ${prefix}_trim_sort_dedup_hg38_minus.bedgraph $hg38chrominfo ${prefix}_trim_sort_dedup_hg38_minus.bw

### using genomecov 
bedtools genomecov -ibam ${prefix}_trim_sort_dedup.bam -bg -strand + > ${prefix}_trim_sort_dedup_plus_genomecov.bedgraph
bedtools genomecov -ibam ${prefix}_trim_sort_dedup.bam -bg -strand - > ${prefix}_trim_sort_dedup_m_genomecov.bedgraph
awk '{$4=$4*-1; print}' ${prefix}_trim_sort_dedup_m_genomecov.bedgraph > ${prefix}_trim_sort_dedup_minus_genomecov.bedgraph
sort -k1,1 -k2,2n ${prefix}_trim_sort_dedup_plus_genomecov.bedgraph > ${prefix}_trim_sort_dedup_plus_genomecov_SROT.bedgraph
sort -k1,1 -k2,2n ${prefix}_trim_sort_dedup_minus_genomecov.bedgraph > ${prefix}_trim_sort_dedup_minus_genomecov_SORT.bedgraph
bedGraphToBigWig ${prefix}_trim_sort_dedup_plus_genomecov_SROT.bedgraph $hg38chrominfo ${prefix}_trim_sort_dedup_plus_genomecov.bw
bedGraphToBigWig ${prefix}_trim_sort_dedup_minus_genomecov_SORT.bedgraph $hg38chrominfo ${prefix}_trim_sort_dedup_minus_genomecov.bw
