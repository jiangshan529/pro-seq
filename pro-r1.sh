## use R1
prefix=PRO6_S5_L008_R1_001
hg38chrominfo=/broad/boxialab/shawn/projects/genome_sequence/human/hg38/hg38.clean.chrom.sizes

# Create a directory for the prefix
mkdir -p ${prefix}
cutadapt -a TGGAATTCTCGGGTGCCAAGG -o ${prefix}/${prefix}_cut.fastq.gz -m 34 -j 8 ${prefix}.fastq.gz
zcat ${prefix}/${prefix}_cut.fastq.gz | fastx_reverse_complement -o ${prefix}/${prefix}_cut_reverse.fastq.gz -z
umi_tools extract --stdin=${prefix}/${prefix}_cut_reverse.fastq.gz --bc-pattern=NNNNNNC --log=${prefix}/${prefix}_cut_reverse.fastq.log --stdout ${prefix}/${prefix}_trim_umi.fastq.gz
cutadapt -u -7 -o ${prefix}/${prefix}_trim_umi_cut.fastq.gz ${prefix}/${prefix}_trim_umi.fastq.gz

# Step 4: Bowtie2 and Samtools (mapping and bam creation)
bowtie2 -p 16 --end-to-end --sensitive --no-unal -U ${prefix}/${prefix}_trim_umi_cut.fastq.gz -x /seq/epiprod02/kdong/SofiaSandbox/NanoNOMe/hg38.chrXYM 2>${prefix}/${prefix}_bowtie2.log | samtools view -@ 16 -Sb -q 20 -o ${prefix}/${prefix}_trim.bam

# Step 5: Sorting BAM
samtools sort -@ 8 ${prefix}/${prefix}_trim.bam -o ${prefix}/${prefix}_trim_sort.bam

# Step 6: Indexing BAM
samtools index ${prefix}/${prefix}_trim_sort.bam

# Step 7: UMI Tools deduplication
umi_tools dedup -I ${prefix}/${prefix}_trim_sort.bam --output-stats=${prefix}/${prefix}_trim_sort_dedup.txt -S ${prefix}/${prefix}_trim_sort_dedup.bam > ${prefix}/${prefix}_dedup.log 2>&1

# Step 8: Convert BAM to BED
bedtools bamtobed -i ${prefix}/${prefix}_trim_sort_dedup.bam > ${prefix}/${prefix}_trim_sort_dedup.bed

# Step 9: Generate bedgraph files for positive strand
awk '$6 == "+"' ${prefix}/${prefix}_trim_sort_dedup.bed | genomeCoverageBed -i stdin -3 -bg -g $hg38chrominfo > ${prefix}/${prefix}_trim_sort_dedup_hg38_plus.bedgraph

# Step 10: Generate bedgraph files for negative strand
awk '$6 == "-"' ${prefix}/${prefix}_trim_sort_dedup.bed | genomeCoverageBed -i stdin -3 -bg -g $hg38chrominfo > ${prefix}/${prefix}_trim_sort_dedup_hg38_m.bedgraph

# Step 11: Multiply values in the negative strand bedgraph by -1
awk '{$4=$4*-1; print}' ${prefix}/${prefix}_trim_sort_dedup_hg38_m.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_hg38_minus.bedgraph

# Step 12: Convert bedgraph to BigWig for positive strand
sort -k1,1 -k2,2n ${prefix}/${prefix}_trim_sort_dedup_hg38_plus.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_hg38_plus_sort.bedgraph
bedGraphToBigWig  ${prefix}/${prefix}_trim_sort_dedup_hg38_plus_sort.bedgraph $hg38chrominfo  ${prefix}/${prefix}_trim_sort_dedup_hg38_plus_sort.bw

# Step 13: Convert bedgraph to BigWig for negative strand
sort -k1,1 -k2,2n ${prefix}/${prefix}_trim_sort_dedup_hg38_minus.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_hg38_minus_sort.bedgraph
bedGraphToBigWig ${prefix}/${prefix}_trim_sort_dedup_hg38_minus_sort.bedgraph $hg38chrominfo  ${prefix}/${prefix}_trim_sort_dedup_hg38_minus_sort.bw

# Using genomecov for further bedgraph creation
# Step 14: Generate genome coverage bedgraph for positive strand
bedtools genomecov -ibam ${prefix}/${prefix}_trim_sort_dedup.bam -bg -strand + > ${prefix}/${prefix}_trim_sort_dedup_plus_genomecov.bedgraph

# Step 15: Generate genome coverage bedgraph for negative strand
bedtools genomecov -ibam ${prefix}/${prefix}_trim_sort_dedup.bam -bg -strand - > ${prefix}/${prefix}_trim_sort_dedup_m_genomecov.bedgraph

# Step 16: Multiply values in the negative strand genomecov bedgraph by -1
awk '{$4=$4*-1; print}' ${prefix}/${prefix}_trim_sort_dedup_m_genomecov.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_minus_genomecov.bedgraph

# Step 17: Sort positive strand bedgraph by coordinates
sort -k1,1 -k2,2n ${prefix}/${prefix}_trim_sort_dedup_plus_genomecov.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_plus_genomecov_SORT.bedgraph

# Step 18: Sort negative strand bedgraph by coordinates
sort -k1,1 -k2,2n ${prefix}/${prefix}_trim_sort_dedup_minus_genomecov.bedgraph > ${prefix}/${prefix}_trim_sort_dedup_minus_genomecov_SORT.bedgraph

# Step 19: Convert sorted positive strand genomecov bedgraph to BigWig
bedGraphToBigWig ${prefix}/${prefix}_trim_sort_dedup_plus_genomecov_SORT.bedgraph $hg38chrominfo ${prefix}/${prefix}_trim_sort_dedup_plus_genomecov.bw

# Step 20: Convert sorted negative strand genomecov bedgraph to BigWig
bedGraphToBigWig ${prefix}/${prefix}_trim_sort_dedup_minus_genomecov_SORT.bedgraph $hg38chrominfo ${prefix}/${prefix}_trim_sort_dedup_minus_genomecov.bw

