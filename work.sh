zcat 23022FL-12-01-09_S9_L001_R2_001.fastq.gz | fastx_clipper -a GATCGTCGGACTGTAGAACTCTGAAC -l 30 -d 1 -c -n -v -M 26 -o fastx.fastq
umi_tools extract --stdin=umi_test.fastq --bc-pattern=NNNNNNXXXXXXXXXXXXXXXXXXXXXXXXXXX --log=processed.log --stdout processed.fastq.gz --3prime
