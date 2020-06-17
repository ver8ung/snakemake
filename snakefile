## define rules for snakemake to execute
SAMPLES = ["A", "B"]                                       ## define input samples

rule all:                                                  ## specify desired target file
    input:
        "plots/quals.svg"

rule bwa_map:                                               ## using BWA to map reads to a reference genome
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"
        
rule samtools_sort:                                         ## using SAMtools to sort the mapped reads
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

rule samtools_index:                                       ## indexing reads for random access
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"

rule bcftools_call:                                       ## call variants using BCFtools
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
  
  rule plot_quals:                                        ## rule for plotting quality score invoking another python script
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/qscore-plot.py"
