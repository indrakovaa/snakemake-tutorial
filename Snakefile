# Define variable for input function expand
SAMPLES = ["A", "B", "C"]
# Better is to use config file
configfile: "config.yaml" # snakemake stores in dictionary config
# learn to use global wildcards


# Here we can have specified target files (files we want to be created)
# as an input files
rule all:
    input:
        expand("counts/{sample}.fastq.count", sample=SAMPLES), # must create variable SAMPLES
        "plots/quals.svg"

## Workflow
# counting reads
rule count_reads:
  input:  "data/samples/{sample}.fastq"
  output: "counts/{sample}.fastq.count"
  shell:
    "echo $(( $(wc -l <{input}) / 4 )) > {output}"

# mapping reads
rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    threads: 8 #can be usefull to specify, 
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

#sorting reads
rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"

# indexing reads
rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"
       
# variant calling
rule variant_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]), 
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
        #snakemake reads from dictionary config
    output:
        "calls/all.vcf"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"

# using script in plots
rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
