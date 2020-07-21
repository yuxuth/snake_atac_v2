shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
bowtie2_INDEX = config['bowtie2_INDEX']

TARGETS = []

## test


## constructe the target if the inputs are fastqs
# ALL_FASTQC  = expand("02_fqc/{sample}_L001_R2_001_fastqc.html", sample = SAMPLES)

bam = expand("03_aln/{sample}.sorted.bam", sample = SAMPLES)

# ALL_QC = ["10multiQC/multiQC_log.html"]
# peak = expand("06_peak_macs2_broad/{sample}_macs2_peaks.narrowPeak", sample = SAMPLES)
flag = expand("00_log/{sample}.sorted.bam.flagstat", sample = SAMPLES)
bam_index = expand("03_aln/{sample}.sorted.bam.bai", sample = SAMPLES)
TARGETS.extend(bam) ##append all list to 
TARGETS.extend(flag)
TARGETS.extend(bam_index)


localrules: all

rule all:
    input: TARGETS




rule bowtie2_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'], #,
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: temp("03_aln/{sample}.sam")
    threads: 12
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bowtie2"
    shell:
        """
        module load bowtie2
        bowtie2  -x {bowtie2_INDEX}  -X 2000 --threads {threads}  -1 {input[0]}  -2 {input[1]} > {output}  2> {log}
        """


rule remove_duplicate:
    input:  "03_aln/{sample}.sam"
    output: temp("03_aln/{sample}.duremovedsam")
    message: "samblaster {input} "
    log: "00_log/{sample}.dup_removed"
    shell:
        """
        samblaster --removeDups -i {input} -o {output} 2>  {log}
        """

rule sam_to_bam:
    input:  ("03_aln/{sample}.duremovedsam")
    output: ( "03_aln/{sample}.tmp.bam")
    message: "samtools {input}: to sam"
    shell:
        """
        module load samtools
        samtools view -Sb -F 4 {input}  > {output}
        """

rule sort_bam:
    input: ("03_aln/{sample}.tmp.bam")
    output: "03_aln/{sample}.sorted.bam"
    message: "samtools  sort {input} "
    threads : 6
    shell:
        """
        module load samtools
        samtools sort -m 6G -@ {threads} -T {output[0]}.tmp -o {output[0]} {input} 
        """

rule index_bam:
    input:  "03_aln/{sample}.sorted.bam"
    output: "03_aln/{sample}.sorted.bam.bai"
    log:    "00_log/{sample}.index_bam"
    message: "index_bam {input}: {threads} threads"
    shell:
        """
        module load samtools
        samtools index {input} 2> {log}
        """

# check number of reads mapped by samtools flagstat, the output will be used for downsampling
rule flagstat_bam:
    input:  "03_aln/{sample}.sorted.bam"
    output: "00_log/{sample}.sorted.bam.flagstat"
    log:    "00_log/{sample}.flagstat_bam"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        module load samtools
        samtools flagstat {input} > {output} 2> {log}
        """

