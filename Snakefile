shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
# BWA_INDEX = config['BWA_INDEX']
bowtie2_INDEX = config['bowtie2_INDEX']

TARGETS = []

## test


## constructe the target if the inputs are fastqs
ALL_FASTQC  = expand("02_fqc/{sample}_L001_R2_001_fastqc.html", sample = SAMPLES)

bam = expand("03_aln/{sample}.sorted.bam", sample = SAMPLES)

ALL_QC = ["10multiQC/multiQC_log.html"]
peak = expand("06_peak_macs2_broad/{sample}_macs2_peaks.narrowPeak", sample = SAMPLES)
flag = expand("00_log/{sample}.sorted.bam.flagstat", sample = SAMPLES)

TARGETS.extend(bam) ##append all list to 
# TARGETS.extend(ALL_FASTQC) ## check later
# TARGETS.extend(ALL_QC)
# TARGETS.extend(peak)
TARGETS.extend(flag)


localrules: all

rule all:
    input: TARGETS


rule fastqc:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: "02_fqc/{sample}_L001_R2_001_fastqc.html" 
    threads: 1
    params : jobname = "{sample}"
    message: "fastqc {input}: {threads}"
    log:   "00_log/{sample}_fastqc"
    shell:
        """
        module load fastqc
        fastqc -o 02_fqc -f fastq --noextract {input}  2> {log}
        """

rule trim_fastqs: ## merge fastq
	input:
		r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
		r2 = lambda wildcards: FILES[wildcards.sample]['R2']
	output:
		("01_trim_seq/{sample}_1.fastq" ), ("01_trim_seq/{sample}_2.fastq")
	log: "00_log/{sample}_trim_adapter.log"
	params:
		jobname = "{sample}"
	threads : 8
	# group: "mygroup"
	message: "trim fastqs {input}: {threads} threads"
	shell:
		"""
		NGmerge  -a -y  -n {threads} -1 {input[0]} -2 {input[1]}  -o 01_trim_seq/{params.jobname} 2> {log} 
		"""
        
# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort
# mm for memory efficiency.  -X  fragment size limitation
rule bowtie2_align:
    input:
        "01_trim_seq/{sample}_1.fastq", 
        "01_trim_seq/{sample}_2.fastq"
    output: temp("03_aln/{sample}.sam")
    threads: 6
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bowtie2"
    shell:
        """
        module load bowtie2
        bowtie2 --mm -x {bowtie2_INDEX}  -X 2000 --threads {threads}  -1 {input[0]}  -2 {input[1]} > {output}  2> {log}
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

## remove the unmapped reads 
rule sam_to_bam:
    input:  ("03_aln/{sample}.duremovedsam")
    output: temp( "03_aln/{sample}.tmp.bam")
    message: "samtools {input}: to sam"
    shell:
        """
        module load samtools
        samtools view -Sb -F 4  -q 30 {input}  > {output}
        """

rule sort_bam:
    input: ("03_aln/{sample}.tmp.bam")
    output: "03_aln/{sample}.sorted.bam"
    message: "samtools  sort {input} "
    threads : 12
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






