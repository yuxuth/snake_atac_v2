shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
BWA_INDEX = config['BWA_INDEX']

TARGETS = []

## test


## constructe the target if the inputs are fastqs
ALL_FASTQC  = expand("02_fqc/{sample}_L001_R2_001_fastqc.html", sample = SAMPLES)

bam = expand("03_aln/{sample}.sorted.bam", sample = SAMPLES)

ALL_QC = ["10multiQC/multiQC_log.html"]
peak = expand("06_peak_macs2_broad/{sample}_macs2_peaks.narrowPeak", sample = SAMPLES)
flag = expand("00_log/{sample}.sorted.bam.flagstat", sample = SAMPLES)

TARGETS.extend(bam) ##append all list to 
TARGETS.extend(ALL_FASTQC) ## check later
TARGETS.extend(ALL_QC)
TARGETS.extend(peak)
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

# get the duplicates marked sorted bam, remove unmapped reads by samtools view -F 4 and dupliated reads by samblaster -r
# samblaster should run before samtools sort

rule bwa_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        r2 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: temp("03_aln/{sample}.sam")
    threads: 12
    message: "bwa {input}: {threads} threads"
    log:
         "00_log/{sample}.bwa"
    shell:
        """
        module load bwa
        bwa mem  -t {threads} {BWA_INDEX} {input}  > {output}  2> {log}
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
    output: temp( "03_aln/{sample}.tmp.bam")
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

rule call_peaks_macs2_narrow:
    input: "03_aln/{sample}.sorted.bam", "03_aln/{sample}.sorted.bam.bai"
    output: bed = "06_peak_macs2_broad/{sample}_macs2_peaks.narrowPeak"  ## case control is defined by the output 
    log: "00_log/{sample}_call_broad_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}"
    message: "call_peaks macs2 narrow {input}: {threads} threads"
    shell:
        """
       module load macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[0]} \
            --keep-dup all -f BAM -g hs --shift 75  --extsize 150  \
            --outdir 06_peak_macs2_broad -n {params.name} -p 0.01  -B --SPMR --nomodel &> {log}
        """


rule multiQC:
    input :
        expand("00_log/{sample}.sorted.bam.flagstat", sample = SAMPLES),
        expand("02_fqc/{sample}_L001_R2_001_fastqc.html", sample = SAMPLES)
    output: "10multiQC/multiQC_log.html"
    log: "00log/multiqc.log"
    message: "multiqc for all logs"
    shell:
        """
        multiqc 02_fqc 00_log -o 10multiQC -d -f -v -n multiQC_log 2> {log}
        """

