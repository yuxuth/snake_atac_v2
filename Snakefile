shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
bowtie2_INDEX = config['bowtie2_INDEX']

TARGETS = []
chromsize = '~/data/Index/hg19.all.chrom.sizes'
## test

# bw = ["07_bw/CHIP-TBID-LTR7-2_vs_CHIP-TBID-EV-2_control_lambda.bw",
# "07_bw/CHIP-TBID-LTR7-2_vs_CHIP-TBID-EV-2_treat_pileup.bw"]


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

# # macs_out =[ "06_peak_macs2_narrow/CHIP-TBID-LTR7-2/CHIP-TBID-LTR7-2_vs_CHIP-TBID-EV-2_control_lambda.bdg",
# "06_peak_macs2_narrow/CHIP-TBID-LTR7-2/CHIP-TBID-LTR7-2_vs_CHIP-TBID-EV-2_treat_pileup.bdg",
# "06_peak_macs2_narrow/CHIP-TBID-LTR7-2/CHIP-TBID-LTR7-2_vs_CHIP-TBID-EV-2_peaks.narrowPeak"]
id = ["CHIP-H9_EV_H3K27ac-1", "CHIP-H9_EV_H3K9me3-1", 
"CHIP-H9_sgYAP-_H3K27ac-1", "CHIP-H9_sgYAP-_H3K27ac-2", 
"CHIP-H9_sgYAP-H3K9me3-2", "CHIP-YAP-1_abcam", "CHIP-YAP-1_nb",
 "CHIP-YAP-2_abcam", "CHIP-YAP-2_nb"]
id2 = [ '_control_lambda.bdg', '_treat_pileup.bdg', '_peaks.narrowPeak']
bdg = expand("06_peak_macs2_narrow/{id}/{id}_vs_input-1{id2}", id=id, id2=id2)

id3 = [ '_control_lambda.bw', '_treat_pileup.bw']
bw= expand("07_bw/{id}_vs_input-1{id3}", id=id, id3=id3)

# print(TARGETS)

localrules: all

rule all:
    input: bdg,bw




rule bowtie2_align:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1'], 
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


## using macs2 to call peaks 

rule call_peaks_macs2_narrow:
    input:  igg = "03_aln/{igg}.sorted.bam", 
            sample = "03_aln/{sample}.sorted.bam"
    output: "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_control_lambda.bdg",
            "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_treat_pileup.bdg",
            "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_peaks.narrowPeak",

    log: "00_log/macs2_{sample}_vs_{igg}_call_narrow_peaks_macs2.log"
    params:
        name = "{sample}_vs_{igg}",
        jobname = "{sample}"
    message: "call_peaks macs2 {input}"
    shell:
        """
       module load macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input.sample} \
            -c {input.igg} --keep-dup all -f BAM -g hs \
            --outdir 06_peak_macs2_narrow/{params.jobname} -n {params.name} -p 1e-2 --SPMR -B --nomodel &> {log}
        """


## convert bdg to bw

rule macs2_bdg_convert1:
    input: "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_control_lambda.bdg"
            # "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_treat_pileup.bdg"
    output: "07_bw/{sample}_vs_{igg}_control_lambda.bw"
            # "07_bw/{sample}_vs_{igg}_treat_pileup.bw",
    message: "bw convert {input}"
    shell:
        """
        bdg2bw  {input} {chromsize}
        mv 06_peak_macs2_narrow/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.igg}_control_lambda.bw   {output}
        """

rule macs2_bdg_convert2:
    input: 
            "06_peak_macs2_narrow/{sample}/{sample}_vs_{igg}_treat_pileup.bdg"
    output: 
            "07_bw/{sample}_vs_{igg}_treat_pileup.bw"
    message: "bw convert {input}"
    shell:
        """
        bdg2bw  {input} {chromsize}
        mv 06_peak_macs2_narrow/{wildcards.sample}/{wildcards.sample}_vs_{wildcards.igg}_treat_pileup.bw   {output}
        """
