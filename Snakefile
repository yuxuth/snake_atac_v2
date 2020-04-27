shell.prefix("set -eo pipefail; echo BEGIN at $(date); ")
shell.suffix("; exitstat=$?; echo END at $(date); echo exit status was $exitstat; exit $exitstat")

configfile: "config.yaml"

FILES = json.load(open(config['SAMPLES_JSON']))

# CLUSTER = json.load(open(config['CLUSTER_JSON']))

SAMPLES = sorted(FILES.keys())
# BWA_INDEX = config['BWA_INDEX']
bowtie2_INDEX = config['bowtie2_INDEX']

TARGETS = []

## map the R1 and R2 seperately. extract the hidIII reads from R1 mapping as the mac2 control. 
## using Mac2 to call peak and building the signaling track

## constructe the target if the inputs are fastqs
ALL_FASTQC_R1  = expand("02_fqc/{sample}_L001_R2_001_fastqc.html", sample = SAMPLES)
bam_r1 = expand("03_aln/{sample}_r1.sorted.bam", sample = SAMPLES)
# bam_r1_updated = expand("03_aln/{sample}_r1_updated.sorted.bam", sample = SAMPLES)
bam_r2 = expand("03_aln/{sample}_r2.sorted.bam", sample = SAMPLES)
ALL_QC = ["10multiQC/multiQC_log.html"]
peak = expand("06_peak_macs2_narrow/{sample}_macs2_peaks.narrowPeak", sample = SAMPLES)
flag_r1 = expand("00_log/{sample}_r1.sorted.bam.flagstat", sample = SAMPLES)
flag_r2 = expand("00_log/{sample}_r2.sorted.bam.flagstat", sample = SAMPLES)
index_r1 = expand("03_aln/{sample}_r1.sorted.bam.bai", sample = SAMPLES)
index_r2 = expand("03_aln/{sample}_r2.sorted.bam.bai", sample = SAMPLES)

"03_aln/{sample}_r1.sorted.bam.bai"


TARGETS.extend(bam_r1) ##append all list to 
TARGETS.extend(bam_r2)
TARGETS.extend(ALL_FASTQC_R1) ## check later
# TARGETS.extend(ALL_QC)
TARGETS.extend(peak)
TARGETS.extend(flag_r1)
TARGETS.extend(flag_r2)


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
# mm for memory efficiency.  -X  fragment size limitation
rule bowtie2_align_r1:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R1']
    output: temp("03_aln/{sample}_r1.sam")
    threads: 16
    message: "bowtie2 {input}: {threads} threads"
    log:
         "00_log/{sample}_r1.bowtie2"
    shell:
        """
        module load bowtie2
        bowtie2  -x {bowtie2_INDEX}  --threads {threads}  -U {input[r1]}  2> {log} \
	| samblaster --removeDups  > {output} ## remove the duplicated
        """

rule bowtie2_align_2:
    input:
        r1 = lambda wildcards: FILES[wildcards.sample]['R2']
    output: temp("03_aln/{sample}_r2.sam")
    threads: 16
    message: "bowtie2 {input}: {threads} threads"
    log:
         "00_log/{sample}_r2.bowtie2"
    shell:
        """
        module load bowtie2
        bowtie2  -x {bowtie2_INDEX}  --threads {threads}  -U {input[r1]}  2> {log} \
	| samblaster --removeDups  > {output}
        """

rule updated_r1_sam:
    input:  "03_aln/{sample}_r1.sam"
    output: sam = temp("03_aln/{sample}_r1.sam_all"),head = temp("03_aln/{sample}_r1.sam_head"), bam = temp("03_aln/{sample}_r1.bam")
    shell:
        """
	module load samtools
        samtools view -q 30 {input}  | awk '$10 ~ /^AGCTT/' > {output[sam]}
	samtools view -H {input} > {output[head]} 
	cat {output[head]} {output[sam]} | samtools -Sb - >  {output[bam]}
        """
	
## remove the unmapped reads 
rule updated_r2_sam:
    input:  "03_aln/{sample}_r2.sam"
    output: temp("03_aln/{sample}_r2.bam")
    shell:
        """
	module load samtools
        samtools view -q 30 -Sb {input} >  {output}
        """

rule sort_bam_r1:
    input: ("03_aln/{sample}_r1.bam")
    output: "03_aln/{sample}_r1.sorted.bam"
    message: "samtools  sort {input} "
    threads : 6
    shell:
        """
        module load samtools
        samtools sort -m 6G -@ {threads} -T {output[0]}.tmp -o {output[0]} {input} 
        """
	
rule sort_bam_r2:
    input: ("03_aln/{sample}_r2.bam")
    output: "03_aln/{sample}_r2.sorted.bam"
    message: "samtools  sort {input} "
    threads : 6
    shell:
        """
        module load samtools
        samtools sort -m 6G -@ {threads} -T {output[0]}.tmp -o {output[0]} {input} 
        """


rule index_bam_r1:
    input:  "03_aln/{sample}_r1.sorted.bam"
    output: "03_aln/{sample}_r1.sorted.bam.bai"
    shell:
        """
        module load samtools
        samtools index {input}
        """
	
rule index_bam_r2:
    input:  "03_aln/{sample}_r2.sorted.bam"
    output: "03_aln/{sample}_r2.sorted.bam.bai"
    shell:
        """
        module load samtools
        samtools index {input} 
        """

rule flagstat_bam_r1:
    input:  "03_aln/{sample}_r1.sorted.bam"
    output: "00_log/{sample}_r1.sorted.bam.flagstat"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 
        """
rule flagstat_bam_r2:
    input:  "03_aln/{sample}_r2.sorted.bam"
    output: "00_log/{sample}_r2.sorted.bam.flagstat"
    threads: 1
    params: jobname = "{sample}"
    message: "flagstat_bam {input}: {threads} threads"
    shell:
        """
        samtools flagstat {input} > {output} 
        """
	
rule call_peaks_macs2_narrow: ## set to large
    input: "03_aln/{sample}_r1.sorted.bam", "03_aln/{sample}_r2.sorted.bam"
    output: bed = "06_peak_macs2_narrow/{sample}_macs2_peaks.narrowPeak"	
    log: "00_log/{sample}_call_broad_peaks_macs2.log"
    params:
        name = "{sample}_macs2",
        jobname = "{sample}"
    message: "call_peaks macs2 narrow {input}: {threads} threads"
    shell:
        """
       module load macs2
       ## for macs2, when nomodel is set, --extsize is default to 200bp, this is the same as 2 * shift-size in macs14.
        macs2 callpeak -t {input[1]} -c  {input[0]} \
            --keep-dup all -f BAM -g hs  --scale-to large  --shift -100 --extsize 200 \
            --outdir 06_peak_macs2_broad -n {params.name} -p 1e-5  -B --SPMR --nomodel &> {log}
        """





