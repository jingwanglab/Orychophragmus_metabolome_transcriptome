
configfile: "config.yaml"

SAMTOOLS=config["SAMTOOLS"]
HISAT_DIR=config["HISAT_DIR"]
FASTP=config["FASTP"]
ALL_SAMPLE = [key for key, value in config["inputFiles"].items()]
genome=config["genome"]
gtf=config["gtf"]
RSCRIPT=config["Rscript"] 


rule all:
    input:
        ["0.fastp/"+str(sample)+"r2.trimmed.fq.gz" for sample in ALL_SAMPLE],
        ["2.align/"+str(sample)+".sorted.bam" for sample in ALL_SAMPLE],
        #["1.index/genome.8.ht2"],
        expand("1.index/genome.{num}.ht2", num=range(1,8)),
        #bam = "2.align/{sample}.sorted.bam
        ["2.align/"+str(sample)+"_align.log" for sample in ALL_SAMPLE],
        ["3.htseq/"+str(sample)+".counts" for sample in ALL_SAMPLE],
        ["4.results/Ov_tpm.txt"]
        #count = "3.htseq/{sample}.counts"


rule fastp:
    input:
        #r1 = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['f1'],
        r1 = lambda w: config["inputFiles"][w.sample]['f1'],
        r2 = lambda w: config["inputFiles"][w.sample]['r1']
        #r2 = lambda w: config["inputFiles"]['{}'.format(w.SPC)]['r1']
    output:
        r1t = "0.fastp/{sample}r1.trimmed.fq.gz",
        r2t = "0.fastp/{sample}r2.trimmed.fq.gz",
        json = "0.fastp/{sample}.json",
        html = "0.fastp/{sample}.html"
    #log:
    #    "0.fastp/{sample}.log"
    threads: 8 
    shell:
        """
        {FASTP} -i {input.r1} -I {input.r2} -o {output.r1t} -O {output.r2t} --length_required=35 --cut_mean_quality 20 -5 20 -3 20 -w {threads} -j {output.json} -h {output.html} 2>/dev/null
        """

rule hisat_index:
    input:
        r2t = ["0.fastp/"+str(sample)+"r2.trimmed.fq.gz" for sample in ALL_SAMPLE],
        r1t = ["0.fastp/"+str(sample)+"r1.trimmed.fq.gz" for sample in ALL_SAMPLE]
    output:
        #ind=["1.index/{g.ht2"]
        expand("1.index/genome.{num}.ht2", num=range(1,8))
    threads: 20
    shell:
        """
        {HISAT_DIR}/hisat2_extract_splice_sites.py {gtf} > 1.index/ss.txt
        {HISAT_DIR}/hisat2_extract_exons.py  {gtf} >  1.index/exon.txt
        {HISAT_DIR}/hisat2-build -p {threads} --ss 1.index/ss.txt --exon 1.index/exon.txt {genome} 1.index/genome
        """


rule hisat_align:
    input:
        r1t = "0.fastp/{sample}r1.trimmed.fq.gz",
        r2t = "0.fastp/{sample}r2.trimmed.fq.gz",
        ind = expand("1.index/genome.{num}.ht2", num=range(1,8))
        #genome=genome
    output:
        bam = "2.align/{sample}.sorted.bam"
    threads: 5
    log:
        "2.align/{sample}_align.log"
    shell:
        """
        {HISAT_DIR}/hisat2 -p {threads} -x 1.index/genome -1 {input.r1t} -2 {input.r2t} -k1 2>{log} | {SAMTOOLS} sort -O BAM  -@ 4 > {output.bam}
        {SAMTOOLS} index {output.bam}
        """

rule htseq:
    input:
        bam = "2.align/{sample}.sorted.bam"
    output:
        counts = "3.htseq/{sample}.counts"
    threads: 1
    shell:    
        """
        htseq-count --mode=union --nonunique=none -s no --secondary-alignments=ignore  -f bam  {input.bam} {gtf} > {output.counts}
        """

rule call:
    input:
        #counts = "3.htseq/"
        counts = ["3.htseq/"+str(sample)+".counts" for sample in ALL_SAMPLE]
        #counts = "3.htseq/"+ALL_SAMPLE[-1]+".counts"
    output:
        tpm = "4.results/Ov_tpm.txt",
        #count= "4.results/Ov_RNA.counts"
    shell:
        """
        python scripts/count_matrix.py {gtf} Ov 3.htseq/ 4.results/Ov_RNA.counts
        {RSCRIPT} scripts/tpm_calculation.R 4.results/Ov_RNA.counts 4.results/Ov_fpkm.txt {output.tpm}
        """

#        python scripts/count_matrix.py {gtf} Ov 3.htseq/ 4.results/Ov_RNA.counts





