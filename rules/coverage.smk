# bigwig and bedGraph coverage file for future use

rule genomecov_bam:
    input:
        "output/mapped/{sample}-{rep}.merge.sort.bam"
    output:
        "output/coverage/{sample}-{rep, [^.]+}.bedGraph"
    log:
        "logs/genomecov/{sample}-{rep}.log"
    params:
        config["genomecov"]
    conda:
        f"{snake_dir}/envs/deeptools.yaml"
    shell:
        "genomeCoverageBed {params} -ibam {input} | sort -k1,1 -k2,2n 1> {output} 2> {log}"

rule bedGraphToBigWig:
    input:
        bedGraph="output/coverage/{sample}-{rep}.bedGraph",
        chromsizes=config["bedGraphToBigWig"]["chrom"]
    output:
        "output/coverage/{sample}-{rep, [^.]+}.bgToBw.bw"
    log:
        "logs/bedGraphToBigWig/{sample}-{rep}.log"
    params:
        config["bedGraphToBigWig"]["params"]
    conda:
        f"{snake_dir}/envs/deeptools.yaml"
    shell:
        """
        bedGraphToBigWig {params} {input.bedGraph} {input.chromsizes} \
            {output} &> {log}        
        """

rule bamCoverage:
    input:
        bam="output/mapped/{sample}-{rep}.merge.sort.bam",
        bai="output/mapped/{sample}-{rep}.merge.sort.bam.bai"
    output:
        "output/coverage/{sample}-{rep, [^.]+}.bamCov.bw"
    log:
        "logs/bamCoverage/{sample}-{rep}.log"
    params:
        config["bamCoverage"]
    threads:
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=config["mem"]
    conda:
        f"{snake_dir}/envs/deeptools.yaml"
    shell:
        """
        bamCoverage --bam {input.bam} --outFileName {output} --outFileFormat bigwig {params} -p {threads}
        """

if checkcontrol(samples):
    rule bamCompare:
        input:
            ip="output/mapped/{sample}-{rep}.merge.sort.bam",
            ip_index="output/mapped/{sample}-{rep}.merge.sort.bam.bai",
            input=f"output/mapped/{samples.loc[samples['condition']=='control', 'sample'].iloc[0]}-{{rep}}.merge.sort.bam",
            input_index=f"output/mapped/{samples.loc[samples['condition']=='control', 'sample'].iloc[0]}-{{rep}}.merge.sort.bam.bai"
        output:
            "output/coverage/{sample}-{rep}.bamCompare.bw"
        params:
            config['bamCompare']
        log:
            "logs/bamCompare/{sample}-{rep}.log"
        threads:
            config['threads']
        resources:
            cpus=config['threads'],
            mem=config['mem']
        conda:
            f"{snake_dir}/envs/deeptools.yaml"
        shell:
            """
            bamCompare -b1 {input.ip} -b2 {input.input} -o {output} -of bigwig \
                {params} -p {threads}
            """

rule TSS_profile:
    input:
        bw="output/coverage/{prefix}.bamCov.bw",
        region=config["computeMatrix"]["region"]
    output:
        matrix="output/profile/{prefix}.tss2kbp.matrix.gz",
        png=report("output/profile/{prefix}.tss2kbp.matrix.heatmap.png", caption="../report/plotHeatmap.rst", category="deeptools")
    params:
        config["computeMatrix"]["params"]
    log:
        "logs/deeptools/{prefix}.computeMatrix.log"
    threads:
        config["threads"]
    resources:
        cpus=config["threads"],
        mem=config["mem"]
    conda:
        f"{snake_dir}/envs/common.yaml"
    shell:
        """
        computeMatrix reference-point -S {input.bw} -R {input.region} -o {output.matrix} {params} -p {threads}
        plotHeatmap -m {output.matrix} -o {output.png}
        """