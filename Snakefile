TMPDIR = "/scratch/tereiter"

rule all:
    input:
        "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv"

rule download_CAMI:
    output: "inputs/CAMI_low.tar"
    threads: 1
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    shell:'''
    wget https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100344/ChallengeDatasets.dir/CAMI_low.tar
    '''

rule decompress_CAMI:
    input: "inputs/CAMI_low.tar"
    output: 
        "inputs/CAMI_low/RL_S001__insert_270.fq.gz",
    threads: 1
    resources: 
        mem_mb = 500,
        tmpdir = TMPDIR
    shell:'''
    tar xvf {input} -C inputs/
    '''

rule fastp:
    """
    The CAMI I and CAMI II challenge both showed that "using read quality trimming or error correction software, such as ...Fastp... impoved assembly quality." https://doi.org/10.1101/2021.07.12.451567
    Set minimum read length to megahit default minimum kmer length, k = 21
    """
    input: "inputs/CAMI_low/RL_S001__insert_270.fq.gz"
    output: 
        r1 = "outputs/fastp/CAMI_low_R1.fastp.fq.gz",
        r2 = "outputs/fastp/CAMI_low_R2.fastp.fq.gz",
        json = "outputs/fastp/CAMI_low.json",
        html = "outputs/fastp/CAMI_low.html"
    resources: 
        mem_mb = 8000,
        tmpdir=TMPDIR
    threads: 8
    benchmark: "benchmarks/fastp_CAMI_low.txt"
    conda: "envs/fastp.yml"
    shell:'''
    fastp -i {input} --interleaved_in -o {output.r1} -O {output.r2} -q 4 -j {output.json} -h {output.html} -R CAMI_low -l 21 -c -w {threads}
    '''

rule abundtrim:
    input:
        r1 = "outputs/fastp/CAMI_low_R1.fastp.fq.gz",
        r2 = "outputs/fastp/CAMI_low_R2.fastp.fq.gz",
    output: "outputs/abundtrim/CAMI_low.abundtrim.fq.gz"
    resources: 
        mem_mb = 64000,
        tmpdir=TMPDIR
    threads: 1
    conda: "envs/khmer.yml"
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule spacegraphcats:
    input:
        fq="outputs/abundtrim/CAMI_low.abundtrim.fq.gz",
        conf="conf/CAMI_low_sgc_conf1.yml"
    output: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv"
    resources: 
        mem_mb = 900000,
        tmpdir=TMPDIR
    threads: 1
    params: outdir = "outputs/spacegraphcats/"
    conda: "envs/spacegraphcats_multifasta_prot.yml"
    shell:'''
    python -m spacegraphcats run {input.conf} build_cdbg_list_by_record_x  --nolock --outdir={params.outdir} --rerun-incomplete 
    '''
