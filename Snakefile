TMPDIR = "/scratch/tereiter"

def checkpoint_separate_cdbg_nodes_by_annot(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.separate_cdbg_nodes_by_annot.get(**wildcards).output[0]    
    file_names = expand("outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.reads.gz",
                        pfam = glob_wildcards(os.path.join(checkpoint_output, "{pfam}_cdbg_nodes.tsv.gz")).pfam)
    return file_names

def checkpoint_separate_cdbg_nodes_by_annot2(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.separate_cdbg_nodes_by_annot.get(**wildcards).output[0]    
    file_names = expand("outputs/gs_read_sets/{pfam}_gs_read_names.tsv",
                        pfam = glob_wildcards(os.path.join(checkpoint_output, "{pfam}_cdbg_nodes.tsv.gz")).pfam)
    return file_names

rule all:
    input:
        "outputs/abundtrim/reads_that_were_trimmed_before_sgc_names.txt",
        checkpoint_separate_cdbg_nodes_by_annot,
        checkpoint_separate_cdbg_nodes_by_annot2

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

rule grab_trimmed_reads:
    input: 
        raw = "inputs/CAMI_low/RL_S001__insert_270.fq.gz",
        abundtrim = "outputs/abundtrim/CAMI_low.abundtrim.fq.gz"
    output: "outputs/abundtrim/reads_that_were_trimmed_before_sgc.fq"
    resources: 
        mem_mb = 16000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    gunzip -c {input.raw} {input.abundtrim} | paste  - - - - | sort |uniq |tr "\t" "\n" > {output}
    '''

rule grab_names_of_trimmed_reads:
    input: "outputs/abundtrim/reads_that_were_trimmed_before_sgc.fq"
    output: "outputs/abundtrim/reads_that_were_trimmed_before_sgc_names.txt" 
    resources: 
        mem_mb = 16000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    grep "@" {input} > {output}
    '''

rule spacegraphcats:
    input:
        fq="outputs/abundtrim/CAMI_low.abundtrim.fq.gz",
        conf="conf/CAMI_low_sgc_conf2.yml"
    output: 
        "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv",
        "outputs/spacegraphcats/CAMI_low_k31/bcalm.unitigs.db"
    resources: 
        mem_mb = 900000,
        tmpdir=TMPDIR
    threads: 1
    benchmark: "benchmarks/sgc_cami_low_k31_r1_mutlifasta_x_gather_nbhd_all_pfam.tsv"
    params: outdir = "outputs/spacegraphcats/"
    conda: "envs/spacegraphcats_prot_gather.yml"
    shell:'''
    python -m spacegraphcats run {input.conf} build_cdbg_list_by_record_x  --nolock --outdir={params.outdir} --rerun-incomplete 
    '''

checkpoint separate_cdbg_nodes_by_annot:
# right now separate by filename; may change with new input format
    input: cdbg_annot = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x/multifasta_x.cdbg_annot.csv"
    output: directory("outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences")
    params: outdir = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/"
    conda: "envs/tidyverse.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_separate_multifasta_results_by_pfam.tsv"
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/separate_multifasta_results_by_pfam.R"
    
rule promote_multifasta_cdbg_nodes_to_neighborhood:
    input: 
        cdbg_nodes = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}_cdbg_nodes.tsv.gz"
    output: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.gz"
    params:
        cdbg_dir="outputs/spacegraphcats/CAMI_low_k31",
        catlas_dir = "outputs/spacegraphcats/CAMI_low_k31_r1",
    conda: "envs/spacegraphcats_prot_gather.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_promote_multifasta_cdbg_nodes_to_nbhd_{pfam}.tsv"
    resources:
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    python -m spacegraphcats.search.extract_neighborhoods_by_cdbg_ids {params.cdbg_dir} {params.catlas_dir} {input.cdbg_nodes} -o {output}
    '''

rule extract_contig_sequences:
    input:
        contigs_db = "outputs/spacegraphcats/CAMI_low_k31/bcalm.unitigs.db",
        cdbg_nbhds = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.gz"
    output: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.fa"
    conda: "envs/spacegraphcats_prot_gather.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_promote_multifasta_cdbg_nodes_to_nbhd_{pfam}.tsv"
    resources:
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    python -m spacegraphcats.search.extract_contigs --contigs-db {input.contigs_db} {input.cdbg_nbhds} -o {output}
    '''

rule make_reads_bgz:
    input:
        reads = "outputs/abundtrim/CAMI_low.abundtrim.fq.gz"
    output: bgz = "outputs/spacegraphcats/CAMI_low/reads.bgz"
    conda: "envs/spacegraphcats_prot_gather.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_promote_reads_bgz.tsv"
    resources:
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    python -Werror -m spacegraphcats.utils.make_bgzf {input.reads} -o {output.bgz}
    '''

rule index_reads:
    input: bgz = "outputs/spacegraphcats/CAMI_low/reads.bgz"
    output: idx = "outputs/spacegraphcats/CAMI_low_k31/reads.bgz.index"
    conda: "envs/spacegraphcats_prot_gather.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_index_reads.tsv"
    resources:
        mem_mb = 200000,
        tmpdir = TMPDIR
    threads: 1
    params: cdbg_dir = "outputs/spacegraphcats/CAMI_low_k31"
    shell:'''
    python -Werror -m spacegraphcats.cdbg.index_reads {params.cdbg_dir} {input.bgz} {output.idx} --expect-paired -k 31
    '''

rule extract_reads:
    input: 
        bgz ="outputs/spacegraphcats/CAMI_low/reads.bgz",
        idx = "outputs/spacegraphcats/CAMI_low_k31/reads.bgz.index",
        cdbg_nbhds = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.gz"
    output: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.reads.gz"
    threads: 1
    conda: "envs/spacegraphcats_prot_gather.yml"
    benchmark: "benchmarks/sgc_cami_low_k31_r1_extract_reads_{pfam}.tsv"
    resources:
        mem_mb = 200000,
        tmpdir = TMPDIR
    shell:'''
    python -Werror -m spacegraphcats.search.extract_reads {input.bgz} {input.idx} {input.cdbg_nbhds} -o {output}
    '''

rule grab_sgc_read_names:
    input: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.reads.gz"
    output: "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.read_names.txt"
    threads: 1
    benchmark: "benchmarks/sgc_cami_low_k31_r1_grab_read_names_{pfam}.tsv"
    resources:
        mem_mb = 10000,
        tmpdir = TMPDIR
    shell:'''
    zcat {input} | grep "@" > {output} 
    '''

rule download_pfam_id_to_name_map:
    output: "inputs/Pfam-A.clans.tsv.gz"
    threads: 1
    resources:
        mem_mb = 10000,
        tmpdir = TMPDIR
    shell:'''
    wget -O {output} ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
    '''

rule grab_gs_read_names:
    input:
        gs = "/home/tereiter/github/2021-cami-annot/outputs/gs_read_annotations/CAMI_low_gs_read_annotations_with_eggnog.tsv",
        pfam_map = "inputs/Pfam-A.clans.tsv.gz",
        pfam_sgc = "outputs/spacegraphcats/CAMI_low_k31_r1_multifasta_x_sequences/{pfam}.nbhds.read_names.txt"
    output: pfam_gs = "outputs/gs_read_sets/{pfam}_gs_read_names.tsv"
    threads: 1
    benchmark: "benchmarks/gs_grab_read_names_{pfam}.tsv"
    resources:
        mem_mb = 64000,
        tmpdir = TMPDIR
    script: "scripts/grab_gs_read_names.R"
    
