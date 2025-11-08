import os

configfile: 'config/config.yml'
ncores = config['ncores']
db_dir = config['db_dir']
gutsmash_dir = config['gutsmash_dir']

INPUT_FILE = config['input_file']
SAMPLES = []
OUTPUTS = []

os.system("chmod -R +x scripts")

samp2path = {}
with open(INPUT_FILE) as f:
    for line in f:
        path = line.strip()
        sample = os.path.basename(path.split(".fa")[0])
        output = os.path.dirname(path)
        SAMPLES.append(sample)
        OUTPUTS.append(output)
        samp2path[sample] = path

for sample in samp2path:
    dirname = os.path.dirname(samp2path[sample])+"/"+sample+"_annotations/logs"
    if not os.path.exists(dirname):
        os.makedirs(dirname)

if not os.path.exists(db_dir+"/amr_finder"):
        os.makedirs(db_dir+"/amr_finder")

if not os.path.exists(db_dir+"/antismash"):
        os.makedirs(db_dir+"/antismash")

rule targets:
    input:
        expand(["{output}/{sample}_annotations/done.txt", "{output}/{sample}_annotations/prokka/{sample}.ffn", "{output}/{sample}_annotations/eggnog.emapper.annotations", "{output}/{sample}_annotations/cazy_results.tsv", "{output}/{sample}_annotations/kegg_orthologs.tsv", "{output}/{sample}_annotations/kegg_modules.tsv", "{output}/{sample}_annotations/amrfinder_results.tsv", "{output}/{sample}_annotations/antismash/summary.tsv", "{output}/{sample}_annotations/gutsmash/summary.tsv"], zip, output=OUTPUTS, sample=SAMPLES)
               
rule prokka:
    input:
        lambda wildcards: samp2path[wildcards.sample]
    output:
        faa = "{output}/{sample}_annotations/prokka/{sample}.faa",
        ffn = "{output}/{sample}_annotations/prokka/{sample}.ffn",
        gff = "{output}/{sample}_annotations/prokka/{sample}.gff",
    params:
        outdir = "{output}/{sample}_annotations/prokka",
    conda:
        "config/envs/annotation.yml"
    resources:
        ncores = ncores
    shell:
        """
        if [[ {input} == *.gz ]]
        then
            gunzip -c {input} > {params.outdir}/{wildcards.sample}.fa
        else
            ln -fs {input} {params.outdir}/{wildcards.sample}.fa
        fi
        prokka --cpus {resources.ncores} {params.outdir}/{wildcards.sample}.fa --outdir {params.outdir} --prefix {wildcards.sample} --force --locustag {wildcards.sample} --rfam
        rm -rf {params.outdir}/*fna {params.outdir}/*err {params.outdir}/*fsa {params.outdir}/*gbk {params.outdir}/*log {params.outdir}/*sqn {params.outdir}/*tbl {params.outdir}/*tsv {params.outdir}/*txt {params.outdir}/{wildcards.sample}.fa
        """

rule eggnog:
    input:
        "{output}/{sample}_annotations/prokka/{sample}.faa"
    output:
        outfile = "{output}/{sample}_annotations/eggnog.emapper.annotations"
    params:
        out = "{output}/{sample}_annotations",
        db = db_dir+"/eggnog"
    conda:
        "config/envs/annotation.yml"
    resources:
        ncores = ncores
    shell:
        """
        rm -rf {params.out}/*tmp_dmd* {params.out}/done.txt
        emapper.py --cpu {resources.ncores} -i {input} -m diamond -o eggnog --output_dir {params.out} --temp_dir {params.out} --data_dir {params.db} --override
        rm -rf {params.out}/eggnog.emapper.hits {params.out}/eggnog.emapper.seed_orthologs
        """

rule cazy:
    input:
        "{output}/{sample}_annotations/prokka/{sample}.faa"
    output:
        parsed = "{output}/{sample}_annotations/cazy_results.tsv"
    params:
        out = "{output}/{sample}_annotations",
        db = db_dir+"/cazy"
    conda:
        "config/envs/annotation.yml"
    resources:
        ncores = ncores
    shell:
        """
        run_dbcan.py --dia_cpu {resources.ncores} --hmm_cpu {resources.ncores} --tf_cpu {resources.ncores} --hotpep_cpu {resources.ncores} {input} protein --db_dir {params.db} --out_dir {params.out}
        python scripts/dbcan_simplify.py {params.out}/overview.txt > {output.parsed}
        rm -rf {params.out}/diamond.out {params.out}/hmmer.out {params.out}/Hotpep.out {params.out}/uniInput {params.out}/overview.txt
        """

rule kofam:
    input:
        "{output}/{sample}_annotations/prokka/{sample}.faa"
    output:
        kofam = "{output}/{sample}_annotations/kegg_orthologs.tsv",
        modules = "{output}/{sample}_annotations/kegg_modules.tsv"
    params:
        out = "{output}/{sample}_annotations",
        db = db_dir+"/kofam/kofam_db.hmm",
        graph = db_dir+"/kofam/pathways/graphs.pkl",
        names = db_dir+"/kofam/pathways/all_pathways_names.txt",
        classes = db_dir+"/kofam/pathways/all_pathways_class.txt"
    conda:
        "config/envs/annotation.yml"
    resources:
        ncores = ncores
    shell:
        """
        python scripts/kofamscan.py -t {resources.ncores} -q {input} -o {params.out} -d {params.db} > /dev/null
        python scripts/give_pathways.py -i {output.kofam} -d {params.out} -g {params.graph} -n {params.names} -c {params.classes} -o {wildcards.sample}
        rm -rf {params.out}/kofam_raw.tsv
        """

rule amrfinder_setup:
    input:
        db_dir+"/amr_finder/"
    output:
        directory(db_dir+"/amr_finder/2023-04-17.1")
    conda:
        "config/envs/annotation.yml"
    shell:
        """
        wget ftp://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/3.11/2023-04-17.1/* -P {output} 
        amrfinder_index {output}
        """

rule amrfinder_run:
    input:
        db = db_dir+"/amr_finder/2023-04-17.1",
        gff = "{output}/{sample}_annotations/prokka/{sample}.gff",
        faa = "{output}/{sample}_annotations/prokka/{sample}.faa",
        fa = lambda wildcards: samp2path[wildcards.sample]
    output:
        amr = "{output}/{sample}_annotations/amrfinder_results.tsv"
    params:
        outdir = "{output}/{sample}_annotations",
        outfa = "{output}/{sample}_annotations/amrfinder.fa"
    conda:
        "config/envs/annotation.yml"
    resources:
        ncores = ncores
    shell:
        """
        if [[ {input.fa} == *.gz ]]
        then
            gunzip -c {input.fa} > {params.outfa}
        else
            ln -fs {input.fa} {params.outfa}
        fi
        grep -w CDS {input.gff} | perl -pe '/^##FASTA/ && exit; s/(\W)Name=/$1OldName=/i; s/ID=([^;]+)/ID=$1;Name=$1/' > {params.outdir}/amrfinder.gff
        amrfinder --threads {resources.ncores} -p {input.faa} -g {params.outdir}/amrfinder.gff -n {params.outfa} -o {output.amr} -d {input.db} --plus
        rm -f {params.outdir}/*.tmp.* {params.outdir}/amrfinder.gff {params.outfa}
        """

rule gutsmash:
    input:
        gff = "{output}/{sample}_annotations/prokka/{sample}.gff",
        fa = lambda wildcards: samp2path[wildcards.sample]
    output:
        tsv = "{output}/{sample}_annotations/gutsmash/summary.tsv",
        main = directory("{output}/{sample}_annotations/gutsmash")
    params:
        gut_exec = gutsmash_dir+"/run_gutsmash.py",
        outfa = "{output}/{sample}_annotations/gutsmash.fa"
    conda:
        "config/envs/gutsmash.yml"
    resources:
        ncores = ncores,
    shell:
        """
        if [[ {input.fa} == *.gz ]]
        then
            gunzip -c {input.fa} > {params.outfa}
        else
            ln -fs {input.fa} {params.outfa}
        fi
        {params.gut_exec} -c {resources.ncores} --cb-knownclusters --genefinding-gff3 {input.gff} --enable-genefunctions {params.outfa} --output-dir {output.main}
        rm -rf {output.main}/gutsmash.gbk {output.main}/gutsmash.zip {output.main}/gutsmash.fa {output.main}/html {output.main}/knownclusterblast* {output.main}/gutsmash.json {output.main}/css {output.main}/images {output.main}/index.html {output.main}/js {output.main}/regions.js {output.main}/svg
        python scripts/gutsmash2tsv.py {output.main} {wildcards.sample} > {output.tsv}
	"""

rule antismash_setup:
    input:
        db_dir+"/antismash"
    output:
        directory(db_dir+"/antismash/clusterblast"),
        directory(db_dir+"/antismash/clustercompare"),
        directory(db_dir+"/antismash/pfam"),
        directory(db_dir+"/antismash/resfam"),
        directory(db_dir+"/antismash/tigrfam")
    conda:
        "config/envs/antismash.yml"
    shell:
        "download-antismash-databases --database-dir {input}"

rule antismash_run:
    input:
        gff = "{output}/{sample}_annotations/prokka/{sample}.gff",
        fa = lambda wildcards: samp2path[wildcards.sample],
        db = db_dir+"/antismash",
        cblast = db_dir+"/antismash/clusterblast",
        ccomp = db_dir+"/antismash/clustercompare",
        pfam = db_dir+"/antismash/pfam",
        resfam = db_dir+"/antismash/resfam",
        tigrfam = db_dir+"/antismash/tigrfam",
    output:
        tsv = "{output}/{sample}_annotations/antismash/summary.tsv",
        json = "{output}/{sample}_annotations/antismash/antismash.json",
        main = directory("{output}/{sample}_annotations/antismash")
    params:
        outfa = "{output}/{sample}_annotations/antismash.fa",
        parent = "{output}/{sample}_annotations"
    conda:
        "config/envs/antismash.yml"
    resources:
        ncores = ncores,
        tmpdir = lambda wildcards: wildcards.output+"/"+wildcards.sample+"_annotations/antismash_tmp"
    shell:
        """
        rm -f {params.parent}/done.txt
        if [[ {input.fa} == *.gz ]]
        then
            gunzip -c {input.fa} > {params.outfa}
        else
            ln -fs {input.fa} {params.outfa}
        fi
        antismash -v -c {resources.ncores} --skip-zip-file --allow-long-headers --databases {input.db} --cc-mibig --cb-general --cb-knownclusters --cb-subclusters --asf --pfam2go --smcog-trees --genefinding-gff3 {input.gff} --output-dir {output.main} --output-basename antismash {params.outfa}
        rm -rf {output.main}/clusterblast* {output.main}/smcogs {output.main}/subcluster* {output.main}_tmp {output.main}/antismash.gbk {output.main}/css {output.main}/images {output.main}/index.html {output.main}/html {output.main}/js {output.main}/regions.js {output.main}/svg {params.outfa}
        python scripts/antismash2tsv.py {output.main} {wildcards.sample} > {output.tsv}
        """

rule clean_up:
    input:
        prokka = "{output}/{sample}_annotations/prokka/{sample}.faa",
        amr = "{output}/{sample}_annotations/amrfinder_results.tsv",
        kegg = "{output}/{sample}_annotations/kegg_orthologs.tsv",
        cazy = "{output}/{sample}_annotations/cazy_results.tsv",
        antis = "{output}/{sample}_annotations/antismash/summary.tsv",
        guts = "{output}/{sample}_annotations/gutsmash/summary.tsv",
        eggnog = "{output}/{sample}_annotations/eggnog.emapper.annotations"
    output:
        "{output}/{sample}_annotations/done.txt"
    params:
        anti_tmp = "{output}/{sample}_annotations/antismash",
        eggnog_tmp = "{output}/{sample}_annotations"
    resources:
        ncores = 1
    shell:
        """
        if [[ -d {params.anti_tmp}"_tmp" ]]
        then
            rm -rf {input.antis} {params.anti_tmp}"_tmp"
        fi

        if [[ $(find {params.eggnog_tmp} -type d -name "emappertmp*" | wc -l)  -gt "0" ]]
        then
            rm -rf {input.eggnog}
        fi

        touch {output}
        """
