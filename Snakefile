############################################################
## Imports
############################################################
import os

############################################################
## Messages
############################################################

############################################################
## Config file and shorthands
############################################################

configfile: "config.json"
workdir: config["WORKFLOW"]["WORKDIR"]

databases_dir = "databases"
tmp_dir = "tmp"
log_dir = "logs"

############################################################
## All rule
############################################################

SAMPLES, = glob_wildcards("samples/{sample_name}.fastq")
DEDUP_STRING = ""
if config["MISC"]["DEDUPLICATE"] == "True":
    DEDUP_STRING = config["EXTENSION"]["DEDUPLICATED"]
else:
    DEDUP_STRING = config["EXTENSION"]["NOT_DEDUPLICATED"]

EXTS = [
    DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"],
    DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"] + config["EXTENSION"]["RESISTOME_RICHNESS"],
    DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"] + config["EXTENSION"]["RESISTOME_DIVERSITY"],
    DEDUP_STRING + "_" + config["MISC"]["MOBILOME_STRATEGY"] + config["EXTENSION"]["MOBILOME"]]

rule all:
    input:
        # Data
        expand("{sample_name}.fastq{ext}", sample_name=SAMPLES, ext=EXTS)

############################################################
## Pipeline Steps
############################################################

############################################################
# Utils

rule read_lengths:
    input:
        reads = "samples/{sample_name}.fastq"

    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8"

    output:
        read_lenghts_json = "{sample_name}.fastq" + config["EXTENSION"]["READS_LENGTH"]

    shell:
        """
        python3 {params.read_lengths_script} {input.reads} > {output.read_lenghts_json}
        """

rule read_lengths_from_workdir:
    input:
        reads = "{sample_name}.fastq"

    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8"

    output:
        read_lenghts_json = "{sample_name}.fastq" + config["EXTENSION"]["READS_LENGTH"]

    shell:
        """
        python3 {params.read_lengths_script} {input.reads} > {output.read_lenghts_json}
        """

ruleorder: read_lengths > read_lengths_from_workdir

############################################################
# Deduplication

rule deduplicate_create_clusters:
    input:
        reads = "samples/{sample_name}.fastq",
        reads_lengths = "{sample_name}.fastq" + config["EXTENSION"]["READS_LENGTH"]

    params:
        num_of_clusters = config["MISC"]["DEDUP_CLUSTERS"],
        clustering_script = workflow.basedir + "/" + config["SCRIPTS"]["CLUSTER_READS"],
        tmp_dir_clusters = tmp_dir + "/tmp_{sample_name}"

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8"

    output:
        clusters = [tmp_dir + "/tmp_{{sample_name}}/{{sample_name}}_{cl_id}.fasta.gz".format(cl_id=cl_id)
                        for cl_id in range(config["MISC"]["DEDUP_CLUSTERS"])]

    shell:
        """
        mkdir -p {params.tmp_dir_clusters}
        python3 {params.clustering_script} -r {input.reads} -o {params.tmp_dir_clusters} \
            -n {params.num_of_clusters} -l {input.reads_lengths}
        """

rule deduplicate_blat:
    input:
        reads_cluster = tmp_dir + "/tmp_{sample_name}/{sample_name}_{cl_id}.fasta.gz"

    params:
        out_pls_dir = tmp_dir + "/tmp_{sample_name}/pls_files"

    threads: config["MISC"]["BLAT_THREADS"]

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "blat/20140318"

    output:
        out_pls = tmp_dir + "/tmp_{sample_name}/pls_files/{cl_id}.pls"

    shell:
        """
        mkdir -p {params.out_pls_dir}
        blat {input.reads_cluster} {input.reads_cluster} {output.out_pls}
        """

rule deduplicate_find_duplicates:
    input:
        pls_file = tmp_dir + "/tmp_{sample_name}/pls_files/{cl_id}.pls"

    params:
        find_duplicates_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_DUPLICATES"],
        similarity_threshold = config["MISC"]["DEDUPLICATION_SIMILARITY_THRESHOLD"]

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8"

    output:
        out_duplicates_csv = tmp_dir + "/tmp_{sample_name}/pls_files/{cl_id}_duplicates.csv"

    shell:
        """
        python3 {params.find_duplicates_script} -p {input.pls_file} -s {params.similarity_threshold} > {output.out_duplicates_csv}
        """

rule deduplicate_merge_duplicates_sets:
    input:
        expand(tmp_dir + "/tmp_{{sample_name}}/pls_files/{cl_id}_duplicates.csv", cl_id=range(config["MISC"]["DEDUP_CLUSTERS"]))

    output:
        out_duplicates_csv = tmp_dir + "/tmp_{sample_name}/merged_duplicates.csv"

    run:
        import csv
        with open(output[0], "w") as out_csv_handle:
            csv_writer = csv.writer(out_csv_handle)
            for cluster_duplicate_file in input:
                with open(cluster_duplicate_file) as csv_handle:
                    csv_reader = csv.reader(csv_handle)
                    for row in csv_reader:
                        csv_writer.writerow(row)

rule deduplicate:
    input:
        reads = "samples/{sample_name}.fastq",
        duplicates_csv = tmp_dir + "/tmp_{sample_name}/merged_duplicates.csv"

    params:
        deduplicate_script = workflow.basedir + "/" + config["SCRIPTS"]["DEDUPLICATE"]

    conda:
        "envs/deduplication.yaml"
    envmodules:
        "python/3.8"

    output:
        deduplicated_reads = "{sample_name}.fastq" + config["EXTENSION"]["DEDUPLICATED"]

    shell:
        """
        python3 {params.deduplicate_script} -r {input.reads} -d {input.duplicates_csv} > {output.deduplicated_reads}
        """

rule not_deduplicated_reads:
    input:
        reads = "samples/{sample_name}.fastq"

    output:
        out_reads = "{sample_name}.fastq" + config["EXTENSION"]["NOT_DEDUPLICATED"]

    shell:
        """
        ln -s {input.reads} {output.out_reads}
        """

############################################################s
# Alignment

rule align_to_megares:
    input:
        reads = "{sample_name}.fastq" + DEDUP_STRING,
        megares_v2_seqs = databases_dir + "/" + "megares_full_database_v2.00.fasta"

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        megares_out_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.megares_v2_seqs} {input.reads} -o {output.megares_out_sam}
        """

rule align_to_mges:
    input:
        reads = "{sample_name}.fastq" + DEDUP_STRING,
        mges_database = databases_dir + "/" + "mges_combined.fa"
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        mges_out_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.mges_database} {input.reads} -o {output.mges_out_sam}
        """

rule align_to_kegg:
    input:
        reads = "{sample_name}.fastq" + DEDUP_STRING,
        kegg_database = databases_dir + "/" + "kegg_genes.fa"

    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]

    conda:
        "envs/alignment.yaml"
    envmodules:
        "python/3.8",
        "minimap/2.21"

    threads: config["MINIMAP2"]["THREADS"]

    output:
        kegg_out_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_KEGG"]

    shell:
        """
        minimap2 -t {threads} {params.minimap_flags} {input.kegg_database} {input.reads} -o {output.kegg_out_sam}
        """

rule pass_config_file:
    output:
        out_config_file = "config.ini"

    run:
        import configparser
        with open(output.out_config_file,'w') as configfile_out:
            config_to_pass = dict(config)
            config_to_pass["DATABASE"]["MEGARES"] = databases_dir + "/" + "megares_full_database_v2.00.fasta"
            config_to_pass["DATABASE"]["MEGARES_ONTOLOGY"] = databases_dir + "/" + "megares_full_annotations_v2.00.csv"
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)

rule resistome_and_mobilome:
    input:
        megares_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_lenght = "{sample_name}.fastq" + config["EXTENSION"]["READS_LENGTH"],
        dedup_reads_lenght = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        config_file = "config.ini"

    params:
        resistome_mobilome_script = workflow.basedir + "/" + config["SCRIPTS"]["GEN_RESISTOME_AND_MOBILOME"],
        output_prefix = "{sample_name}.fastq" + DEDUP_STRING

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        resistome_richness = "{sample_name}.fastq" + DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                    + config["EXTENSION"]["RESISTOME_RICHNESS"],
        resistome_diversity = "{sample_name}.fastq" + DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                     + config["EXTENSION"]["RESISTOME_DIVERSITY"],
        mobilome = "{sample_name}.fastq" + DEDUP_STRING + "_" + config["MISC"]["MOBILOME_STRATEGY"] + config["EXTENSION"]["MOBILOME"]

    shell:
        """
        python3 {params.resistome_mobilome_script} \
            -r {wildcards.sample_name}.fastq \
            -a {input.megares_sam} \
            -m {input.mges_sam} \
            -c {input.config_file} \
            -o {params.output_prefix}
        """

rule find_colocalizations:
    input:
        reads = "{sample_name}.fastq" + DEDUP_STRING,
        megares_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        kegg_sam = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["A_TO_KEGG"],
        reads_lenght = "{sample_name}.fastq" + config["EXTENSION"]["READS_LENGTH"],
        dedup_reads_lenght = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        config_file = "config.ini"

    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_COLOCALIZATIONS"],
        output_directory = os.getcwd()

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        colocalizations = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"]

    shell:
        """
        python3 {params.find_colocalizations_script} \
            -r {input.reads} \
            --arg {input.megares_sam} \
            --mge {input.mges_sam} \
            -k {input.kegg_sam} \
            -c {input.config_file} \
            -o {params.output_directory} \
            > {output.colocalizations}
        """

rule colocalization_richness:
    input:
        colocalizations = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"

    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["COLOCALIZATIONS_RICHNESS"]

    conda:
        "envs/pipeline.yaml"
    envmodules:
        "python/3.8"

    output:
        colocalizations_richness = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"]

    shell:
        """
        python3 {params.find_colocalizations_script} \
            -i {input.colocalizations} \
            -c {input.config_file} \
            > {output.colocalizations_richness}
        """


############################################################
## Plots
############################################################

rule read_lengths_plot:
    input:
        reads = "{sample_name}.fastq" + DEDUP_STRING

    params:
        num_of_bins = 100,
        std_deviations = 4,
        read_lengths_script = workflow.basedir + "/" + "src/plot_read_lengths.py"

    output:
        out_plot_name = "{sample_name}" + "_deduplicated_read_lengts_hist.pdf"

    conda:
        "envs/plots.yaml"
    envmodules:
        "python/3.8"

    shell:
        """
        python3 {params.read_lengths_script} -i {input.reads} -s {params.std_deviations} -b {params.num_of_bins} \
            -o {output.out_plot_name} --title {wildcards.sample_name}_deduplicated
        """

rule violin_plots_notebook:
    input:
        megares_db = databases_dir + "/megares_full_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_full_annotations_v2.00.csv",
        data = expand("{sample_name}.fastq{ext}",sample_name=SAMPLES,ext=EXTS)

    params:
        samples_list = SAMPLES,
        config_dict  = config

    output:
        out_plot_name = "violin_plot_all_samples.pdf"

    conda:
        "envs/plots.yaml"
    envmodules:
        "python/3.8"

    notebook:
        "src/plots_notebooks/violin_notebook.py.ipynb"


rule heatmap_notebook:
    input:
        megares_db = databases_dir + "/megares_full_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_full_annotations_v2.00.csv",
        data = expand("{sample_name}.fastq{ext}",sample_name=SAMPLES,ext=EXTS)

    params:
        samples_list = SAMPLES,
        config_dict  = config

    output:
        out_plot_name = "heatmap_all_samples.pdf"

    conda:
        "envs/plots.yaml"
    envmodules:
        "python/3.8"

    notebook:
        "src/plots_notebooks/heatmap_notebook.py.ipynb"


rule colocalization_visualizations_notebook:
    input:
        megares_db = databases_dir + "/megares_full_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_full_annotations_v2.00.csv",
        mges_db = databases_dir + "/mges_combined.fa",
        kegg_db = databases_dir + "/kegg_genes.fa",
        dedup_reads_lenght = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        colocalizations = "{sample_name}.fastq" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"]

    params:
        config_dict  = config

    output:
        out_plot_name = "{sample_name}_colocalizations_plot.pdf"

    conda:
        "envs/plots.yaml"
    envmodules:
        "python/3.8"

    notebook:
        "src/plots_notebooks/colocalizations_notebook.py.ipynb"

############################################################
## Databases
############################################################

rule get_megares_v2:
    output:
        megares_v2_seqs = os.path.join(databases_dir,"megares_full_database_v2.00.fasta"),
        megares_v2_ontology = os.path.join(databases_dir,"megares_full_annotations_v2.00.csv")

    shell:
        """
        mkdir -p {databases_dir}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_database_v2.00.fasta -O {output.megares_v2_seqs}
        wget http://megares.meglab.org/download/megares_v2.00/megares_full_annotations_v2.00.csv -O {output.megares_v2_ontology}
        """

rule get_MGEs_DBs:
    output:
        mges_combined_db = os.path.join(databases_dir,"mges_combined.fa")

    params:
        mges_path = config["DATABASE"]["MGES"]

    shell:
        """
        mkdir -p {databases_dir}
        cp {params.mges_path} {output.mges_combined_db}
        """

rule get_KEGG_Prokaryotes_DBs:
    output:
        kegg_prokaryotes_db = os.path.join(databases_dir,"kegg_genes.fa")

    params:
        kegg_path = config["DATABASE"]["KEGG"]

    shell:
        """
        mkdir -p {databases_dir}
        cp {params.kegg_path} {output.kegg_prokaryotes_db}
        """

############################################################
## Cleans
############################################################

rule clean:
    shell:
        """
        rm -f *.csv *.sam *.json *.pdf *_deduplicated.fastq
        rm -rf {databases_dir} {tmp_dir}
        rm -f config.ini
        """
