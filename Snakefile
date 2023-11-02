# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

################################################################################
#                                Imports                                       #
################################################################################

import os, ast, glob

################################################################################
#                       Config File & Shorthands                               #
################################################################################

configfile: "config/config.json"
workdir: config["WORKFLOW"]["WORKDIR"]

databases_dir = config["WORKFLOW"]["DATABASES_DIR"]
samples_dir = "samples"
tmp_dir = "tmp"
log_dir = "logs"

################################################################################
#                                 All Rule                                     #
################################################################################

#DEDUP_STRING = ""
if config["WORKFLOW"]["LONG_READS"].upper() == "TRUE":
    SAMPLES, = glob_wildcards(samples_dir + "/{sample_name}.fastq")
    if config["WORKFLOW"]["DEDUPLICATE"].upper() == "TRUE":
        DEDUP_STRING = config["EXTENSION"]["DEDUPLICATED"]
    else:
        DEDUP_STRING = config["EXTENSION"]["NOT_DEDUPLICATED"]
else:
    SAMPLES, = glob_wildcards(samples_dir + "/{sample_name}.fastq") # CHANGE THIS TO FASTA, ALSO CONSIDER READ PAIRS
    DEDUP_STRING = config["EXTENSION"]["ASSEMBLY"]

EXTS = [
    DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
    DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"],
    DEDUP_STRING + config['EXTENSION']['OVERLAP'],
    DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"] + config["EXTENSION"]["RESISTOME_RICHNESS"],
    DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"] + config["EXTENSION"]["RESISTOME_DIVERSITY"],
    DEDUP_STRING + "_" + config["MISC"]["MOBILOME_STRATEGY"] + config["EXTENSION"]["MOBILOME"]]

rule all:
    input:
        expand("{sample_name}{ext}", sample_name=SAMPLES, ext=EXTS)

################################################################################
#                             Long Read Steps                                  #
################################################################################

if config["WORKFLOW"]["LONG_READS"].upper() == "TRUE":

########################### Deduplication ######################################
    
    if config["WORKFLOW"]["DEDUPLICATE"].upper() == "TRUE":
        rule bin_reads_by_length:
            input:
                samples_dir + "/{sample_name}.fastq"
            output:
                touch(tmp_dir + "/{sample_name}.bin.reads.done")
            params:
                outdir = tmp_dir + "/{sample_name}_read_bins",
                bin_script = workflow.basedir + "/" + config["SCRIPTS"]["BIN_READS"]
            conda:
                "workflow/envs/deduplication.yaml"
            shell:
                "mkdir -p {params.outdir}; "
                "{params.bin_script} "
                "--infile {input} "
                "--outdir {params.outdir}"

        rule cluster_reads:
            input:
                tmp_dir + "/{sample_name}.bin.reads.done"
            output:
                touch(tmp_dir + "/{sample_name}.cluster.reads.done")
            params:
                indir = tmp_dir + "/{sample_name}_read_bins",
                outdir = tmp_dir + "/{sample_name}_read_clusters",
                cluster_script = workflow.basedir + "/" + config["SCRIPTS"]["CLUSTER_READS"]
            conda:
                "workflow/envs/deduplication.yaml"
            shell:
                "mkdir -p {params.outdir}; "
                "{params.cluster_script} "
                "--indir {params.indir} "
                "--outdir {params.outdir}; "
                "rm -rf {params.indir}; "
                "rm {input}"

        rule blat_clustered_reads:
            input:
                tmp_dir + "/{sample_name}.cluster.reads.done"
            output:
                touch(tmp_dir + "/{sample_name}.blat.done")
            params:
                rc = tmp_dir + "/{sample_name}_read_clusters/",
                o = tmp_dir + "/{sample_name}_psl_files/",
                blat_script = workflow.basedir + "/" + config["SCRIPTS"]["RUN_BLAT"]
            threads:
                32
            conda:
                "workflow/envs/deduplication.yaml"
            shell:
                "mkdir -p {params.o}; "
                "{params.blat_script} "
                "--outdir {params.o} "
                "--threads {threads} "
                "--read_clusters {params.rc}; "
                "rm -rf {params.rc}; "
                "rm {input}"

        rule find_duplicates:
            input:
                tmp_dir + "/{sample_name}.blat.done"
            output:
                touch(tmp_dir + "/{sample_name}.find.duplcates.done")
            params:
                similarity_threshold = config["MISC"]["DEDUPLICATION_SIMILARITY_THRESHOLD"],
                pls_dir = tmp_dir + "/{sample_name}_psl_files/",
                outdir = tmp_dir + "/{sample_name}_duplicate_txts/",
                find_dups_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_DUPLICATES"]
            conda:
                "workflow/envs/deduplication.yaml"
            shell:
                "mkdir -p {params.outdir}; "
                "{params.find_dups_script} "
                "{params.outdir} "
                "{params.pls_dir} "
                "{params.similarity_threshold}; "
                "rm -rf {params.pls_dir}; "
                "rm {input}"

        rule merge_duplicates_lists:
            input:
                tmp_dir + "/{sample_name}.find.duplcates.done"
            output:
                tmp_dir + "/{sample_name}.duplicates.txt"
            params:
                indir = tmp_dir + "/{sample_name}_duplicate_txts/"
            shell:
                "cat {params.indir}/* > {output}; "
                "rm -rf {params.indir}"

        # change output so that it is NOT gzipped *********
        rule deduplicate:
            input:
                reads = samples_dir + "/{sample_name}.fastq",
                duplicates_list = tmp_dir + "/{sample_name}.duplicates.txt"
            output:
                reads =  "{sample_name}" + config["EXTENSION"]["DEDUPLICATED"],
                dupes = "{sample_name}.dup.reads.fastq"
            params:
                dedup_script = workflow.basedir + "/" + config["SCRIPTS"]["DEDUPLICATE"]
            conda:
                "workflow/envs/deduplication.yaml"
            shell:
                "{params.dedup_script} "
                "--reads {input.reads} "
                "--duplicates {input.duplicates_list} "
                "--out_reads {output.reads} "
                "--out_dupes {output.dupes}"
        
        READS = "{sample_name}" + config["EXTENSION"]["DEDUPLICATED"]
    else:
        rule not_deduplicated_reads:
            input:
                samples_dir + "/{sample_name}.fastq"
            output:
                "{sample_name}" + config["EXTENSION"]["NOT_DEDUPLICATED"]
            shell:
                "ln -s {input} {output}"
        
        READS = "{sample_name}" + config["EXTENSION"]["NOT_DEDUPLICATED"]

################################################################################
#                           Short Reads Steps                                  #
################################################################################

elif config["WORKFLOW"]["LONG_READS"].upper() == "FALSE":

############################## Assembly ########################################

    rule meta_spades_assembly:
        input:
            f_reads = samples_dir + "/{sample_name}_R1.fastq",
            r_reads = samples_dir + "/{sample_name}_R2.fastq"
        params:
            kmers = confgig["SPADES"]["KMERS"],
            memory = config["SPADES"]["MEMORY"],
            phred = config["SPADES"]["PHRED"]
        conda:
            "workflow/envs/assembly.yaml"
        threads:
            config["SPADES"]["THREADS"]
        output:
            temp(tmp_dir + "/spades_files/{sample_name}/")
        shell:
            "spades.py --meta "
            "--phred-offset {params.phred} "
            "-t {threads} "
            "-m {params.memory} "
            "-k {params.kmers} "
            "-1 {input.f_reads} "
            "-2 {input.r_reads} "
            "-o {output}"

    READS = tmp_dir + "/spades_files/{sample_name}/scaffolds.fasta"

############################
# Read Lengths #
############################## Read Lengths ####################################

rule read_lengths:
    input:
        READS
    output:
        # "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"]
         "{sample_name}" + config["EXTENSION"]["READS_LENGTH"]
    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]
    conda:
        "workflow/envs/deduplication.yaml"
    shell:
        "python {params.read_lengths_script} {input} > {output}"

################################################################################
#                                Alignment                                     #
################################################################################

rule align_to_megares:
    input:
        reads = READS,
        megares_seqs = ancient(databases_dir + "/" + "megares_modified_database_v2.00.fasta")
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"]
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]
    conda:
        "workflow/envs/alignment.yaml"
    threads:
        config["MINIMAP2"]["THREADS"]
    shell:
        "minimap2 -Y "
        "-t {threads} "
        "{params.minimap_flags} "
        "{input.megares_seqs} "
        "{input.reads} "
        "-o {output}"

rule align_to_mges:
    input:
        reads = READS,
        mges_database = ancient(databases_dir + "/" + "mges_combined.fasta")
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"]
    params:
        minimap_flags = config["MINIMAP2"]["ALIGNER_PB_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_ONT_OPTION"] + " "
                        + config["MINIMAP2"]["ALIGNER_HIFI_OPTION"]
    conda:
        "workflow/envs/alignment.yaml"
    threads:
        config["MINIMAP2"]["THREADS"]
    shell:
        "minimap2 -Y "
        "-t {threads} "
        "{params.minimap_flags} "
        "{input.mges_database} "
        "{input.reads} "
        "-o {output}"

rule pass_config_file:
    output:
        out_config_file = "config.ini"
    run:
        import configparser
        with open(output.out_config_file,'w') as configfile_out:
            config_to_pass = dict(config)
            config_to_pass["DATABASE"] = dict()
            config_to_pass["DATABASE"]["MEGARES"] = databases_dir + "/" + "megares_modified_database_v2.00.fasta"
            config_to_pass["DATABASE"]["MEGARES_ONTOLOGY"] = databases_dir + "/" + "megares_modified_annotations_v2.00.csv"
            config_to_pass["DATABASE"]["MGES"] = databases_dir + "/" + "mges_combined.fasta"
            config_parser = configparser.ConfigParser()
            config_parser.read_dict(config_to_pass)
            config_parser.write(configfile_out)

rule overlap:
    input:
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_length = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"], # I have no idea where these come from
        config_file = "config.ini"
    params:
        overlap_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_OVERLAP"],
        output_prefix = "{sample_name}" + DEDUP_STRING
    conda:
        "workflow/envs/pipeline.yaml"
    output:
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP']
    shell:
        "python {params.overlap_script} "
        "-r {wildcards.sample_name}.fastq "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-c {input.config_file} "
        "-o {params.output_prefix}"

rule merge_overlap_info:
    input:
        expand("{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'], sample_name = SAMPLES)
    output: 
        merged_info = 'merged_overlaped_mges_info.csv'
    run:
        import csv
        merged_overlaped_mges_info = set()
        for overlap_info_filename in input:
            with open(overlap_info_filename) as overlap_info_file:
                overlap_reader = csv.reader(overlap_info_file, delimiter=',')
                for row in overlap_reader:
                    merged_overlaped_mges_info.add(row[0])
        with open(output[0], 'w') as merged:
            merged_writer = csv.writer(merged, delimiter=',')
            for mge in merged_overlaped_mges_info:
                merged_writer.writerow([mge])

rule resistome_and_mobilome:
    input:
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_length = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'],
        config_file = "config.ini"
    output:
        resistome_richness = "{sample_name}" + DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                    + config["EXTENSION"]["RESISTOME_RICHNESS"],
        resistome_diversity = "{sample_name}" + DEDUP_STRING + "_" + config["MISC"]["RESISTOME_STRATEGY"]
                                                     + config["EXTENSION"]["RESISTOME_DIVERSITY"],
        mobilome = "{sample_name}" + DEDUP_STRING + "_" + config["MISC"]["MOBILOME_STRATEGY"] + config["EXTENSION"]["MOBILOME"]
    params:
        resistome_mobilome_script = workflow.basedir + "/" + config["SCRIPTS"]["GEN_RESISTOME_AND_MOBILOME"],
        output_prefix = "{sample_name}" + DEDUP_STRING
    conda:
        "workflow/envs/pipeline.yaml"
    shell:
        "python {params.resistome_mobilome_script} "
        "-r {wildcards.sample_name}.fastq "
        "-a {input.megares_sam} "
        "-m {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_prefix}"

rule find_colocalizations:
    input:
        reads = READS,
        megares_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MEGARES"],
        mges_sam = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["A_TO_MGES"],
        reads_length = "{sample_name}" + config["EXTENSION"]["READS_LENGTH"],
        #dedup_reads_lenght = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        overlap = "{sample_name}" + DEDUP_STRING + config['EXTENSION']['OVERLAP'],
        config_file = "config.ini"
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"]
    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["FIND_COLOCALIZATIONS"],
        output_directory = os.getcwd()
    conda:
        "workflow/envs/pipeline.yaml"
    shell:
        "python {params.find_colocalizations_script} "
        "-r {input.reads} "
        "--arg {input.megares_sam} "
        "--mge {input.mges_sam} "
        "-s {input.overlap} "
        "-c {input.config_file} "
        "-o {params.output_directory} "
        "> {output}"

rule colocalization_richness:
    input:
        colocalizations = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"
    output:
        "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS_RICHNESS"]
    params:
        find_colocalizations_script = workflow.basedir + "/" + config["SCRIPTS"]["COLOCALIZATIONS_RICHNESS"]
    conda:
        "workflow/envs/pipeline.yaml"
    shell:
        "python {params.find_colocalizations_script} "
        "-i {input.colocalizations} "
        "-c {input.config_file} "
        "> {output}"


############################################################
## Plots
############################################################

rule read_lengths_plot:
    input:
        "{sample_name}" + DEDUP_STRING
    output:
        "{sample_name}" + "_deduplicated_read_lengts_hist.pdf"
    params:
        num_of_bins = 100,
        std_deviations = 4,
        read_lengths_script = workflow.basedir + "/" + "scripts/plot_read_lengths.py"
    conda:
        "workflow/envs/plots.yaml"
    shell:
        "python {params.read_lengths_script} "
        "-i {input} "
        "-s {params.std_deviations} "
        "-b {params.num_of_bins} "
        "-o {output} "
        "--title {wildcards.sample_name}_deduplicated"

rule violin_plots_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta  ",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        config_file = "config.ini",
        data = expand("{sample_name}.fastq{ext}",sample_name=SAMPLES,ext=EXTS)
    output:
        "violin_plot_all_samples.pdf"
    params:
        samples_list = SAMPLES
    conda:
        "workflow/envs/plots.yaml"
    notebook:
        "workflow/notebooks/violin_notebook.py.ipynb"


rule heatmap_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        config_file = "config.ini",
        data = expand("{sample_name}{ext}",sample_name=SAMPLES,ext=EXTS)
    output:
        "heatmap_all_samples.pdf"
    params:
        samples_list = SAMPLES
    conda:
        "workflow/envs/plots.yaml"
    notebook:
        "workflow/notebooks/heatmap_notebook.py.ipynb"

# Going to have to change code here as well
rule colocalization_visualizations_notebook:
    input:
        megares_db = databases_dir + "/megares_modified_database_v2.00.fasta",
        megares_annotation = databases_dir + "/megares_modified_annotations_v2.00.csv",
        mges_db = databases_dir + "/mges_combined.fasta",
        dedup_reads_lenght = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["READS_LENGTH"],
        colocalizations = "{sample_name}" + DEDUP_STRING + config["EXTENSION"]["COLOCALIZATIONS"],
        config_file = "config.ini"
    output:
        "{sample_name}_colocalizations_plot.pdf"
    conda:
        "workflow/envs/plots.yaml"
    notebook:
        "workflow/notebooks/colocalizations_notebook.py.ipynb"

############################################################
## Databases
############################################################

rule get_megares:
    output:
        megares_seqs = os.path.join(databases_dir,"megares_modified_database_v2.00.fasta"),
        megares_ontology = os.path.join(databases_dir,"megares_modified_annotations_v2.00.csv")
    params:
        megares_seqs = "https://www.meglab.org/downloads/megares_v3.00/megares_database_v3.00.fasta",
        megares_ann = "https://www.meglab.org/downloads/megares_v3.00/megares_annotations_v3.00.csv"
    conda:
        "workflow/envs/download_databases.yaml"
    shell:
        "mkdir -p {databases_dir}; "
        "wget {params.megares_seqs} -O {output.megares_seqs}; "
        "wget {params.megares_ann} -O {output.megares_ontology}"

# There has to be a better way to get the plasmid finder db...

rule get_plasmid_finder_db:
    params:
        git_repo = "https://bitbucket.org/genomicepidemiology/plasmidfinder_db.git",
        commit = config["MISC"]["PLASMID_FINDER_COMMIT"]
    output:
        temp(databases_dir + "/plasmid_finder_db.fasta")
    conda:
        "workflow/envs/download_databases.yaml"
    shell:
        "mkdir -p {tmp_dir}; "
        "mkdir -p {databases_dir}; "
        "cd {tmp_dir}; "
        "git clone {params.git_repo}; "
        "cd plasmidfinder_db; "
        "git checkout {params.commit}; "
        "cd ../..; "
        "cat tmp/plasmidfinder_db/*.fsa > {output}"

rule get_aclame_db:
    output:
        temp(databases_dir + "/aclame_db.fasta")
    params:
        aclame = "https://drive.google.com/file/d/1ipiFom9ha_87k-bNzK_nIM4ABO1l_WP9/view?usp=sharing"
    conda:
        "workflow/envs/download_databases.yaml"
    shell:
        "mkdir -p {databases_dir}; "
        "gdown {params.aclame} --fuzzy -O {output}"

rule get_iceberg_db:
    output:
        temp(databases_dir + "/iceberg_db.fasta")
    params:
        iceberg = "https://drive.google.com/file/d/14MNE_738gIgmQU68y-529DitW_zqxIXy/view?usp=sharing"
    conda:
        "workflow/envs/download_databases.yaml"
    shell:
        "mkdir -p {databases_dir}; "
        "gdown {params.iceberg} --fuzzy -O {output}"

rule get_MGEs_DBs:
    input:
        databases_dir + "/plasmid_finder_db.fasta",
        databases_dir + "/aclame_db.fasta",
        databases_dir + "/iceberg_db.fasta"
    output:
        databases_dir + "/mges_combined.fasta"
    conda:
        "workflow/envs/download_databases.yaml"
    shell:
        "cat {input} > {output}"

############################################################
## Cleans
############################################################

# Need to run a test to see if this will work WILL UPDATE
# Also need to change some things with the config file before deleting all "tmp" files

# onsuccess:
#     shell("rm -f *.csv *.sam *.json *.pdf *_deduplicated.fastq")
#     shell("rm -rf {tmp_dir}")
#     shell("rm -f config.ini")
#     shell("rm -rf {databases_dir}")

# rule clean:
#     shell:
#         """
#         rm -f *.csv *.sam *.json *.pdf *_deduplicated.fastq
#         rm -rf {tmp_dir}
#         rm -f config.ini
#         """

# rule clean_databases:
#     shell:
#         """
#         rm -rf {databases_dir} {tmp_dir}
#         """