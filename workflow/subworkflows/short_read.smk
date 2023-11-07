# Copyright (c) Boucher Lab. All rights reserved.
# Licensed under the GNU license.
# See LICENSE file in the repository root for full license information.

################################################################################
#                                Imports                                       #
################################################################################

import os, ast, glob

################################################################################
#                           Short Reads Steps                                  #
################################################################################

############################## Assembly ########################################

rule meta_spades_assembly:
    input:
        f_reads = samples_dir + "/{sample_name}_R1.fasta",
        r_reads = samples_dir + "/{sample_name}_R2.fasta"
    output:
        "{sample_name}" + config["EXTENSION"]["DEDUPLICATED"]
    params:
        kmers = confgig["SPADES"]["KMERS"],
        memory = config["SPADES"]["MEMORY"],
        phred = config["SPADES"]["PHRED"],
        outdir = tmp_dir + "/spades_files/{sample_name}/",
        reads = tmp_dir + "/spades_files/{sample_name}/scaffolds.fasta"
    conda:
        workflow.basedir + "/" + config["CONDA"]["ASSEMBLY"]
    threads:
        config["SPADES"]["THREADS"]
    shell:
        "spades.py --meta "
        "--phred-offset {params.phred} "
        "-t {threads} "
        "-m {params.memory} "
        "-k {params.kmers} "
        "-1 {input.f_reads} "
        "-2 {input.r_reads} "
        "-o {params.outdir}; "
        "cp {params.reads} {output}"

############################## Read Lengths ####################################

rule read_lengths:
    input:
        "{sample_name}" + config["EXTENSION"]["DEDUPLICATED"]
    output:
         "{sample_name}" + config["EXTENSION"]["READS_LENGTH"]
    params:
        read_lengths_script = workflow.basedir + "/" + config["SCRIPTS"]["READS_LENGTH"]
    conda:
        workflow.basedir + "/" + config["CONDA"]["ASSEMBLY"]
    shell:
        "python {params.read_lengths_script} {input} > {output}"