# TELSeq
A workflow for contextualization of antibiotic resistance in microbiomes.

### Requirements 
We manage our dependencies trough [`conda`](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Install Snakemake and Clone the Repository

Create the environment for `telseq` and install snakemake.

```bash
conda create -c conda-forge -n mamba_base mamba
conda activate mamba_base
mamba create -c conda-forge -c bioconda -n telseq snakemake git git-lfs
mamba activate telseq
```

Clone the repository.

```bash
mamba activate telseq
git lfs install
git clone https://github.com/marco-oliva/argmobrich_analysis.git
```

### Usage on local desktop

The `telseq` workflows assumes that the `fastq` files will be stored in a directory called `samples` in the working directory. Here we show the usage and the directories structure that can be used with the default `config.json` file.

```bash
cd argmobrich_analysis
mamba activate telseq

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq[.gz] work_dir/samples

# Run the workflow
snakemake -j <number of threads available> --configfile config.json --use-conda --conda-frontend mamba
```

#### Plots
```bash
# Spawn the jupyter notebook for the violin plot
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook violin_plot_all_samples.pdf

# Spawn the jupyter notebook for the resitome heatmap
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook heatmap_all_samples.pdf

# Spawn the jupyter notebook for the colocalizations plot for a specific sample
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook <your sample name>_colocalizations_plot.pdf"
```

### Usage on slurm cluster

Edit the `cluster.json` file in order to fit your resources.

```bash
cd argmobrich_analysis

# Create the directories structure
mkdir -p work_dir/samples work_dir/logs 

# Move the fastq files in the samples directory
mv <your_data>.fastq[.gz] work_dir/samples

# Run the workflow
sbath run.sh
```

#### Plots
Here remember that this will run on a node in the cluster so you will need an ssh tunnel to the node to be able to connect to the jupyter notebook.

```bash
# Spawn the jupyter notebook for the violin plot
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook violin_plot_all_samples.pdf

# Spawn the jupyter notebook for the resitome heatmap
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook heatmap_all_samples.pdf

# Spawn the jupyter notebook for the colocalizations plot for a specific sample
snakemake --configfile config.json -j1 --use-conda --conda-frontend mamba --edit-notebook <your sample name>_colocalizations_plot.pdf
```

