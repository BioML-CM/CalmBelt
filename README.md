# CalmBelt
Rapid SARS‑CoV‑2 genome characterisation for outbreak tracking

![2021-05 COVID figures](https://user-images.githubusercontent.com/76929527/127766996-cdd82bb1-4e2c-49cd-b413-822ecf254eb5.png)

CalmBelt's [demo](https://calmbelt.mtms.dev)

Contact: suphavilaic@gis.a-star.edu.sg and hatairat.y@cmu.ac.th


## How to run CalmBelt
### 1. Create conda environment
In this example, we use the lightweight `miniconda`.
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda create -n covid python=3.7
```

### 2. Installing required packages
A precompiled binary for `iqtree2` can be downloaded from https://github.com/iqtree/iqtree2/releases/tag/v2.1.2. The remaining packages can be installed as follow:
```bash
conda activate covid 
conda install -c bioconda blast
conda install mafft 

git clone https://github.com/cov-lineages/pangolin.git 
conda env update -n covid --file pangolin/environment.yml 
cd pangolin 
pip install .
cd .. 

git clone https://github.com/neherlab/treetime.git
pip install .
cd ..

pip install -r requirement.txt 
conda install -c plotly python-kaleido
```
Toytree will require `gs`, so you might need to install ghostscript as well `sudo apt-get install ghostscript`.

### 3. Preprocessing your input genome set
This step would require users to prepare input files as follow:
- `genome_dir` contain fasta files
- `genome_metadata.tsv` contains region, clade, and lineage information
- `world_metadata.tsv` contains continent, clade, and lineage information of samples around the world (can be obtained from GISAID).
```bash
python preprecess.py genome_dir reference_fname reference_gene_loci predefined_clade predefined_label preprocess_dir world_metadata_fname country_name subsampling_n_samples
```
For example (an example dataset is available at `calmbelt/preprocess.zip`),
```bash
python preprocess.py genome_dir reference.fasta gene.tsv clade.tsv name_by_who.tsv preprocess world_metadata.tsv Singapore -1
```

### 4. Starting CalmBelt web application

```bash
python covid_app.py preprocess 
```
CalmBelt is based on Dash/Flask framework (default port = 8050), and users can access the website at http://localthost:8050. Additional information for deploying CalmBelt on a remote server is avilable upon request.

