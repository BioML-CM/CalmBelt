# CalmBelt
Rapid SARS‑CoV‑2 genome characterisation for outbreak tracking

![2021-05 COVID figures](https://user-images.githubusercontent.com/76929527/127766996-cdd82bb1-4e2c-49cd-b413-822ecf254eb5.png)

CalmBelt's [demo](https://calmbelt.mtms.dev)

Contact: suphavilaic@gis.a-star.edu.sg and hatairat.y@cmu.ac.th


# preparation
### 1. install anaconda and create environment <br>
- wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh <br>
- bash Miniconda3-latest-Linux-x86_64.sh <br>
- conda create -n covid python=3.7 <br>

### 2. install blastn <br>
- sudo apt-get update <br>
- sudo apt-get install ncbi-blast+ <br>

### 3. install iqtree <br>
- `iqtree2` precompiled binary can be downloaded from https://github.com/iqtree/iqtree2/releases/tag/v2.1.2 <br>

### 4. install mafft, pango and requirement package <br>
- conda activate covid <br>
- conda install mafft <br>
- git clone https://github.com/cov-lineages/pangolin.git <br>
- conda env update -n covid --file pangolin/environment.yml <br>
- cd pangolin <br>
- pip install . <br>
- cd .. <br>
- pip install -r requirement.txt <br>
- conda install -c plotly plotly-orca <br>

# preprocessing
- download folder `calmbelt` <br>
- create genome_dir folder in calmbelt folder that contain fasta file and genome_metadata(tsv)  <br>
- python preprecess.py genome_dir genome_metadata reference_fname preprocess_dir [predefined_clade] [predefined_label] [subsampling_n_samples]

# run app
- python covid_app.py preprocess_dir port[8050] <br>
- localthost:8050

