# CalmBelt
Rapid SARS‑CoV‑2 genome characterisation for outbreak tracking

![2021-05 COVID figures](https://user-images.githubusercontent.com/76929527/127766996-cdd82bb1-4e2c-49cd-b413-822ecf254eb5.png)

CalmBelt's [demo](https://calmbelt.mtms.dev)

Contact: suphavilaic@gis.a-star.edu.sg and hatairat.y@cmu.ac.th


# preparation
1. install anaconda and create environment
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \n
bash Miniconda3-latest-Linux-x86_64.sh
conda create -n covid python=3.7

2. install blastn
sudo apt-get update
sudo apt-get install ncbi-blast+

3. install iqtree
`iqtree2` precompiled binary can be downloaded from https://github.com/iqtree/iqtree2/releases/tag/v2.1.2

4. install mafft, pango and requirement package
conda activate covid
conda install mafft

git clone https://github.com/cov-lineages/pangolin.git
conda env update -n covid --file pangolin/environment.yml
cd pangolin
pip install .
cd ..
pip install -r requirement.txt
conda install -c plotly plotly-orca
pip install gunicorn
pip install notebook
