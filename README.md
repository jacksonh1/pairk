pairk
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/jacksonh1/pairk/workflows/CI/badge.svg)](https://github.com/jacksonh1/pairk/actions?query=workflow%3ACI)
<!-- [![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pairk/branch/main/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/pairk/branch/main) -->


motif conservation in IDRs through pairwise k-mer alignment

This work was supported by the National Institutes of Health under Award Number R35GM149227. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

### Features
Quantify the relative conservation of a small sequence motif in intrinsically disordered regions (IDRs) of proteins, without the need for a multiple sequence alignment (MSA).

The pairk method:
![pairk method](docs/images/fragpair_cartoon_v2.png)


Example - PairK vs MSA conservation:
<p align="center">
  <img src="docs/images/f1-example_MSA_problems.png" width="500">
</p>


See the demo/tutorial jupyter notebook here: [demo/pairk_tutorial.ipynb](demo/pairk_tutorial.ipynb)

<!-- ![PairK vs MSA conservation](docs/images/f1-example_MSA_problems.png) -->

### Installation
<!-- not yet published to pypi -->

<!-- #### Current recommended installation: -->
<!-- to ensure that you have the correct dependencies (specically the correct version of biopython), -->
<!-- we recommend installing in a conda environment with the provided environment.yml file with the following commands: -->
<!-- ```bash -->
<!-- git clone https://github.com/jacksonh1/pairk.git -->
<!-- cd pairk -->
<!-- conda env create -f=environment.yml -->
<!-- ``` -->
<!-- Then activate the environment with: -->
<!-- ```bash -->
<!-- conda activate pairk -->
<!-- ``` -->
<!-- then install pairk with: -->
<!-- ```bash -->
<!-- pip install . -->
<!-- ``` -->
<!-- or for an editable install that you can modify: -->
<!-- ```bash -->
<!-- pip install -e . -->
<!-- ``` -->

<!-- pip install pairk@git+git://github.com/jacksonh1/pairk.git -->

<!-- #### very near future installation instructions (after publication to pypi): -->

```bash
pip install pairk
```
or for an editable install that you can modify:
```bash
git clone https://github.com/jacksonh1/pairk.git
cd pairk
pip install -e .
```

#### virtual environment installation:

We suggest using a virtual environment to install pairk, such as conda or venv. 
You can create a new environment and just install pairk as above, or you can 
use the provided environment.yml file to create a new environment with the 
necessary dependencies like so:
```bash
git clone https://github.com/jacksonh1/pairk.git
cd pairk
conda env create -f=environment.yml
```
Then activate the environment with:
```bash
conda activate pairk
```
and install pairk with either:
```bash
pip install .
```
or for an editable install that you can modify:
```bash
pip install -e .
```


### Documentation
see the [pairk documentation](https://pairk.readthedocs.io/en/latest/).

Also see our jupyter notebook tutorial in the `demo` folder.


### Copyright
Copyright (c) 2024, Jackson Halpin


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.1.
