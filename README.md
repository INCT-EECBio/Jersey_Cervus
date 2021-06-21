# Quantitative genetics of extreme insular dwarfing: The case of red deer on Jersey

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains codes of the [paper](https://onlinelibrary.wiley.com/doi/10.1111/jbi.14109) published in the Journal of Biogeography. 

## How to cite

> Diniz-Filho, J.A.F., Santos, A.M.C., Barreto, E., Naves, F., Santos, W., Souza, K.S., Santos-Silva, R., Dobrovolski, R., Soares, T.N., Tidon, R., Spigoloni, Z.A., Rangel, T.F., Raia, P., Hortal, J. and Jardim, L. (2021), Quantitative genetics of extreme insular dwarfing: The case of red deer on Jersey. J Biogeogr. [https://doi.org/10.1111/jbi.14109](https://doi.org/10.1111/jbi.14109)

## Data  

This study simulates an individual-based quantitative genetics model parameterized with red deer life-history data gathered from the literature ([more details](https://onlinelibrary.wiley.com/doi/10.1111/jbi.14109)). All data are available in [DRYAD](https://datadryad.org/stash/dataset/doi:10.5061/dryad.47d7wm3cf).

## Analyses

The folders contain files to reproduce each simulation scenarios (1) Baseline, (2) constant_k, (3) no_plasticity, (4) no_recolonization and (5) no_reproduction_system. There is also a file named island_area.txt, which contains data of the island area (Kisl) and isolation (migr) changes through time.

### Repository structure

```bash
Jersey_cervus/
├── baseline  
# The baseline model includes non-constant carrying capacity,phenotypic plasticity, island recolonization,
# and sexual dimorphism.
│   ├── AdaptSS_cervus_K.R
│   └── parallel_function_cervus.R
│   
│   
├── constant_k
# Scripts to run baseline models with constant carrying capacity through time.
│   ├── AdaptSS_cervus_K.R
│   └── parallel_function_cervus.R
│   
├── no_plasticity
# Scripts to run baseline models without phenotypic plasticity.
│   ├── AdaptSS_cervus_K_no_plasticity.R
│   └── parallel_function_cervus.R
│   
├── no_recolonization
# Scripts to run baseline models without island recolonization.
│   ├── AdaptSS_cervus_K_no_recolonization.R
│   └── parallel_function_cervus.R
│   
├── no_reproduction_system
# Scripts to run baseline models without sexual dimorphism.
│   ├── AdaptSS_cervus_no_reproduction_system.R
│   └── parallel_function_cervus.R
│   
├── island_area.txt
├── LICENSE
└── README.md

```

* AdaptSS_cervus_(scenario).R - Scripts to run adaptation through generations. The function AdaptSS runs a simulation of trait adaptation for each generation step and the function run_generation updates simulation's parameters for each generation and calls AdaptSS.

* parallel_function_cervus.R - Set simulation parameters and run simulations.

### Glossary

* h2 - heritability

* cv - phenotypic coefficient of variation

* vm - mutational variance

* Ancestral - Initial trait

* meanK - mean carrying capacity

* Ni - initial population size 

* Nrecol - number of recolonization

* Precol - probability of colonization

* w2 - Length of adaptive landscape

* time_adap - time until adaptation

* selgrad - selection gradient

* meanG - mean genotype

* varG - genotypic variance

* N_end - Final population size

* F_Peak - adaptive landscape optimum (optimum trait)

* meanP - Final mean trait

* bp.max - maximum phenotypic plasticity

* bp - phenotypic plasticity

## Contributing

Contributions are welcome. Open an issue to discuss the changes. 

## Funding

* Ministério da Ciência, Tecnologia e Inovações 

* Conselho Nacional de Desenvolvimento Científico e Tecnológico

* Fundação de Amparo à Pesquisa do Estado de Goiás

* Spanish MICIU Juan de la Cierva-Incorporación

