# VeloSim

**VeloSim** is a simulation package that generates single-cell RNA Sequencing data for continuous cell developmental process. Given the cell developmental trajectory backbone, e.g. tree, cycle, cycle-tree, etc, **VeloSim** is able to simulate the whole dynamics of mRNA molecule generation, produces **unspliced mRNA count matrix**, **spliced mRNA count matrix** and **RNA velocity** at the same time. VeloSim can be easily used to benchmark **trajectory inference** and **RNA velocity inference** methods.



### Installation

You can install VeloSim using following command:

```r
devtools::install_github("PeterZZQ/VeloSim")
```



### Vignettes

To learn more about using VeloSim for simulation, please check out three example vignettes below:

* Simulation for cell cycle trajectory
* Simulation for tree-like trajectory
* Simulation for Cycle-tree trajectory

