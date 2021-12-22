# TITAN2D

TITAN2D is a free software application developed by the [Geophysical Mass Flow Group](http://www.gmfg.buffalo.edu/index.php) at the State University of New York (SUNY) at Buffalo. 
TITAN2D was developed for the purpose of simulating granular flows (primarily geological mass flows such as debris avalanches and landslides) over digital elevation models (DEM)s of natural terrain. 
The code is designed to help scientists and civil protection authorities assess the risk of, and mitigate, hazards due to dry debris flows and avalanches. 
TITAN2D combines numerical simulations of a flow with digital elevation data of natural terrain supported through a Geographical Information System (GIS) interface such as GRASS.

TITAN2D was desined for parallel computing on multiple processes/threads, which effectively increases computational power, decreases computing time, and allows for the use of large data sets.

Adaptive Mesh Refienemnt (AMR) allows for the concentration of computing power on regions of special interest.
Mesh refinement captures the complex flow features that occur at the leading edge of a flow, as well as locations where rapid changes in topography induce large mass and momentum fluxes. 
Mesh unrefinement is applied where solution values are relatively constant or small to further improve computational efficiency.

**Reference papers:**
- [Software application](https://www.sciencedirect.com/science/article/pii/S0377027304002288)
- [Modern architecture of TITAN2D](https://link.springer.com/chapter/10.1007/978-3-030-34356-9_10)
- [Modeling & assumptions (a)](https://www.frontiersin.org/articles/10.3389/feart.2020.00275/full)
- [Modeling & assumptions (b)](https://link.springer.com/chapter/10.1007/978-3-319-93701-4_57)


## Modeling Features

TITAN2D currently supports the following rheology models as modeling modules:

- *Mohr-Coulomb*
- *Two-fluid Pitman-Le*
- *Pouliquen-Forterre*
- *Voelmy-Salm*

Please see the [`User Guide`](https://github.com/TITAN2D/titan2d/blob/master/doc/Titan2D_User_Guide.pdf) to learn how to use these modeling features.


## Installation

See [`INSTALL.MD`](https://github.com/TITAN2D/titan2d/blob/master/INSTALL.MD) for compilation instruction
See [`User Guide`](https://github.com/TITAN2D/titan2d/blob/master/doc/Titan2D_User_Guide.pdf) for installation and usage


## Execution

In order to run TITAN2D in `openMP` mode (prefered at this time):

    source <path_to_titan>/bin/titanvars.sh
    <path_to_titan>/bin/titan -nt <number of threads> <input python script>

In order to run TITAN2D in hybrid `MPI/OpenMP` mode:

    source <path_to_titan>/bin/titanvars.sh
    mpirun -np <number of processes> <path_to_titan>/bin/titan -nt <number of threads> <input python script>

## Cite TITAN2D

If you use TITAN2D for your academic research, you are highly encouraged to cite the following publications:

```
@article{titan2d2005,
  title={Parallel adaptive numerical simulation of dry avalanches over natural terrain},
  author={Patra, Abani K and Bauer, Andrew C and Nichita, CC and Pitman, E Bruce and Sheridan, Michael F and Bursik, M and Rupp, Byron and Webber, A and Stinton, AJ and Namikawa, LM and others},
  journal={Journal of Volcanology and Geothermal Research},
  volume={139},
  number={1-2},
  pages={1--21},
  year={2005},
  publisher={Elsevier}
}
```

```
@inproceedings{titan2d2019,
  title={Modernizing Titan2D, a parallel AMR geophysical flow code to support multiple rheologies and extendability},
  author={Simakov, Nikolay A and Jones-Ivey, Renette L and Akhavan-Safaei, Ali and Aghakhani, Hossein and Jones, Matthew D and Patra, Abani K},
  booktitle={International Conference on High Performance Computing},
  pages={101--112},
  year={2019},
  organization={Springer}
}
```

```
@article{titan2d2020,
  title={Comparative Analysis of the Structures and Outcomes of Geophysical Flow Models and Modeling Assumptions Using Uncertainty Quantification},
  author={Patra, Abani and Bevilacqua, Andrea and Akhavan-Safaei, Ali and Pitman, E Bruce and Bursik, Marcus and Hyman, David},
  journal={Frontiers in Earth Science},
  volume={8},      
  pages={275},     
  doi={10.3389/feart.2020.00275},      
  issn={2296-6463}, 
  year={2020},
  publisher={Frontiers}
}
```
