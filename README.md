# Firn thickness, grain size, accumulation model


Code to solve the nondimensional firn compaction model presented in the preprint article *Grain-size evolution controls the accumulation dependence of modeled firn thickness*. 

This paper is unpublished/unreviewed, but will be availible in preprint form soon and a link to the paper will be provided here when it is. 

All matlab code needed to run the model and produce the figures and save them as pngs are included. All the resulting pngs are also provided in the directory *figures*.

All code runs in Matlab 9.7 (and probably most other matlab versions) and does not require any toolboxes.


---

### Downloading and running the code
The fastest way to download and run this code is to clone the repo, go into the directory and run the main plotting script ([plottingScripts.m](https://github.com/ldeo-glaciology/firn-compaction/blob/d54619ee266640860a0d52bda4bb078c01e34a70/plottingScripts.m)), as follows.

In Matlab run

```
! git clone git@github.com:ldeo-glaciology/firn_accumulation_model.git
cd firn_accumulation_model
plottingScripts

```
This takes about 2 minutes on my Mac book air.

By default [plottingScripts.m](https://github.com/ldeo-glaciology/firn-compaction/blob/d54619ee266640860a0d52bda4bb078c01e34a70/plottingScripts.m) will run some of the simulations again, but many of the results will be loaded from .mat files saved in the directory *savedResults*. 

If instead you want to rerun all the simulations change `rerun = 0` to `rerun = 1` in [this](https://github.com/ldeo-glaciology/firn-compaction/blob/d54619ee266640860a0d52bda4bb078c01e34a70/plottingScripts.m#L7) line in [plottingScripts.m](https://github.com/ldeo-glaciology/firn-compaction/blob/d54619ee266640860a0d52bda4bb078c01e34a70/plottingScripts.m), then run the script.

Note that, depending on how you have git setup on your machine you may need to replace the first line above with 

```
! git clone https://github.com/ldeo-glaciology/firn_accumulation_model.git
```

Please let us know if you have any questions or suggestions through issues.

---
### Acknowledgements

Thank you to <https://github.com/foxelas> for `table2latex.m`:

Eleni Aloupogianni (2021). TABLE2LATEX (https://github.com/foxelas/Matlab-assisting-functions/releases/tag/v1.1), GitHub. Retrieved August 18, 2021.
<https://github.com/foxelas/Matlab-assisting-functions/blob/master/table2latex.m>

Thank you to Simon Kerschbaum for `two_point_upwind_uni_D1.m`:

Kerschbaum, Simon. (2020). Backstepping Control of Coupled Parabolic Systems with Varying Parameters: A Matlab Library (1.0). Zenodo. <https://doi.org/10.5281/zenodo.4274740>

Thank you to the US National Science Foundation's Office of Polar Programs for funding this work through grants OPP 19-35438 and PSU 5861-CU-NSF-8934.