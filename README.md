# Firn thickness, grain size, accumulation model

Code to solve the nondimensional firn compaction model presented in the preprint article *Grain-size evolution controls the accumulation dependence ofmodeled isothermal firn thickness*. 

This paper is unpublished/unreviewed, but will be availible in preprint form soon and a link to the paper will be provided here when it is. 

All code needed to produce the figures from the manuscript is contained in the directory *model*. All the resulting figures are provided as pngs.

All code runs in Matlab 9.7 (and probably most other matlab versions) and does not require any toolboxes.

---
Run `plottingScripts.m` to plot and save new versions of all the figures. The code was 

For example, Fig. 2:
![figure 2](exampleFig2.png)
---
###Acknowledgements
Thank you to <https://github.com/foxelas> for `table2latex.m`:

Eleni Aloupogianni (2021). TABLE2LATEX (https://github.com/foxelas/Matlab-assisting-functions/releases/tag/v1.1), GitHub. Retrieved August 18, 2021.
<https://github.com/foxelas/Matlab-assisting-functions/blob/master/table2latex.m>

Thank you to Simon Kerschbaum for `two_point_upwind_uni_D1.m`:

Kerschbaum, Simon. (2020). Backstepping Control of Coupled Parabolic Systems with Varying Parameters: A Matlab Library (1.0). Zenodo. <https://doi.org/10.5281/zenodo.4274740>
