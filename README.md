## Topoclimate buffers floristic diversity from macroclimate in temperate mountain forests.

[![DOI](https://zenodo.org/badge/652591507.svg)](https://zenodo.org/doi/10.5281/zenodo.12626860)

### Abstract:

Microclimates strongly influence the composition and diversity of forest plant communities. Recent studies have highlighted the role of tree canopies in shaping understory thermal conditions at small spatial scales, especially in lowland forests. In mountain forests, however, the influence of topography in environmental conditions (e.g. topoclimate) is ought to also influence plants’ perceived temperature. Understanding how topography and canopies interactively affect understory temperature is key to identifying stable refugia that could shelter cold-adapted forest specialist species under climate change.

Here we report on growing season understory temperatures using 48 loggers in contrasting topographic features of a mid-range mountain valley spanning from 475 m.a.s.l. to 1203 m.a.s.l. in the Vosges Mountains (NE France). We disentangle the relative importance and the effects of topography vs. canopies in determining local temperatures. We then evaluate how topography and canopy-induced variation in temperature drive plant community composition and richness in 306 floristic surveys distributed across the studied mountain valley.

Our results show that topography outweighed canopy cover in explaining growing season understory temperatures. Regardless of canopy, the daily mean temperature of the growing season in south-facing ridges was 1.5 °C (CI: ± 0.88 °C) warmer than shaded valley bottoms, while dense canopies cooled temperatures by 0.5 °C (CI: ± 0.48 °C) compared to open canopies. Topoclimate explained community composition as much as elevation and was the only significant predictor of species richness. Cold topoclimates harbored 30% more species than the average species richness across our plots. This increase in species richness was explained by an increase of cold-adapted species, both forest specialist and generalist species.

*Synthesis*. Our findings highlight a stronger role of topography compared to canopy cover on community composition in mountain forests via topoclimatic cooling of north-facing slopes and valley bottoms. The importance of topographic features to explain temperature cooling and diversity underpins their role as present and future microrefugia.

**DOI:** [10.24072/pcjournal.519](https://doi.org/10.24072/pcjournal.519)


This repo contains the raw dataset and the code used in the study "Topoclimate buffers floristic diversity from macroclimate in temperate mountain forests" currently *under review*

You can clone this R project to reproduce the Analysis, the results and the figures, we advise you to use R studio and its Git connection.

To rerun the analysis, launch the project thanks to the .rproj file and run the `Main_script.r` in the `scripts` folder. There is also a way to reproduce the random sampling and compute the climatic mean of the study region with the scritp in the `additional_scripts` folder.

Tutorial to getting stated with R studio and git :https://jennybc.github.io/2014-05-12-ubc/ubc-r/session03_git.html

The dataset contains microclimatic data also ssubmitted to SoilTemp, remember to credit this study and the SoilTemp database if you plan on using them.
The dataset contains community data originating from AgroParisTech studies, please contact the authors if you plan on using them.

We included subsets of the rasters in order to reproduce the figures
The analysis was run using R `4.1.1`
