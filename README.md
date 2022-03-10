
# WimtrapWeb

**Wimtrap** is a R shiny application for the prediction of
condition-specific transcription factor binding sites (TFBS) in plant
species. The tool locates the potential TF binding-sites by
pattern-matching and annotate them with genomic features that
characterize their genomic context (opening of the DNA, chromatin marks,
DNA conservation, digital genomic footprints,…). Then, the potential TF
binding sites likely to be CARE in the condition of interest are
selected based on a ‘decision rules’ model. This model is otained by
maching learning using ChIP-seq data as reference.

Prebuilt-models and genomic features integrated into the website allow
to quickly obtain prediction of CARE and related potential gene targets
in 10 different conditions *Arabidopsis thaliana* and 2 different
conditions for *Solanum lycopersicum*. Go to the ‘**Query**’ panel to
take advantage of these ressources.

The ‘**Build**’ and ‘**Predict**’ panels allow you to build your own
models and use them for any other condition and organism.

## Installation

Wimtrap is an R Shiny Application that requires the last version of R (R
4.0.4), BiocManager and remotes to be installed.

In R, type the following lines:

``` r
    if(!require("remotes", quietly = TRUE)){  
        install.packages("remotes")
        }
    if(!require("BiocManager", quietly = TRUE)){  
        install.packages("BiocManager")
        }
```

Then, you can enter:

``` r
  BiocManager::install("RiviereQuentin/WimtrapWeb",                     
    dependencies = TRUE,                     
    build_vignettes = TRUE,
    force = TRUE)
```

## Example

To start to use WimtrapWeb, please run the app as follows

``` r
library(WimtrapWeb)
run_app()
#> PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.
```

<div style="width: 100% ; height: 400px ; text-align: center; box-sizing: border-box; -moz-box-sizing: border-box; -webkit-box-sizing: border-box;" class="muted well">Shiny applications not supported in static R Markdown documents</div>

## More informations

Please follow the tutorial video to start using WimtrapWeb:
![](https://www.youtube.com/watch?v=6371fN7dkak)

Methodological details are described in the user guide [user
guide](https://htmlpreview.github.io/?https://github.com/RiviereQuentin/Wimtrap/blob/main/vignettes/Wimtrap.html)
of the cognate Wimtrap R package and in the manual pages of the
functions, which can be accessed by entering in R:

        library(Wimtrap)
        ?importGenomicData()
        ?getTFBSdata()
        ?buildTFBSmodel()
        ?predictTFBS()
        ?plotPredictions()

\##Authors

This web-application was developed at the [LPGMP](https://lpgmp.ulb.be/)
and [MLG](https://mlg.ulb.ac.be/wordpress) groups of the Université
Libre de Bruxelles by Quentin Rivière <qri@hotmail.be>, Madalina
Ciortan, Massimiliano Corso, Grégoire Noël, Nathalie Verbruggen and
Matthieu Defrance <Matthieu.Defrance@ulb.be>.
