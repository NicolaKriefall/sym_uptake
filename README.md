#sym_uptake
Scripts used in Ali et al. (in review): "Reef sediment serves as an essential source of algal symbionts for horizontally-transmitting coral recruits"

About the "Total script.R" file:
We began with .sff files from 454 Roche sequencing output of the ITS2 region (see Ali et al. paper https://www.biorxiv.org/content/biorxiv/early/2018/09/20/421339.full.pdf for primer sequences), which we then converted to .fastq files for further analysis.

The dada2 pipeline used in our script was modified from:
https://benjjneb.github.io/dada2/tutorial.html

After our analysis was complete, the dada2 authors released an ITS-specific pipeline, found here for your benefit:
https://benjjneb.github.io/dada2/ITS_workflow.html

The output file from the dada2 analysis is provided in this repository, and further output files can be found in the Ali et al. paper supplements. The input files for the uptake plot & heat map are also provided in this repository. 

Please contact Nicola Kriefall [thenicolakriefall(at)gmail.com] with any questions!
