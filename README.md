# MCF7-ROS-RNAseq
This is a script I used to make an R Shiny application that allows for browsing the time-course gene expression signal, and associated statistics, for an experiment where MCF-7 cells were exposed to oxidative stress-inducing compounds.  This dataset was published in *Free Radical Biology and Medicine*, article found here: [Levings et al. *FRBM*. 2021](https://doi.org/10.1016/j.freeradbiomed.2021.05.016). 


To use this app-

1) Clone this GitHub repository to [desired location].
2) Make sure you have the R programming language and the following packages installed: *magrittr*, *ggplot2*, *ggrepel*, *ggtext*, *gridExtra*, *scales*, *shiny*, *tidyverse*, *xtable*
3) Run the following code in your terminal: ```R -e "shiny::runApp('~/[desired location]/MCF7-ROS-RNAseq/app.R')"```
4) R will run the app, and show what address the plot will be output at.  It should look something like: ```Listening on http://127.0.0.1:6880``` 
5) Copy the address and paste it into your web browser to use the app.

---------------------------------------------------------------------------------------------------

**Raw data sources:**

Gene Expression Omnibus ([GEO]( http://www.ncbi.nlm.nih.gov/geo/)): [GSE173664](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173664)

---------------------------------------------------------------------------------------------------

Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
