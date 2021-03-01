# iSLS9 Workshop Lipidomics Data Analysis

Sphingolipids and T2DM Risk

## Summary

This repository contains datasets and R codes used for 9th International Singapore Lipid Symposium [iSLS9](https://sling.sg/news-events/isls/) workshop held March 4th, 2021. We will be using the data published in *Chew et al., 2019, JCI Insight 5(13)* [10.1172/jci.insight.126925](https://doi.org/10.1172/jci.insight.126925) as an example dataset for this workshop.

In the first part, we will look at the structure of lipidomics datasets and discuss lipid nomenclature systems. We then look at some aspects and challenges of the analytical data processing and quality control workflow. The practical part will include reading of lipidomics data and metadata into R, parsing of the data into different formats for visualization and downstream statistical analyses (see below), converting lipid nomenclature, and creating selected quality control plots.

In the second part, we practice how to reorganize and inspect the overall data trends from both sample meta data and lipidomics data via visualization and dimension reduction. We probe the correlation structure between the sphingolipids and other competing risk factors collected at baseline, and test associations between the lipids and the risk of DM incidence, using logistic regression (binary outcome) and Cox regression analysis (time-to-event analysis).

## Prerequisites

Download the R Project containing the scripts and data used in this workshop from this repository. Alternatively you can clone this repository in RStudio/git. Following packages are needed:

-   Part 1: Following packages will be used: `tidyverse`, `broom` and `rgoslin`. `rgoslin`, used to convert lipid nomenclatures,is only available via github ((<https://github.com/lifs-tools/rgoslin>). See the R notebook (.Rmd) in folder Part_1 on how to install it. Please note that on Windows, you will need [rtools](https://cran.r-project.org/bin/windows/Rtools/), or on macOS [XCode](https://apps.apple.com/sg/app/xcode/id497799835?mt=12), installed in order to install this package.   
    We suggest you to install the currently newest version of RStudio (1.4 or later), which offers a visual RMarkdown/Rnotebook editor.

-   Part 2: TBA

## Authors

-   Bo Burla - [Singapore Lipidomics Incubator \@ National University of Singapore](https://sling.sg)

-   Hyungwon Choi - [Computational & Statistical Systems Biology Laboratory \@ National University of Singapore](https://www.cssblab.org)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

-   Wee Siong Chew
