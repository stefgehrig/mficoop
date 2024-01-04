In this repository, I share the scientific manuscript **"Societal norms of cooperation and microcredit repayment under group lending"** and all data and analysis code. 
The study looks at the relationship between repayment rates reported by microfinance institutions, their level of group lending and the norms of cooperation prevailing in a given country. 

<p align="center">
<img src="https://github.com/stefgehrig/mficoop/blob/main/outputs/fig_params_main.png" width="600">
</p>


<p align="center">
<img src="https://github.com/stefgehrig/mficoop/blob/main/outputs/fig_predcoop_main.png" width="600">
</p>

<p align="center">
<img src="https://github.com/stefgehrig/mficoop/blob/main/outputs/fig_varcomp_main.png" width="600">
</p>

The manuscript can be found in the `manuscript.pdf` file. The work is a long-term project that was supported by co-authors ([thkyritsis](https://github.com/thkyritsis) and others) and other people who gave valuable advise. 
The decision was taken to not bring this manuscript to peer-reviewed publication. The major reason is that, over the years and through the input of critical reviewers at *World Development* (R&R, withdrawn) and *Oxfod Open Economics* (R&R, withdrawn), 
I have come to realize that with the aggregated data and relatively weak study design at hand, 
it just too difficult to make strong conclusions about an effect of societal norms on repayment behavior of microcredit clients. Multiple sources of uncontrolled confounding
and selection biases could be at work. This is beside our best efforts to specify rich hierarchical models that incorporate many sources of (co-)variation. 
Rather than drawing potentially unwarranted conclusions about policy design and human behavior in the
scientific record, here we simply present the research project's analysis and results for public access. 

## R scripts

To reproduce the analyses, the `modeling.R` script has to be run as a first step. It imports the required data and the helper functions stored in `functions.R`, and runs all Bayesian models. 
The `produce_outputs.R` script uses the modeling outputs to create figures and tables. The `manuscript.Rmd` script renders the final manuscript.
