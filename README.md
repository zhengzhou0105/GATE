# GATE

Genotoxicity Adaptive Evaluation and Testing (GATE) is a repository for:

  - the analysis and visualization code used in the manuscript,

  - the deployed Shiny application for adaptive genotoxicity decision support.

The repository is organized to support two different audiences:

  - users who want to run or explore the deployed app without installation,

  - researchers who want to reproduce the manuscript analysis and outputs locally.

## Related manuscript

This repository is associated with a manuscript on adaptive genotoxicity testing and decision support.

  - Manuscript title: Optimizing An In Vitro Battery for Genotoxicity Testing Using a Probabilistic Network To Support Adaptive Decision Making

  - Status: Under Review

  - Manuscript link: [DOI when available]

- Recommended citation (_bibtex_ format):

   @article{Zhou2026, 
   
    title={Optimizing In vitro Battery for Genotoxicity Testing Using a Probabilistic Network To Support Decision making}, 
    
    author={Zhou, Zheng and Pennings, Jeroen L.A. and Sahlin, Ullrika}, 
    
    year={2026},
    
    Journal = {Environmental and Molecular Mutagenesis},
    
    doi = {}
    
    }


## Using the deployed app

The easiest way to use GATE is through the public Shiny app, which does not require local installation.

Open the app here:
- https://7t2hvl-zheng-zhou.shinyapps.io/GATE/

### App workflow

1. Open the app in your browser.

2. Decide if you want to select starting assays or the next follow-up assays.

3. Update the prior belief if needed.

4. Review the posterior probabilities for the relevant genotoxicity endpoints.

5. Use the app recommendations to identify the most informative next assay when evidence is inconclusive.

This route is intended for:

- interactive use,

- reviewers,

- users interested in decision support,

- exploratory analysis without a local R setup.

## Reproducing the manuscript analysis

This section is for researchers who want to reproduce the analysis, figures, tables, and manuscript outputs locally.

### System requirements

- R version 4.2 or newer

- package installation access

- stan to perform Bayesian inference (`cmdstanr` preferred)

- optional: `renv` for reproducible package management

- recommended: a modern operating system with Rtools installed on Windows

### Installation

Clone the repository:

```bash
git clone https://github.com/zhengzhou0105/GATE.git
cd GATE
```

### Running the analysis

The manuscript-specific workflow is expected to live under the `manuscript/` directory.
Supporting information of the manuscript is available under the `manuscript/docs/` directory.

A typical reproduction workflow may follow the script order below:

```bash
Rscript manuscript/analysis/main.R
```

To reproduce the figures:

```bash
Rscript manuscript/visualization/plotting.R
```

## Contribute and contact

For questions, collaboration, or bug reports, please contact the repository maintainer or the corresponding author of the manuscript at: (zheng.zhou@ruokavirasto.fi)
