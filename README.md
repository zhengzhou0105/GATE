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

- Status: Under review by *Environmental and Molecular Mutagenesis*

- Manuscript link: [DOI when available]

- Recommended citation (BibTeX):

  ``` bibtex
  @article{Zhou2026,
    title   = {Optimizing In Vitro Battery for Genotoxicity Testing Using a Probabilistic Network To Support Decision Making},
    author  = {Zhou, Zheng and Pennings, Jeroen L.A. and Sahlin, Ullrika},
    year    = {2026},
    journal = {Environmental and Molecular Mutagenesis},
    doi     = {}
  }
  ```

## Using the deployed app

The easiest way to use GATE is through the public Shiny app, which requires no local installation.

Open the app here: - <https://7t2hvl-zheng-zhou.shinyapps.io/GATE/>

### What the app returns

The app returns the posterior probabilities of the two genotoxicity endpoints — gene mutation, P(GM), and chromosome aberration, P(CA) — and evaluates each of them against a conclusive decision threshold. Optionally, it also returns the probability of rodent carcinogenicity, P(RC). P(RC) is provided as a complementary output rather than as a primary decision target, and is flagged in the app because it is currently less stable than the endpoint-specific probabilities (see Section 4.3.1 of the manuscript).

### Two entry points

The app supports two entry points that mirror the two scenarios of the adaptive framework in the manuscript:

- **Starting assays.** When no data have been collected for the five in vitro assays (Ames, MLA and HPRT for gene mutation; MNvt and CAvt for chromosome aberration), the app recommends a starting assay for each endpoint and returns the endpoint probabilities conditional on the user-entered outcomes.

- **Existing data and follow-up.** When some assay results are already available, the app returns the endpoint probabilities given the current evidence, and — when at least one endpoint remains inconclusive — recommends follow-up assays that would most strongly resolve the remaining uncertainty.

The two entry points connect: after entering the outcomes of the starting assays, the user can move directly to the follow-up view, which carries the starting results forward.

### Independent recommendations

The app produces three independent recommendations, one per decision target:

- **GM recommendation.** Considers only untested GM assays (Ames, MLA, HPRT).

- **CA recommendation.** Considers only untested CA assays (MNvt, CAvt).

- **RC recommendation.** Considers all five assays, and is shown only when P(RC) is requested.

For each recommendation the app returns two candidate assays:

- the assay that would most strongly **confirm** the endpoint if the result were positive, and

- the assay that would most strongly **reject** the endpoint if the result were negative,

together with the resulting endpoint probability in each case. This pairwise framing lets the user choose between a test designed to establish hazard and a test designed to rule it out, according to the decision at hand.

The GM and CA recommendations are endpoint-specific: a GM assay updates P(GM), and a CA assay updates P(CA); neither recommendation proposes an assay targeting the other endpoint. A cross-endpoint recommendation arises only in the RC recommendation, where the assay that most improves P(RC) may target either endpoint.

### Adaptive features

The app implements the adaptive logic of the manuscript rather than a fixed battery:

- **Conclusive thresholds are user-set.** The default is P \> 0.975 for a positive conclusion and P \< 0.025 for a negative conclusion. Lowering or raising the threshold changes the stopping point and the recommended follow-up.

- **The prior for P(RC) is user-set in the existing-data path.** A compound-specific prior can be supplied to reflect in silico, read-across or expert evidence (see Section 4.1 of the manuscript). The endpoint-specific posteriors P(GM) and P(CA) do not depend on this prior.

- **Recommendations update with each new result.** As the user enters additional assay outcomes, the endpoint probabilities and the recommended follow-up are recomputed, so the testing strategy adapts to the evidence accumulated so far.

### App workflow at a glance

1.  Open the app.

2.  Choose an entry point: starting assays, or existing data and follow-up.

3.  Set the conclusive threshold (or keep the default) and, in the existing-data path, optionally supply a prior for P(RC).

4.  Review P(GM) and P(CA) — and P(RC) if requested — together with the conclusive verdict for each.

5.  When an endpoint is inconclusive, use the independent recommendations to identify the most informative follow-up assay.

This route is intended for:

- interactive use,

- reviewers,

- users interested in decision support,

- exploratory analysis without a local R setup.

## Worked examples

The three examples below illustrate the workflow. Each example records the inputs entered, the posterior probabilities returned, the threshold verdict, and the recommendations produced. Values shown are from the current app build; entering the same inputs in the deployed app returns the same numbers.

### Example 1 — Starting from scratch (no data)

- **Input.** No assay has been performed; the user selects "No data has been collected from these five assays".

- **Starting recommendation.** The app recommends HPRT for gene mutation and CAvt for chromosome aberration, as the pair that maximizes the information gained across the two endpoints.

- **Outcomes entered.** HPRT positive, CAvt positive.

- **Returned probabilities.** The app returns P(GM) and P(CA) conditional on these two positive results. Because HPRT and CAvt are single assays per endpoint, the returned probabilities typically remain between the conclusive thresholds, so both endpoints are reported as **inconclusive**.

- **Recommendations.** The GM recommendation proposes a further GM assay (Ames or MLA) as the next step, with two candidates — one to confirm, one to reject. The CA recommendation proposes the second CA assay (MNvt). This illustrates the adaptive logic: starting assays open the analysis, and follow-up narrows the uncertainty that remains.

### Example 2 — A single positive Ames result

- **Input.** Ames positive; all other assays NA.

- **Returned probabilities.** P(GM) = 0.805 and P(CA) = 0.759.

- **Threshold verdict.** Both endpoints remain between 0.025 and 0.975 and are therefore reported as **inconclusive**. A single Ames result is not sufficient to resolve either endpoint.

- **Recommendations.** The GM recommendation proposes a further GM assay (MLA or HPRT) to confirm or reject gene mutation. Because no CA assay has been tested, the CA recommendation proposes a starting CA assay (MNvt or CAvt) as the most informative next step for that endpoint. The two recommendations are independent: the GM follow-up addresses P(GM), and the CA follow-up addresses P(CA). If P(RC) is requested, the RC recommendation considers all five assays and may point to either endpoint.

### Example 3 — Conflicting endpoint evidence

- **Input.** MLA negative; MNvt positive; all other assays NA.

- **Returned probabilities.** P(GM) = 0.638 and P(CA) = 0.680.

- **Threshold verdict.** Both endpoints remain between the thresholds and are therefore reported as **inconclusive**.

- **Recommendations.** The GM recommendation proposes HPRT, the GM assay that would most strongly change P(GM) from its current value. The CA recommendation proposes CAvt, the second CA assay, to strengthen or overturn the CA evidence from MNvt. Each recommendation is specific to its endpoint, reflecting the mechanistic separation between gene mutation and chromosome aberration used throughout the framework.

### A note on conclusive outcomes

If the entered results are sufficient to push P(GM) above the positive threshold (or below the negative threshold), or P(CA) to the same effect, the app reports the endpoint as confirmed (or rejected) and, when both endpoints are conclusive, offers a "Finish" action rather than a follow-up recommendation. The stopping point depends on the threshold chosen at step 3, which is the feature that makes the workflow adaptive rather than fixed.

## Reproducing the manuscript analysis

This section is for researchers who want to reproduce the analysis, figures, tables and manuscript outputs locally.

The files needed for reproduction are organized as follows:

```         
manuscript/
├── utilities.R
├── data.R
├── modeling.R
├── decision_analysis.R
├── plotting.R
├── net_fullbayes.stan
├── run.R
└── data/
    └── Luijten25_TRAID09.xlsx
```

### System requirements

- R version 4.2 or newer

- package installation access

- Stan to perform Bayesian inference (`cmdstanr` preferred)

- optional: `renv` for reproducible package management

- recommended: a modern operating system with Rtools installed on Windows

### Installation

Clone the repository:

``` bash
git clone https://github.com/zhengzhou0105/GATE.git
cd GATE
```

### Running the analysis

A typical reproduction workflow follows the script order below:

``` bash
Rscript manuscript/run.R
```

The user may also execute the workflow step by step, running `data.R`, `modeling.R`, `decision_analysis.R` and `plotting.R` in that order.

## Contribute and contact

For questions, collaboration or bug reports, please contact the repository maintainer or the corresponding author of the manuscript at:

- Zheng Zhou: [zheng.zhou\@ruokavirasto.fi](mailto:zheng.zhou@ruokavirasto.fi)
