# Proteomics final project: six-network comparison

This project extends the summer poster by Allison Burnett and Colin Dempsey,
supervised by Dr. Daniel Malinsky. It compares Fisher-z PC with the original
Pearson + WGCM PC pipeline at alpha = 0.01, 0.05, and 0.10.

## Start in RStudio

1. Extract the entire ZIP. Open the extracted `proteomics_project` folder as an
   RStudio project (File > New Project > Existing Directory).
2. Run `source("setup.R")` once in the Console to install dependencies.
3. Run `source("analysis.R")`. This fits all six networks, saves each completed
   fit, and writes comparison CSV files. WGCM may take time; runtime has not been
   measured in this environment. Do not start a second run while one is active.
4. Open `index.qmd` and click Render. The report reuses saved matching fits and
   creates `docs/index.html`, including three visualizations.
5. Review `results/model_warnings.csv` and the figures before writing the final
   discussion. Send the error message if a fit fails. Completed fits are retained.

RStudio must have Quarto available. Setup requires internet access; compiling R
packages on Windows may require Rtools. No Rgraphviz installation is needed by
this report. Quarto: https://quarto.org/docs/get-started/ .

## Important method detail

Both supplied WGCM scripts use `cor.test` (Pearson) when S is empty, and use
`wgcm.fix` with GAM regression, `weight.num = 7`, `weight.meth = "sign"` when
conditioning variables exist. The report preserves and explicitly labels this
hybrid pipeline. A nonlinear but uncorrelated relationship can be removed at the
Pearson screening stage.

For each fit, the seed is reset to 15; WGCM's simulation count is explicitly 499.
All protein columns are standardized first. This, the per-fit seed reset, and
different installed package versions can cause differences from the poster.

## Files and outputs

- `index.qmd`: Quarto report and three plots, with explanations for a general reader.
- `analysis.R`: six fits, checkpoints, rankings, edge counts, and overlap tables.
- `R/network_helpers.R`: adapted original WGCM wrapper with p-value validation,
  plus consistent undirected-edge calculations.
- `setup.R`: dependency installation.
- `check_helpers.R`: optional focused checks of the wrapper and graph metrics,
  using only base R. Run `source("check_helpers.R")` before the full analysis.
- `data/`: supplied Excel file, unchanged.
- `reference_scripts/`: all six supplied summer scripts, unchanged; not executed.
- `results/`: generated network objects, six adjacency CSVs, degree/rank tables,
  same-alpha overlap table, warning audit, and three exported PNG figures.

Checkpoints are reused only if data, analysis/helper code, settings, R version,
and key package versions match. Changing the report's prose or plot styling does
not refit the networks. To force all fits to run again, delete only the generated
`results/checkpoints` folder and rerun `analysis.R`.

## Finish and publish later

After reviewing actual results, write the final interpretation in `index.qmd`
and render again. Keep attribution to the summer collaborators and describe the
new comparisons as the class project's contribution. The site is not published.
The generated `docs` folder is configured for GitHub Pages. After a final render,
commit the source and `docs` to a repository, then configure Pages to deploy from
your default branch's `/docs` folder. Only publish data you may redistribute.

## Validation status

The uploaded scripts and method settings were reviewed. Input dimensions, numeric
types, missingness, finite values, and nonconstant columns were checked previously.
R and Quarto are unavailable in this environment. The updated R code and the
included helper checks have not been executed here, and no fitted results or
rendered website are included. Conclusions must follow a successful local run.
