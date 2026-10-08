# origin

Forward-time simulation of growing cell populations carrying lineage
recorders. `origin` grows a population in continuous time with an exact
Gillespie process, then replays recorder state along the branches of the
resulting lineage: base editors, Cas9 nuclease arrays, prime-editing mark
alphabets, PALINCODE cBits, mitochondrial variants and non-Mendelian ecDNA.
It writes the true tree alongside the character matrices a reconstruction
method would actually see, which is what makes it useful for asking how much
of a lineage a given recorder can recover.

## Install

The package lives in the `origin/` subdirectory of this repository, not at its
root, and needs nothing beyond `Matrix`.

```r
install.packages("remotes")
remotes::install_github("mckennalab/origin", subdir = "origin")
```

The repository is private, so `remotes` needs a GitHub token; see
`?gitcreds::gitcreds_set`. From a clone instead:

```bash
git clone https://github.com/mckennalab/origin.git
cd origin
R CMD INSTALL origin
```

Running `R CMD INSTALL origin` anywhere other than the repository root gives
`normalizePath(path) : path[1]="origin": No such file or directory`.

## A first simulation

Grow a small population, record on three barcode integrations, and draw the
true lineage. This runs in about a second.

```r
library(origin)

params <- origin_params(sim_length = 6, founders = 3, death = 0.3,
                        barcode_length = 60, integrations = 3,
                        be_rate = 0.15, targets = 12)

run_gillespie_lineage_pipeline(
  params, params_path = NULL, output_dir = "run",
  overrides = list(progress = FALSE, seed = 1, modalities = "barcode"))

# The true lineage, and the characters a recorder hands to tree building.
tree <- ape::collapse.singles(ape::read.tree("run/physicell_lineage_sampled.nwk"))
scores <- as.matrix(data.table::fread("run/barcode_binary_score_matrix.csv.gz")[, -1])
informative <- sum(apply(scores, 2, function(x) length(unique(x))) > 1)

cat(sprintf("%d cells, %d informative characters, %.1f edits per cell\n",
            ape::Ntip(tree), informative, mean(rowSums(scores != 0))))
#> 221 cells, 29 informative characters, 12.3 edits per cell

plot(tree, show.tip.label = FALSE, type = "fan", edge.color = "#1F3B4D")
```

`ape` and `data.table` are used for reading and plotting here; the package
itself does not require them.

The run directory holds the true tree (`physicell_lineage_sampled.nwk`), the
per-cell character matrix, and a manifest of what was simulated. Reconstructing
a tree from `scores` and comparing it with `tree` is the whole experiment the
rest of this repository automates.

## Seeing the edits on the tree

The point of a recorder is that edits are inherited, so cells sharing an edit
should sit together on the lineage. Drawing the character matrix beside the
tree, with rows in the order the tree draws them, shows whether that holds.
Continuing from the run above:

```r
tree <- ape::ladderize(tree)

scores <- data.table::fread("run/barcode_binary_score_matrix.csv.gz")
cells <- scores[[1]]
scores <- as.matrix(scores[, -1])
rownames(scores) <- cells

# Keep the characters that vary; an all-zero or all-one column says nothing.
varies <- apply(scores, 2, function(x) length(unique(x))) > 1
informative <- scores[, varies, drop = FALSE]

# Tips in the order the tree draws them, so the two panels line up.
tips <- tree$edge[, 2]
on_plot <- tree$tip.label[tips[tips <= ape::Ntip(tree)]]
edits <- informative[on_plot, , drop = FALSE] != 0

layout(matrix(1:2, nrow = 1), widths = c(1, 1.6))
par(mar = c(4, 1, 2, 0))
plot(tree, show.tip.label = FALSE, edge.color = "#1F3B4D")
title("lineage", adj = 0)

par(mar = c(4, 0, 2, 2))
image(x = seq_len(ncol(edits)), y = seq_len(nrow(edits)), z = t(edits),
      col = c("#F2F2F2", "#C05746"), axes = FALSE, xlab = "character", ylab = "")
title("edits", adj = 0)
```

The edits fall into blocks that follow the clades: an edit acquired in one cell
is carried by all of its descendants. That block structure is the lineage
signal, and how cleanly it separates is what determines whether a reconstruction
method can recover the tree. Where blocks overlap without nesting, the same
character was edited independently in two places -- homoplasy, which is what
limits most real recorders.

## Documentation

| | |
|---|---|
| [`SIMULATION_GUIDE.md`](SIMULATION_GUIDE.md) | what the simulator models, every recorder, the JSON parameter file, output layout |
| `vignette("gillespie-quickstart", package = "origin")` | building parameters, registering a recorder, and the two relations the parameter values do not predict |
| [`REPO_LAYOUT.md`](REPO_LAYOUT.md) | what lives where, and which entry point to use |
| [`FUNCTION_REFERENCE.md`](FUNCTION_REFERENCE.md) | the contract of every named function |
| [`analysis/runs/README.md`](analysis/runs/README.md) | reproducing the large-scale benchmark runs |
| [`BUG_AUDIT.md`](BUG_AUDIT.md) | defect audit and remaining model-level limitations |

Two things are worth reading before trusting a calibration: the configured
editing rate is nominal rather than realised, and population growth is steeply
non-linear in the death probability. Both are measurable with
`origin_rate_response()` and `origin_growth_response()`, and both are covered
in the vignette.
