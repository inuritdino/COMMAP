# Communication Map based on the single-cell RNA-sequencing data

ComMap builds an inter-cellular communication map (Interactome) based
on the scRNA-seq data and pathway enrichment.

## Quick start

In the directory where the package files are type at the R prompt:
```
source("commap.R")
```
or use the full path to the `commap.R` file, e.g. `source("/home/user1/lib/R/commap/commap.R")`

To run pathway enrichment, given the data `scData` and experiment in mice:
```
enrich.res <- run.enrichment(scData, org="MMU")
```

#### Note on data
`scData` must be in a matrix format (genes x cells) with `rownames`
being gene names and `colnames` being the cell type annotation (R
allows for duplicated colnames, so the same cell type must be
duplicated for each cell).

#### Note on enrichment output
`enrich.res` is a recursive list structure. Each cell-type list,
contains list for each cell (named as in `scData`). The cell-level
list contains lists corresponding to the pathways
(e.g. KEGG_JAK_STAT_PATHWAY). The pathway-level lists contain lists
with parameters and gene names.

```
CELL-TYPE1
|->CELL1
  |->PATHWAY1
    |->r
    |->p
    |->if.TF
    |-> other pathway-specific fields (see below)
  |->PATHWAY2
  |-> ...
|->CELL2
|-> ...
CELL-TYPE2
```

##### Pathway-specific fields
- `r`: vector of receptors belonging to the pathway and expressed in
this cell
- `p`: P-value of enrichment
- `if.TF`: vector of interface transcription factors (TF), that are
expressed and belong to the pathways
- `n.compat.t`: number of compatible targets for each if.TF
above. Target is under immediate effect of one of the if.TF's
according to the TF-TF interaction network. Target is compatible when
an if.TF activates it and it is expressed, or when it is inhibited and
is down.
- `n.pos.compat.t`: # of target TF's that are compatible with
  activation from an if.TF.
- `n.all.t`: # of all target TF's
- `compat.t`: a list with compatible targets for each if.TF
- `pos.compat.t`: a list with "positively" compatible targets for each
  if.TF (subset of `compat.t`).
- `cell.id`: cell id (ordering number).
- `p.adj`: adjusted P-value for enrichment

##### Example enrichment

We will take a data set (GEO: GSE60749) where different conditions for
mouse embryonic stem cells were used. This data set does not contain interacting populations
so we cannot get an interaction map, but compare the enrichment results.
While in the source directory:

```
d <- readRDS("mESC_expr_all_conditions.Rds")
enrich <- run.enrichment(d, org="MMU")
## Get all if.TF's in the LIF+serum media (cell-type mESC)
ifs.gr <- unlist(lapply(enrich$mESC,function(x) lapply(x,function(y) y$if.TF)))
## Get all if.TF's in the 2i+LIF media (cell-type TwoiLIF)
ifs.2i <- unlist(lapply(enrich$TwoiLIF,function(x) lapply(x,function(y) y$if.TF)))
## Difference in TF's
setdiff(ifs.gr,ifs.2i)
# [1] "Myc"    "Snai1"  "Aire"   "Jun"    "Gli1"   "Id2"    "Smad7"  "Id1"   
# [9] "Id4"    "Smad1"  "Id3"    "Rxrg"   "Nfyb"   "Rfxank" "Nfya"   "Nfyc"
```

As we can see that Id1-4 and Myc factors of stemness are differential
enriched for the LIF+serum media as was noted previously (Ghimire et
al. Sci Rep 8, 5884, 2018 https://doi.org/10.1038/s41598-018-24051-5)



## Credits and special thanks

1. TF-TF interaction network is from
[SigHotSpotter](https://academic.oup.com/bioinformatics/article/36/6/1963/5614427)

2. Ligand-Receptor interaction were collected by PhD candidates (31.08.2020) Kartikeya Singh and Chrysovalantou Kalaitzidou at
University of Luxembourg.

3. Ligand-Receptor matching using mean-field approximation was written by Sascha Jung (Center of Excellence Severo Ochoa, Bilbao, Spain).

