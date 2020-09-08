# Communication Map based on the single-cell RNA-sequencing data

ComMap builds an inter-cellular communication map (Interactome) based
on the scRNA-seq data and pathway enrichment. The final interactome
contains receptors connected via the intra-cellular signaling pathways
to a set of transcription factors (TF) and their immediate targets
(TF-TF regulatory network) that are maintained throughout a population
of single cells. Thus, persistent TF pattern emulates a "phenotype"
for each population. Receptors signaling to maintain the phenotype of
each cell-type population are then paired with their cognate ligands
via a database of protein-protein interactions.

## Quick start

We have scRNA-seq data (GEO: GSE96981) from co-culturing of three
individual cell lines (mesenchymal (MC), hepatic (HE), and endothelial
(EC)) found in liver followed by an event of self-organization into a
Liver Bud (LB). We could study the change in putative inter-cellular
interactions between the cells.  In the directory where this package
is located, type:
```
source("commap.R")

``` 
This loads all functions from the package. Next, we will load the data
and run the whole pipeline, the result of which is the communication
map. In this experiment, one has the three cell lines (MC, HE, and EC)
individually (prior to LB) and in LB.

```
## Prepare the data
d.indiv <- readRDS("LiverLB_indiv.Rds")
d.lb <- readRDS("LiverLB_lb.Rds")
## Run the pipeline (use genes that are expressed in at least 10% of the population
## Use "old" ligand-receptor DB (this is a more curated one)
cm.indiv <- sc.commap(d.indiv, org="HSA", conserv.thr=0.1, LR.db="old")
cm.lb <- sc.commap(d.lb, org="HSA", conserv.thr=0.1, LR.db="old")
```

One thing that you could observe is that VEGFA->KDR signalling is more
prominent when in LB phase as was observed in the original paper (Camp
et al. Nature 546, 533–538, 2017):
```
> subset(cm.indiv, Receptor == "KDR" & Ligand == "VEGFA")
    Lig.pop    L.frac Ligand Receptor    R.frac Rec.pop
223      EC 0.1621622  VEGFA      KDR 0.3108108      EC
589      HE 0.6548673  VEGFA      KDR 0.3108108      EC
976      MC 0.6346154  VEGFA      KDR 0.3108108      EC
> subset(cm.lb, Receptor == "KDR" & Ligand == "VEGFA")
     Lig.pop    L.frac Ligand Receptor    R.frac Rec.pop
237    EC.LB 0.8679245  VEGFA      KDR 0.4716981   EC.LB
684    HE.LB 0.8703704  VEGFA      KDR 0.4716981   EC.LB
1138   MC.LB 0.9850746  VEGFA      KDR 0.4716981   EC.LB
```
In the LB case one can see that the relative expression
`L.frac/R.frac` is much stronger for both the ligand (`VEGFA`) and
receptor (`KDR`). Also note that all three populations signal toward
endothelial cells.

## Quick start (Enrichment only)

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
duplicated for each cell of that type).

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
al. Sci Rep 8, 5884, 2018 https://doi.org/10.1038/s41598-018-24051-5,
Fig.6E therein)


## Credits and special thanks

1. TF-TF interaction network is from
[SigHotSpotter](https://academic.oup.com/bioinformatics/article/36/6/1963/5614427)

2. Ligand-Receptor interaction were collected by PhD candidates (31.08.2020) Kartikeya Singh and Chrysovalantou Kalaitzidou at
University of Luxembourg.

3. Ligand-Receptor matching using mean-field approximation was written by Sascha Jung (Center of Excellence Severo Ochoa, Bilbao, Spain).

