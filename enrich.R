### Pathway enrichment specific to the commap package

## Original RTF pairing function: purely R code with parallelism
sc.RTF.pairing <- function(expTable.b){
### Get the list of meaningful Receptors, that is expressed, pathway-enriched etc.. in a pop
### expTable -- boolean gene expression table for a single population
    suppressPackageStartupMessages({
        require(parallel)
    })
    ## Take only KEGG for the moment
    ##PTHWGENES_ <- PTHWGENES[grep("REAC",names(PTHWGENES))]
    PTHWGENES_ <- PTHWGENES[grep("KEGG",names(PTHWGENES))]
    ##PTHWGENES_ <- PTHWGENES[union(grep("KEGG",names(PTHWGENES)),grep("REAC",names(PTHWGENES)))]
    ##PTHWGENES_ <- PTHWGENES[1:20]
    ## PTHWGENES_ <- PTHWGENES[grep("KEGG_ERBB_SIGNALING_PATHWAY",names(PTHWGENES),fixed=T)]
    ## PTHWGENES_ <- PTHWGENES[grep("KEGG_WNT_SIGNALING_PATHWAY",names(PTHWGENES),fixed=T)]
    ## PTHWGENES_ <- PTHWGENES[grep("KEGG_PENTOSE_PHOSPHATE_PATHWAY",names(PTHWGENES),fixed=T)]
    ## PTHWGENES_<- PTHWGENES[grep("KEGG_VEGF_SIGNALING_PATHWAY",names(PTHWGENES),fixed=T)]
    g.list <- rownames(expTable.b)
    cell.names <- colnames(expTable.b)
    empty.ret.list <- list(PLACEHOLDER_PATHWAY=list(r=NA,p=NA,if.TF=NA,n.expd.t.TF=NA,n.t.TF=NA,cell.id=NA))
    enrich.str <- parallel::mclapply(as.list(1:dim(expTable.b)[2]),function(cell){
        cell.exp <- expTable.b[,cell] ## this has all gene names set already
        R.expd <- g.list[cell.exp & (g.list %in% RLIST)] # expressed R's
        TF.expd <- g.list[cell.exp & (g.list %in% TFLIST)] # expressed TF's
        ## names(cell.exp) <- g.list #important to properly subset cell.exp array
        ## cell.exp - bool array, named
        ## g.list is names for cell.exp
        ## g.list[cell.exp] -- genes expressed
        ## cell.exp["GENE"] -- is "GENE" expressed?
        ##cat("#R expd:",length(R.expd),head(R.expd),"#TF expd:",length(TF.expd),head(TF.expd),"\n")
        r.lst.pthw <- lapply(PTHWGENES_,function(pthw.genes){
            ## Some PTHWGENES sets have duplicated entries, e.g. KEGG_VEGF_SIGNALING_PATHWAY
            ## This affects the p-values, so remove duplicates:
            pthw.genes <- unique(pthw.genes)
            if(any(R.expd %in% pthw.genes) && any(TF.expd %in% pthw.genes)){
                pthw.genes.idx <- g.list %in% pthw.genes #excludes NA: cf. pthw.genes %in% g.list
                non.pthw.genes.idx <- !(g.list %in% pthw.genes)
                ## Cell-by-cell enrichment
                ## NA's should be appropriately accounted for, when gene is not found in data
                ## 1. cell.exp[pthw.genes.idx] -- expd genes of pthw, excludes NA's
                ## 2. !cell.exp[pthw.genes] -- pthw.gene not-expr including NA's
                ## convert NA to TRUE then use => non-present genes considered not expressed,
                ## it is more conservative position as with the bulk pipeline (see RTF.pairing.commap)
                ## 3. cell.exp[non.pthw.genes.idx] -- non.pthw.gene expr, excluding NA's
                ## 4. !cell.exp[non.pthw.genes.idx] -- non-pthw.gene not-expr excluding NA's,
                ## non.pthw.genes' NA's are not available in principle.
                ## ========
                ## Treat not-expressed in pthw genes separately, to substitute NA->TRUE
                not.exprd.in.pthw <- !cell.exp[pthw.genes]
                not.exprd.in.pthw[is.na(not.exprd.in.pthw)] <- TRUE
                pval <- fisher.test(matrix(c(sum(cell.exp[pthw.genes.idx]), # expd in pthw
                                          sum(not.exprd.in.pthw), # not-expd in pthw
                                          sum(cell.exp[non.pthw.genes.idx]), # expd in outside of pthw
                                          sum(!cell.exp[non.pthw.genes.idx])), # not-expd in outside of pthw
                                          ncol=2),
                                    alternative="greater")$p.value
                if(pval < 0.05){## not adjusted here! Pre-selection on enrichment
                    pthw.expd.R <- R.expd[R.expd %in% pthw.genes] ## R's of the pthw
                    pthw.expd.TF <- TF.expd[TF.expd %in% pthw.genes] ## interface TFs!
                    if(length(pthw.expd.TF) == 0){#if yes, no need to proceed
                        return(NA)
                    }
                    
                    ## +++++++++++++++++++++++++++++++++++++++++++++++
                    ## Check targets for compatibility and select iTFs
                    targets <- sapply(pthw.expd.TF,function(x)
                        as.character(TFTF$Target[TFTF$Source == x]),simplify=FALSE)
                    effects <- sapply(pthw.expd.TF,function(x)
                        as.character(TFTF$Effect[TFTF$Source == x]),simplify=FALSE)
                    ## ++++++++++++++++++++ 
                    ## As a result of the meeting(11.10.2019) we will
                    ## consider all compatible targets
                    n.pos.compat.t <- vector("integer",length(pthw.expd.TF)) # num + compat.t
                    n.neg.compat.t <- vector("integer",length(pthw.expd.TF)) # num - compat.t
                    n.all.compat.t <- vector("integer",length(pthw.expd.TF)) # all compat t
                    n.all.t <- vector("integer",length(pthw.expd.TF)) # num all t
                    compat.t <- vector("list",length(pthw.expd.TF)) # list of compat t
                    pos.compat.t <- vector("list",length(pthw.expd.TF)) # list of + compat t
                    for(i in 1:length(pthw.expd.TF)){
                        pos.compat.id <- (cell.exp[ targets[[i]] ]) & # expressed
                            (effects[[i]] == "1") & # target is to be activated
                            (!is.na(cell.exp[ targets[[i]] ])) # it is available in data
                        neg.compat.id <- (!cell.exp[ targets[[i]] ]) & # not expressed
                            (effects[[i]] == "-1") & # target is to be inhibited
                            (!is.na(cell.exp[ targets[[i]] ])) # available in data
                        n.pos.compat.t[i] <- sum( pos.compat.id )
                        n.neg.compat.t[i] <- sum( neg.compat.id )
                        n.all.compat.t[i] <- n.pos.compat.t[i] + n.neg.compat.t[i]
                        n.all.t[i] <- sum(!is.na(cell.exp[ targets[[i]] ])) # all available in data
                        compat.t[[i]] <- targets[[i]][ pos.compat.id | neg.compat.id ]
                        pos.compat.t[[i]] <- targets[[i]][ pos.compat.id ]
                    }
                    names(compat.t) <- pthw.expd.TF
                    names(pos.compat.t) <- pthw.expd.TF
                    ## +++++++++++++++++++++++++++++++
                    return(list(r=pthw.expd.R,p=pval,
                                if.TF=pthw.expd.TF,n.compat.t=n.all.compat.t,
                                n.pos.compat.t=n.pos.compat.t,
                                n.all.t=n.all.t,compat.t=compat.t,
                                pos.compat.t=pos.compat.t,cell.id=cell))
                }
                else
                    return(NA)
            }
            else
                return(NA)
        })
        ## Remove NA's: not significant pathways
        r.lst.pthw <- r.lst.pthw[sapply(r.lst.pthw,is.list)]
        ## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ## Adjust the pthw p-value and SELECT pathways: ENRICHMENT SELECTION
        p.adj <- p.adjust(unlist(lapply(r.lst.pthw,function(x){
            x$p
        })),"fdr")
        ## Select based on enrichment p-value
        idx <- p.adj < 0.05
        r.lst.pthw <- r.lst.pthw[idx]
        p.adj <- p.adj[idx]
        ## Keep the names
        n <- names(r.lst.pthw)
        ## Report adjusted p-values
        r.lst.pthw <- lapply(seq_along(r.lst.pthw),function(i){
            r.lst.pthw[[i]]$p.adj <- unname(p.adj[i])
            r.lst.pthw[[i]]
        })
        names(r.lst.pthw) <- n
        if(length(r.lst.pthw) == 0){
            cat("Warning: not found any R/Pthw/TF pairing in a whole cell...\n")
            return(empty.ret.list)
        }
        else
            return(r.lst.pthw)
    })
    ## Cell names
    names(enrich.str) <- cell.names
    return(enrich.str)
}

## ++++++++++++++++++++++++++++
## Target compatibility
## ++++++++++++++++++++++++++++

compat.t.iTF <- function(iTF,cell.exp){
    ## Return compatible targets of an iTF for a single cell (gene-named array of exp)
    targets <- TFTF$Target[TFTF$Source == iTF]
    effects <- TFTF$Effect[TFTF$Source == iTF]
    idx <- ((cell.exp[ targets ]) & (!is.na(cell.exp[ targets ])) & (effects == "1")) |
        ((!cell.exp[ targets ]) & (!is.na(cell.exp[ targets ])) & (effects == "-1"))
    return(targets[ idx ])
}

compat.frac.iTF <- function(iTF,cell.exp){
    ## Return compatibility counts of the targets of an iTF
    targets <- TFTF$Target[TFTF$Source == iTF]
    effects <- TFTF$Effect[TFTF$Source == iTF]
    n.compat.t <- sum((cell.exp[ targets ]) & (!is.na(cell.exp[ targets ])) & (effects == "1")) +
        sum((!cell.exp[ targets ]) & (!is.na(cell.exp[ targets ])) & (effects == "-1"))
    n.all.t <- sum(!is.na(cell.exp[ targets ])) ## all available in data
    return(c(n.compat.t,n.all.t - n.compat.t))
}

compat.frac.t.iTF <- function(t,iTFs,cell.exp){
    ## Return compatibility counts for a paired array of iTFs and t's
    ## t and iTFs of the same size
    effects <- sapply(seq_along(t),function(i) TFTF$Effect[(TFTF$Source == iTFs[i])
                                                           & (TFTF$Target == t[i])])
    n.compat <- sum((cell.exp[t]) & (!is.na(cell.exp[t])) & (effects == "1")) +
        sum((!cell.exp[t]) & (!is.na(cell.exp[t])) & (effects == "-1"))
    return(c(n.compat,length(t) - n.compat))
}

compat.prob <- function(t,cell.exp){
    ## Compatibility probability estimation
    ## Simple compatibility probability estimation
    regs <- as.character(TFTF$Source[TFTF$Target == t])
    effs <- as.character(TFTF$Effect[TFTF$Target == t])
    all.regs <- regs[!is.na(cell.exp[regs])] # all regulators available in data
    all.exp.regs <- regs[(cell.exp[regs]) & (!is.na(cell.exp[regs]))] # all regulators available in data
    actrs.a <- regs[(effs == "1")  & (cell.exp[regs]) & (!is.na(cell.exp[regs]))] # working activators
    inhrs.a <- regs[(effs == "-1") & (cell.exp[regs]) & (!is.na(cell.exp[regs]))] # working inhibitors
    actrs.i <- regs[(effs == "1")  & (!cell.exp[regs]) & (!is.na(cell.exp[regs]))] # not working activators
    inhrs.i <- regs[(effs == "-1") & (!cell.exp[regs]) & (!is.na(cell.exp[regs]))] # not working inhibitors
    if(cell.exp[t]){## target is on
        return((length(actrs.a))/( length(all.exp.regs) ))
    }
    else if(!cell.exp[t]){## target is off
        return((length(inhrs.a))/( length(all.exp.regs) ))
    }
    else {## is.na, e.g.
        return(NA)
    }
}

compat.frac.back <- function(cell.exp,comp=1,t.set.size=10,n.perm=1000){
    ## Background of compatibility
    PTHWGENES_ <- PTHWGENES[grep("KEGG",names(PTHWGENES))]  # KEGG
    g.list <- names(cell.exp) # Gene names in data
    TF.expd <- g.list[cell.exp & (g.list %in% TFLIST)] # expressed TF's
    ## ================================================
    ## 1. Background could come from all pathways (iTFs)
    all.iTF <- unique(unlist(lapply(PTHWGENES_,function(pthw.genes) # all expressed iTF from data
        TF.expd[TF.expd %in% pthw.genes])))
    ##cat("# iTFs expressed:",length(all.iTF),"\n")
    ## 2. Background could come from all expressed TF's (iTF's and not)
    ## all.iTF <- TF.expd
    ## All targets 
    all.t <- unique(TFTF$Target[(TFTF$Target %in% g.list) & (TFTF$Source %in% all.iTF)])

    ## ================================================
    ## 1. Compatibility of iTF
    if(comp == 1){
        y <- c(0,0)
        for(it in seq_along(all.iTF)){
            y <- compat.frac.iTF(all.iTF[it],cell.exp) + y
        }
    }
    ## ================================================
    ## 2. Compatibility of t's
    if(comp == 2){
        y <- matrix(NA,nrow=n.perm,ncol=2)
        for(i in 1:n.perm){
            rnd.t <- sample(all.t,t.set.size,replace=TRUE)
            ## cat("rnd.t:",head(rnd.t),"\n")
            rnd.iTFs <- vector("character",length(rnd.t))
            for(j in 1:length(rnd.t)){
                ## if(!is.character(TFTF$Source[TFTF$Target == rnd.t[j]])){
                ##     cat("Non-character:",TFTF$Source[TFTF$Target == rnd.t[j]],"\n")
                ## }
                ## cat("rnd iTFs:",head(TFTF$Source[TFTF$Target == rnd.t[j]]),"\n")
                rnd.iTFs[j] <- sample(TFTF$Source[TFTF$Target == rnd.t[j]],1)
            }
            y[i,] <- compat.frac.t.iTF(rnd.t,rnd.iTFs,cell.exp)
        }
    }
    ## ==============================================
    ## 3. Probabilistic compatibility of t's
    if(comp == 3){
        y <- c(0,0)
        for(t in seq_along(all.t)){
            p <- compat.prob(all.t[t],cell.exp)
            print(p)
            if(runif(1) < p)
                y[1] <- y[1] + 1
            else
                y[2] <- y[2] + 1
        }
    }
    return(y)
}

compat.stat <- function(enrich.str,expTable.b,cell.id,comp=1,verb=0){
    ## A single population input is required! That is:
    ## enrich.str is for a single population!
    ## expTable.b is for a single population!
    ## Calculate significance of iTF-t's interactions
    cell.exp <- expTable.b[,cell.id]
    enrich.pthws.lst <- enrich.str[[cell.id]]
    ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## 1. Fisher test: compatibility vs. background iTF-effectors
    ## Background can be formed iTF-by-iTF (comp=1), i.e. by counting multiple
    ## times all targets for all iTFs or target-by-target (comp=3), i.e. assigning
    ## a target using a random number according to its prob to be compatible.
    if(comp == 1 || comp == 3){
        enrich.pthw.TF <- unique(unlist(lapply(enrich.pthws.lst,function(x) x$if.TF[x$n.all.t != 0])))
        pvals <- vector("numeric",length(enrich.pthw.TF))
        pvals.2 <- vector("numeric",length(enrich.pthw.TF))
        y <- compat.frac.back(cell.exp,comp=comp) ## background counts
        for(it in seq_along(enrich.pthw.TF)){
            x <- compat.frac.iTF(enrich.pthw.TF[it],cell.exp)
            pvals[it] <- 1.0 - phyper(x[1]-1,y[1],y[2],x[1]+x[2]) ## x[1]-1 to get P(X >= x[1])
            if(verb==1){
                comp.mat <- matrix(c(x[1], # compatible for the given iTF
                                     y[1]-x[1], # compatible in the background iTFs
                                     x[2], # not-compat for the given iTF
                                     y[2]-x[2]), # not-compat in the background iTFs
                                   ncol=2)
                cat("iTF",enrich.pthw.TF[it],":",comp.mat,"\n")
            }
            ## pvals.2[it] <- fisher.test(comp.mat,
            ##                            alternative="greater")$p.value

        }
        names(pvals) <- enrich.pthw.TF
        ## names(pvals.2) <- enrich.pthw.TF
        if(verb == 1)
            print(pvals[order(pvals)])
        ## cat("===\n")
        ## print(pvals.2[order(pvals.2)])
        return(names(pvals)[pvals < 0.05])
    }
    ## ++++++++++++++++++++++++++++++++++++++++++++
    ## 2. Permutation test: Take random set (same size as the given) of targets from
    ## the background and estim probability to see the same compat.frac or higher as in this
    ## null distribution.
    if(comp == 2){
        ## First collect all unique iTFs from the pathways to avoid multiple perm test on same iTFs
        ## Take if.TF with targets
        enrich.pthw.TF <- unique(unlist(lapply(enrich.pthws.lst,function(x) x$if.TF[x$n.all.t != 0])))
        pvals <- vector("numeric",length(enrich.pthw.TF))
        for(it in seq_along(enrich.pthw.TF)){
            x <- compat.frac.iTF(enrich.pthw.TF[it],cell.exp)
            y <- compat.frac.back(cell.exp,t.set.size=sum(x),n.perm=100,comp=2)
            ecdFn <- ecdf(y[,1])
            pvals[it] <- 1.0 - ecdFn(x[1])
            cat("*** iTF:",enrich.pthw.TF[it],":",c(x[1],x[2]),":",pvals[it],":")
            cat(" ",paste(names(enrich.pthws.lst)[
                         sapply(enrich.pthws.lst,function(x)
                             enrich.pthw.TF[it] %in% x$if.TF)],
                         collapse=", "),"\n")

        }
        names(pvals) <- enrich.pthw.TF
        cat(" === Significance:\n")
        print(pvals[order(pvals)])
    }
    ## ++++++++++++++++++++++++++++++++++++++++++++++++
    ## 3. Persistence of certain iTF-target interaction
    ## throughout the cells
    if(comp == 4){
        enrich.pthw.TF <- unique(unlist(lapply(enrich.pthws.lst,function(x) x$if.TF[x$n.all.t != 0])))
        interacs <- c()
        for(it in seq_along(enrich.pthw.TF)){
            interacs <- c(interacs,paste(enrich.pthw.TF[it],compat.t.iTF(enrich.pthw.TF[it],cell.exp),sep="_"))
            ## cat(it,"/",length(enrich.pthw.TF),":",setdiff(unique(unlist(lapply(enrich.pthws.lst,function(x) x$compat.t[[enrich.pthw.TF[it]]]))),
            ##                                               compat.t.iTF(enrich.pthw.TF[it],cell.exp)),"\n")
        }
        return(interacs)
    }
}

compat.interacs.pop.ui <- function(enrich.str,expTable.b,pop,comp=4,n.repeat=1){
    ## User-interface function to get compatible interacs for a pop
    ## Full enrich str and exp table are assumed
    erch.str <- enrich.str[[pop]]
    expTable.b.pop <- expTable.b[,colnames(expTable.b) == pop]
    interacs.tbl <- compat.interacs.pop(erch.str,expTable.b.pop,comp=comp,n.repeat=n.repeat)
    return(sort(interacs.tbl,decreasing=TRUE))
}

compat.interacs.pop <- function(enrich.str,expTable.b,comp=1,n.repeat=1){
    ## Produce iTFs/iTF-t's and print sorted significance thereof
    ## Population specific parameters are inputs here, that is:
    ## enrich.str is for a single population
    ## expTable.b is for a single population
    ## This function is not intended for standalone use.
    tfs <- c()
    for(i in 1:length(enrich.str)){
        if(comp == 3)
            for(r in 1:n.repeat)
                tfs <- c(tfs,compat.stat(enrich.str,expTable.b,i,comp=comp))
        else
            tfs <- c(tfs,compat.stat(enrich.str,expTable.b,i,comp=comp))
    }
    if(comp == 3)
        print(sort(table(tfs)/(n.repeat*length(enrich.str)),decreasing=TRUE))
    else if(comp == 4)
        tf.tbl <- sort(table(tfs)/length(enrich.str))
    else
        print(sort(table(tfs)/length(enrich.str),decreasing=TRUE))
    if(comp == 4)
        return(tf.tbl)
    else
        return(tfs)
}

compat.prune <- function(enrich.str,expTable.b,min.q=0.95,pop="all"){
    ## Prune based on the compatibility criteria
    ## Enrich. struct for all populations or a single pop is allowed (pop param)
    if(pop == "all"){
        pops <- names(enrich.str)
        for(p in pops){
            enrich.str[[p]] <- compat.prune.pop(enrich.str[[p]],expTable.b[,colnames(expTable.b) == p,drop=FALSE],
                                                min.q=min.q)
        }
    }
    else {
        if(length(unique(names(enrich.str))) != 1){
            stop("Does not look like a single population enrich-struct presented.\n")
        }
        enrich.str <- compat.prune.pop(enrich.str,expTable.b[,colnames(expTable.b) == pop,drop=FALSE],
                                       min.q=min.q)
    }
    return(enrich.str)
}

compat.prune.pop <- function(enrich.str.pop,expTable.b.pop,min.q=0.95,verb=1){
    ## Prune a single population enrichment struct
    tf.tbl <- compat.interacs.pop(enrich.str.pop,expTable.b.pop,comp=4)
    if(verb==1){
        print(sort(tf.tbl[tf.tbl > quantile(tf.tbl,probs=min.q)],decreasing=TRUE)[1:10])
    }
    sign.interacs <- names(tf.tbl)[tf.tbl > quantile(tf.tbl,probs=min.q)] # quantile based filter
    prune.rate <- vector("numeric",length(enrich.str.pop))
    cat(unique(names(enrich.str.pop)),":",sep="")
    for(i in 1:length(enrich.str.pop)){## over cells in pop
        mark.arr <- vector("logical",length(enrich.str.pop[[i]]))
        for(j in 1:length(enrich.str.pop[[i]])){## over pathways in a cell
            mark.arr[j] <- check.pthw(enrich.str.pop[[i]][[j]],sign.interacs)
        }
        ## print(names(enrich.str.pop[[i]])[mark.arr])
        enrich.str.pop[[i]] <- enrich.str.pop[[i]][mark.arr]
        prune.rate[i] <- sum(!mark.arr)/length(mark.arr)
        ## print(prune.rate[i])
    }
    cat(" avg pruning rate:",round(mean(prune.rate),1),"\n")
    return(enrich.str.pop)
}

check.pthw <- function(pthw,sign.iTF.t.interacs){
    ## Check a pathway for the iTF-t compatibility
    interacs <- c()
    for(itf in pthw$if.TF){
        interacs <- c(interacs,paste(itf,pthw$compat.t[[itf]],sep="_"))
    }
    if(any(interacs %in% sign.iTF.t.interacs)){
        return(TRUE)
    }
    else{
        return(FALSE)
    }
}
## ++++++++++++++++++++++++++++
## Analyze enrichment structure
## ++++++++++++++++++++++++++++

phenoTF.maintain.rec <- function(){
    
}

## ++++++++++++++++++++++++++++++++++++++++++
## Standard ways to analyze enrichment struct
## ++++++++++++++++++++++++++++++++++++++++++

## Atomic functions to inquire a cell (collection of pathway-structs) on:
## 1. iTF-tTF interactions
## 2. pathways
## 3. receptors
##
has.cell.compat.t <- function(cell,iTF,tTF,verb=1,any.logic=TRUE){
    ## any.logic=T: any combinations of iTF-tTF
    ## any.logic=F: all combinations of iTF-tTF
    if(is.null(iTF) || is.null(tTF))
        return(NULL)
    test.interacs <- paste(iTF,tTF)
    if(verb)
        print(test.interacs)
    interac.found <- rep(FALSE,length(test.interacs))
    for(ipth in seq_along(cell)){## over pathways in the cell
        pth <- cell[[ipth]]
        target.interacs <- unname(unlist(sapply(pth$if.TF,function(x) paste(x,pth$compat.t[[x]]))))
        interac.found <- ifelse(test.interacs %in% target.interacs,
                                test.interacs %in% target.interacs,
                                interac.found)
    }
    if((any.logic && any(interac.found)) || (!any.logic && all(interac.found)))
        return(TRUE)
    else
        return(FALSE)
}

has.cell.pthw <- function(cell,pthw,verb=1){
    if(is.null(pthw))
        return(NULL)
    pth.idx <- grep(toupper(pthw),names(cell),fixed=TRUE)
    if(verb)
        print(names(cell)[pth.idx])
    if(any(pth.idx))
        return(TRUE)
    else
        return(FALSE)
}

has.cell.rec <- function(cell,rec,any.logic=TRUE){
    ## any.logic=T: any of the rec's
    ## any.logic=F: all of the rec's
    if(is.null(rec))
        return(NULL)
    rec.found <- rep(FALSE,length(rec))
    for(ipth in seq_along(cell)){## over pathways in the cell
        pth <- cell[[ipth]]
        rec.found <- ifelse(rec %in% pth$r,rec %in% pth$r,rec.found)
    }
    if((any.logic && any(rec.found)) || (!any.logic && all(rec.found)))
        return(TRUE)
    else
        return(FALSE)
}

has.cell.multi <- function(cell,rec=NULL,pthw=NULL,iTF=NULL,tTF=NULL,any.logic=TRUE,verb=1){
    ## Multiple requests for a single cell
    ## NOTE: R Pthw and TF-TF may not be necessarily connected physically, i.e. belong
    ## to the same path-way. Just co-occurrence in a cell
    cell.id <- c()
    rec.b <- has.cell.rec(cell,rec,any.logic=any.logic)
    pthw.b <- has.cell.pthw(cell,pthw,verb=verb)## here pthw not necessarily connected to R/TF
    tf.b <- has.cell.compat.t(cell,iTF,tTF,verb=verb,any.logic=any.logic)
    cell.b <- c(rec.b,pthw.b,tf.b) ## Note, NULL's will not be added here naturally
    if(verb)
        cat("R Pthw TF-TF:",cell.b,"\n")
    if(all(cell.b))
        return(TRUE)
    else
        return(FALSE)
}

get.cell.id <- function(enrich.str,pop=NULL,rec=NULL,pthw=NULL,iTF=NULL,tTF=NULL){
    ## Cell id's of the cells with conditions (only non-NULL applies)
    ## AND logic is maintaned for multiple requests: all cells having rec, pthw, and TF-TF interac
    ## Also see notes in has.cell.multi() and other has.cell.? functions
    ## ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ## The output list, size = # pops in enrich.str
    out.cell.id <- vector("list",length(unique(names(enrich.str))))
    names(out.cell.id) <- unique(names(enrich.str))
    if(length(out.cell.id) == 1){ ## one pop case, ignore pop=non-NULL
        out.cell.id[[1]] <- rep(FALSE,length(enrich.str))
        for(ic in 1:length(enrich.str)){## over cells
            cell <- enrich.str[[ic]] ## cell pthw-struct
            if(has.cell.multi(cell,rec=rec,pthw=pthw,iTF=iTF,tTF=tTF,verb=0))
                out.cell.id[[1]][ic] <- TRUE
        }
    }
    else { ## multiple pop case
        for(ip in 1:length(enrich.str)){## over pops
            out.cell.id[[ip]] <- rep(FALSE,length(enrich.str[[ip]]))
            if((!is.null(pop) && any(pop == names(enrich.str)[ip])) |
               (is.null(pop))){## either all or selected pops
                for(ic in 1:length(enrich.str[[ip]])){## over cells
                    cell <- enrich.str[[ip]][[ic]]
                    if(has.cell.multi(cell,rec=rec,pthw=pthw,iTF=iTF,tTF=tTF,verb=0))
                    out.cell.id[[ip]][ic] <- TRUE
                }
            }
        }
    }
    return(lapply(out.cell.id,which))## to return actual indices
}

get.pthw.params <- function(enrich.str,pthw.name,pop,rec.name=NULL){
    ## NOTE: enrich.str is a list of enrich. structures for each population!
    ## This list is supposed to be saved during sc.commap calls, just for convenience
    ## The other functions work on individual single pop data structures
    r.lst.ext <- enrich.str[[pop]]
    n.cells <- length(r.lst.ext)
    out <- list(n.cells=n.cells,
                params=data.frame(p.adj=rep(NA,times=n.cells),if.TF=NA,
                                  n.compat.t=NA,n.all.t=NA,cell.id=NA))
    for(c in 1:length(r.lst.ext)){
        x <- r.lst.ext[[c]] # a cell's info (pthw etc.)
        target.pthw.idx <- grep(toupper(pthw.name),names(x),fixed=TRUE)
        if(length(target.pthw.idx)){## Found the pathway
            if((!is.null(rec.name)) && (rec.name %in% x[[target.pthw.idx]]$r)){
                first.na <- which(is.na(out$params$p.adj))[1]
                out$params$p.adj[first.na] <- x[[target.pthw.idx]]$p.adj
                out$params$if.TF[first.na] <- paste(x[[target.pthw.idx]]$if.TF,collapse=",")
                out$params$n.compat.t[first.na] <- paste(x[[target.pthw.idx]]$n.compat.t,collapse=",")
                out$params$n.all.t[first.na] <- paste(x[[target.pthw.idx]]$n.all.t,collapse=",")
                out$params$cell.id[first.na] <- paste(x[[target.pthw.idx]]$cell.id,collapse=",")
            }
            if(is.null(rec.name)){
                first.na <- which(is.na(out$params$p.adj))[1]
                out$params$p.adj[first.na] <- x[[target.pthw.idx]]$p.adj
                out$params$if.TF[first.na] <- paste(x[[target.pthw.idx]]$if.TF,collapse=",")
                out$params$n.compat.t[first.na] <- paste(x[[target.pthw.idx]]$n.compat.t,collapse=",")
                out$params$n.all.t[first.na] <- paste(x[[target.pthw.idx]]$n.all.t,collapse=",")
                out$params$cell.id[first.na] <- paste(x[[target.pthw.idx]]$cell.id,collapse=",")
            }
        }
    }
    out$params <- out$params[!is.na(out$params[,1]),]
    if(dim(out$params)[1] == 0) return(out)
    rownames(out$params) <- seq(1,dim(out$params)[1])
    return(out)
}

get.pthw.stat <- function(enrich.str,pop,rec.name=NULL){
    ## NOTE: enrich.str is a list of enrich. structures for each population!
    ## This list is supposed to be saved during sc.commap calls, just for convenience
    ## The other functions work on individual single pop data structures
    out <- data.frame(Pthw=rep(NA,times=length(PTHWGENES)),Cell.frac=NA)
    r.lst.ext <- enrich.str[[pop]]
    n.cells <- length(r.lst.ext)
    ## cat("Total num cells:",n.cells,"\n")
    lapply(r.lst.ext,function(x) { ## per cell
        pthw.n <- names(x)
        lapply(seq_along(x), function(i) { ## per pthw in a cell
            if((!is.null(rec.name)) && (rec.name %in% x[[i]]$r)){
                idx <- which(out$Pthw == pthw.n[i])
                if(length(idx))## pthw found
                    out$Cell.frac[idx] <<- out$Cell.frac[idx] + 1
                else{## pthw not found in the table
                    first.na <- which(is.na(out$Pthw))[1]
                    out$Pthw[first.na] <<- pthw.n[i]
                    out$Cell.frac[first.na] <<- 1
                }
            }
            if(is.null(rec.name)){
                idx <- which(out$Pthw == pthw.n[i])
                if(length(idx))
                    out$Cell.frac[idx] <<- out$Cell.frac[idx] + 1
                else{
                    first.na <- which(is.na(out$Pthw))[1]
                    out$Pthw[first.na] <<- pthw.n[i]
                    out$Cell.frac[first.na] <<- 1
                }
            }
        })
    })
    out <- out[!is.na(out[,1]),1:2]
    if(dim(out)[1] == 0) return(out)
    out$Cell.frac <- out$Cell.frac / (1.0*n.cells)
    out <- out[order(out$Cell.frac,decreasing=TRUE),1:2]
    rownames(out) <- seq(1,dim(out)[1])
    return(out)
}

########################
### Phenotypic TF tables
########################

get.all.pheno.TF.tables <- function(cm,enrich.str,thresPath=0.0,min.n.expd.tTF=1,min.tTF.frac=0.0){
    r.pop.pairs <- get.unique.rec.pop(cm)
    pops <- unique(cm$Rec.pop)
    ph.tbl.all <- vector("list",length(r.pop.pairs))
    for(i in 1:length(r.pop.pairs)){
        ph.tbl <- get.pheno.TF.table(enrich.str, r.pop.pairs[[i]][2], r.pop.pairs[[i]][1],
                                     thresPath=thresPath,min.n.expd.tTF=min.n.expd.tTF,
                                     min.tTF.frac=min.tTF.frac)
        ph.tbl.all[[i]] <- ph.tbl
    }
    names(ph.tbl.all) <- sapply(r.pop.pairs,paste,collapse="_")
    return(ph.tbl.all)
}

get.pheno.TF.table <- function(enrich.str, pop, recptr, thresPath = 0.0, min.n.expd.tTF = 1, min.tTF.frac = 0.0){
    ## Produce a phenotype table based on the if.TF and t.TF information
    pthws <- get.pthw.stat(enrich.str,pop,recptr)
    kept.pthws <- pthws[pthws$Cell.frac > thresPath,1]
    n.cells <- length(enrich.str[[pop]])
    ## Merged params table from all pathways
    params.merged <- c()
    for(ptw in kept.pthws)
        params.merged <- rbind(params.merged,get.pthw.params(enrich.str,ptw,pop,recptr)$params)
    ## print(params.merged)
    ## Extract values from the merged table
    ## cat(pop,"-",recptr,"\n")
    all.iTF <- strsplit(params.merged$if.TF,",")
    n.expd.tTF <- strsplit(params.merged$n.expd.t.TF,",")
    n.all.tTF <- strsplit(params.merged$n.t.TF,",")
    ## print(head(all.iTF))
    ## print(head(n.expd.tTF))
    ## print(n.all.tTF)
    u.iTF <- unique(unlist(all.iTF))
    ## print(u.iTF)
    n <- length(u.iTF)
    ## Phenotype table
    phenoTF <- list(nc=n.cells,tbl=data.frame(if.TF=rep(NA,n),cell.id=rep(NA,n),cell.frac=rep(NA,n),
                                              t.TF.frac=rep(NA,n),n.expd.t.TF=rep(NA,n)))
    no.match <- c()
    lapply(seq_along(u.iTF),function(j){
        phenoTF$tbl$if.TF[j] <<- u.iTF[j] ## if.TF
        ## idx <- sapply(all.iTF,function(x) which(u.iTF[j] == x)[1]) ## row & position in params.merged
        ## row (if !is.na) & position (value) in params.merged
        idx <- sapply(seq_along(all.iTF),function(i) which((u.iTF[j] == all.iTF[[i]]) &
                                                           (as.numeric(n.expd.tTF[[i]]) >= min.n.expd.tTF) &
                                                           (as.numeric(n.expd.tTF[[i]]) / as.numeric(n.all.tTF[[i]]) >
                                                            (min.tTF.frac-1e-6)))[1])
        ##print(idx)
        if(all(is.na(idx))){## No matched iTF found
            no.match <<- c(no.match,j)
        }
        c.idx <- params.merged$cell.id[!is.na(idx)] ## cell indices where the if.TF is found
        phenoTF$tbl$cell.id[j] <<- paste(unique(c.idx),collapse=",") ## take unique cell indices
        phenoTF$tbl$cell.frac[j] <<- length(strsplit(phenoTF$tbl$cell.id[j],",")[[1]])/(1.0*n.cells)
        phenoTF$tbl$t.TF.frac[j] <<- paste(sapply(which(!is.na(idx)),function(k)
            round(as.numeric(n.expd.tTF[[k]][idx[k]]) / as.numeric(n.all.tTF[[k]][idx[k]]),2)),collapse=",")
        phenoTF$tbl$n.expd.t.TF[j] <<- paste(sapply(which(!is.na(idx)),function(k) as.numeric(n.expd.tTF[[k]][idx[k]])),collapse=",")
    })
    if(length(no.match) > 0){
        phenoTF$tbl <- phenoTF$tbl[-no.match,]
        ## print(dim(phenoTF$tbl))
        if(dim(phenoTF$tbl)[1] > 0)
            rownames(phenoTF$tbl) <- 1:dim(phenoTF$tbl)[1]
    }
    return(phenoTF)
}


