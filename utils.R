#########
### Utils
#########

## ==========================================
## Helper functions comparing 2 receptor sets
## ==========================================
get.rec.set.diff <- function(enrich.str,enrich.str.2,conserv.thr){
    r.lst <- get.conserved.receptors(enrich.str,conserv.thr)
    r.lst.pruned <- get.conserved.receptors(enrich.str.2,conserv.thr)
    x <- lapply(seq_along(r.lst),function(i) setdiff(r.lst[[i]],r.lst.pruned[[i]]))
    names(x) <- names(r.lst)
    return(x)
}

get.rec.set.same <- function(enrich.str,enrich.str.2,conserv.thr){
    r.lst <- get.conserved.receptors(enrich.str,conserv.thr)
    r.lst.pruned <- get.conserved.receptors(enrich.str.2,conserv.thr)
    x <- lapply(seq_along(r.lst),function(i) intersect(r.lst[[i]],r.lst.pruned[[i]]))
    names(x) <- names(r.lst)
    return(x)
}


## =================================
## Meta-analysis on the commap table
## =================================
get.unique.rec.pop <- function(cm){
    ## cm -- comMap table
    rec.pop.pairs <- unique(paste(cm$Receptor,cm$Rec.pop,sep="_"))
    spl <- strsplit(rec.pop.pairs,"_")
    return(spl)
}

rec.list <- function(ph.tbl.like.list){
    return(unlist(lapply(strsplit(names(ph.tbl.like.list),"_"),`[`,1)))
}

pop.list <- function(ph.tbl.like.list){
    return(unlist(lapply(strsplit(names(ph.tbl.like.list),"_"),`[`,2)))
}

rec.list.ct <- function(ph.tbl.like.list,ct){
    return(rec.list(ph.tbl.like.list)[pop.list(ph.tbl.like.list) == ct])
}

pop.list.rec <- function(ph.tbl.like.list,recs){
    ## Cell populations where all recs are present
    recs <- as.list(recs)
    pop.lst <- lapply(recs,function(x) pop.list(ph.tbl.like.list)[rec.list(ph.tbl.like.list) == x])
    return(Reduce(intersect,pop.lst))
}

num.arr <- function(x,sep=","){
    ## x -- sep separated character array of values, e.g. c("2,3,4,5","0.2,0.4,0.9,0.0,0.1")
    lapply(strsplit(x,sep),as.numeric)
}

char.arr <- function(x,sep=","){
    lapply(strsplit(x,sep),as.character)
}

two.close.factor <- function(n){
    ## Find tow-factor factorization with approximate equality
    x <- ceiling(n^(0.5))
    while(n %% x != 0)
        x <- x - 1
    if(x == 1){ ## Prime number
        return(two.close.factor(n+1))
    }
    return(c(x,n/x))
}

##############################################################################################
## Tests for getting similar if.TF sets with two separate thresholding: R.frac and if.TF.cfrac
rec.set.by.R.frac <- function(cm,ct,rf.thr){
    x <- get.unique.rec.pop(cm[cm$R.frac > rf.thr,])
    y <- sapply(x,function(y) {if(y[2] == ct){ return(y[1])} else {return(NA)}})
    return(y[!is.na(y)])
}

rec.set.by.iTF.cfrac <- function(ph.tbls,ct,cf.thr=0.2){
    if.TFs <- lapply(ph.tbls,function(x) x$tbl$if.TF[x$tbl$cell.frac > cf.thr])
    ##c.fracs <- lapply(ph.tbls,function(x) x$tbl$cell.frac[x$tbl$cell.frac > cf.thr])
    c.types <- pop.list(ph.tbls)
    ##u.c.types <- unique(c.types)
    ct.iTFs <- if.TFs[c.types == ct]
    ##ct.cFracs <- c.fracs[c.types == ct]
    ## How many Rs left after applying if.TF.cfrac cutoff
    flt.ct.iTFs <- ct.iTFs[sapply(ct.iTFs,length) != 0]
    ##flt.ct.cFracs <- ct.cFracs[sapply(ct.iTFs,length) != 0]
    ##return(unlist(lapply(strsplit(names(flt.ct.iTFs),"_"),function(x) x[1])))
    return(rec.list(flt.ct.iTFs))
}

## Set distance
set.dist <- function(s1,s2){
    max.set.l <- length(union(s1,s2))
    if(max.set.l == 0)
        return(1.0)
    inter.set.l <- length(intersect(s1,s2))
    ## Dist 1
    return( (max.set.l - inter.set.l) / (1.0*max.set.l) );
    ## Dist 2 (becomes infinity and hclust does not like it)
    ##return( max.set.l / (1.0*inter.set.l) );
}

multi.set.dist <- function(set.lst){
    dist.mat <- matrix(NA,nrow=length(set.lst),ncol=length(set.lst))
    for(i in 1:length(set.lst)){
        for(j in 1:length(set.lst)){
            dist.mat[i,j] <- set.dist(set.lst[[i]],set.lst[[j]])
        }
    }
    rownames(dist.mat) <- rec.list(set.lst)
    colnames(dist.mat) <- rec.list(set.lst)
    return(as.dist(dist.mat,diag=T))
}

## The test function
rec.set.dist <- function(ph.tbls,cm,ct,cf.thr,rf.thr){
    out1 <- rec.set.by.iTF.cfrac(ph.tbls,ct,cf.thr)
    out2 <- rec.set.by.R.frac(cm,ct,rf.thr)
    return(set.dist(out1,out2))
}

rec.set.dist.image <- function(ph.tbls,cm,ct){
    thr.seq <- seq(0.0,1.0,by=0.05)
    xy <- expand.grid(1:length(thr.seq),1:length(thr.seq))
    dist.mat <- matrix(NA,nrow=length(thr.seq),ncol=length(thr.seq))
    for(i in 1:dim(xy)[1]){
        dist.mat[xy[i,1],xy[i,2]] <- rec.set.dist(ph.tbls,cm,ct,thr.seq[xy[i,1]],thr.seq[xy[i,2]])
    }
    rownames(dist.mat) <- as.character(thr.seq)
    colnames(dist.mat) <- as.character(thr.seq)
    image(t(dist.mat),ylab="iTF.cell.frac",xlab="R.frac",main=ct)
    lines(c(0.0,1.0),c(0.0,1.0),col='black',lwd=4)
}

multi.rec.set.dist.image <- function(ph.tbls,cm){
    cts <- unique(pop.list(ph.tbls))
    par(mfrow=two.close.factor(length(cts)))
    for(ct in cts){
        cat(ct,"...\n")
        rec.set.dist.image(ph.tbls,cm,ct)
    }
}

## =====================
## Cell-by-cell routines
## =====================

### Booleanization
make.exp.bool.bm <- function(expTable.bm,exp.thr=1){
    ## Booleanization for filematrix objects
    expTable.b <- matrix(FALSE,nrow=dim(expTable.bm)[1],ncol=dim(expTable.bm)[2])
    for(i in 1:dim(expTable.bm)[2]){
        expTable.b[,i] <- expTable.bm[,i] > exp.thr
    }
    colnames(expTable.b) <- colnames(expTable.bm)
    rownames(expTable.b) <- rownames(expTable.bm)
    return(expTable.b)
}

make.exp.bool <- function(expTable,exp.thr=1){
    return(expTable > exp.thr)
}

select.pop.exp <- function(fullExpTable,p){
    popExpTable <- fullExpTable[,colnames(fullExpTable) == p,drop=FALSE]
    colnames(popExpTable) <- rep(p,dim(popExpTable)[2])
    return(popExpTable)
}

### This function is used within sc.commap and is applicable to both extended R.list and L.list ###
get.all.rl.tab <- function(rl.lst.ext,pop,mode="r"){
    ## NOTE: this function still uses pop parameter, to select pop if all are included (as in L case)
    ## See examples in sc.commap for R.list.ext.pop and L.list.ext.all separately
    rl.lst <- get.rl.list(rl.lst.ext,mode=mode)
    idx <- names(rl.lst) == pop ## select only cells from pop
    n.cells <- sum(idx)
    return( sort(table(unlist(rl.lst[idx])) / n.cells,decreasing=TRUE) )
}

get.rl.list <- function(lst.ext,mode="r"){
    ## Remove Pathway and other(p.vals etc.) information
    if(mode == "r"){
        lst <- lapply(lst.ext,function(x){ ## by cells
            ## Unique R's in each cell
            unique(unlist(lapply(x,function(y) y$r),use.names=FALSE)) ## by pathways
        })
        lst <- lst[sapply(lst,function(x) length(x) != 0)] ## remove empty pathways (after pruning)
        lst <- lst[!sapply(lst,function(x) all(is.na(x)))] ## remove PLACEHOLDER_PATHWAY NA's
    }
    else if(mode == "l"){
        lst <- lapply(lst.ext,function(x){
            ## For L.lists one less recursive list, as no pathway info
            unique(x$l)
        })
    }
    else{
        stop("Mode is unknown: 'r' or 'l' ?\n")
    }
    ## Returned is still a list or R/L in each cell
    return(lst)
}

### L-struct 
sc.L.list.ext.bm <- function(expTable.b.bm){
    g.list <- rownames(expTable.b.bm)
    L.list.ext <- vector("list",length=dim(expTable.b.bm)[2])
    names(L.list.ext) <- colnames(expTable.b.bm)
    for(i in 1:dim(expTable.b.bm)[2]){
        L.list.ext[[i]] <- list(l=g.list[expTable.b.bm[,i] & (g.list %in% LLIST)])
    }
    return(L.list.ext)
}

sc.L.list.ext <- function(expTable.b){
### Get the list of expressed Ligands, cell-by-cell
    g.list <- rownames(expTable.b)
    L.list.ext <- apply(expTable.b,2,function(cell.exp){
        L.expd <- g.list[cell.exp & (g.list %in% LLIST)] # expressed L's
        ## Could be more stuff added into a returning list, e.g., parameters
        ## ...
        return(list(l=L.expd))
    })
    return(L.list.ext)
}


################################
### Prepare data from the format
################################

prepare.data <- function(expr,ph,placeholder=NULL,
                         ph.cell.id.col="cell.id",ph.cell.type.col="cell.type"){
    ## Prepares the data of a certain format (see below) to be supplied to the commap.
    ## Format:
    ## rownames(expr) have a format ENSMBLID_GENESYMBOL
    ## colnames(expr) are cell id's
    ## phenotype (ph) data have cell id column (ph.cell.id.col)
    ##  and cell type column (ph.cell.type.col)
    ## Remove ENSMBLID part
    ## print(grep("HLA",rownames(expr),value=T))
    ## =======================================================
    ## If data is formatted without ENSMBLID field, then we take on gene symbol:
    ## check for "_" (underscore) in 10 randomly selected names.
    if(all(grepl("_",sample(rownames(expr),10))))
        new.rownames <- sapply(strsplit(rownames(expr),"_"),function(x) x[2])
    else
        new.rownames <- rownames(expr)

    if(!is.null(placeholder)){
        ## Just remove rows with un-identified gene symbols marked with the placeholder
        idx.to.remove <- new.rownames == placeholder
        expr <- expr[!idx.to.remove,]
        rownames(expr) <- new.rownames[!idx.to.remove]
        ## it will have a warning here automatically about the duplicated rownames
    }
    else {
        ## Just remove duplicates in the rownames, leaving potentially a single unknown placeholder
        idx.to.remove <- duplicated(new.rownames)
        cat(length(which(idx.to.remove)),"duplicated entries\n")
        expr <- expr[!idx.to.remove,]
        rownames(expr) <- new.rownames[!idx.to.remove]
    }
    ## Substitue cell.id with cell.type information
    if(!(ph.cell.id.col %in% colnames(ph)) || !(ph.cell.type.col %in% colnames(ph)))
        stop("Indicate the proper names of cell.id and cell.type columns")
    new.colnames <- sapply(colnames(expr),function(x) {
        make.names(as.character(ph[[ph.cell.type.col]][x == as.character(ph[[ph.cell.id.col]])]),
                   unique=FALSE)
    })
    colnames(expr) <- new.colnames
    return(expr)
}
