source("utils.R")
source("enrich.R")
source("synergy.R")
source("lr_pairing.R")

## Get Ligand-Receptor data
LRDAT0 <- list(HSA=read.delim("LR.pairs.HSA.iRef.txt",stringsAsFactors=FALSE),
               MMU=read.delim("LR.pairs.MMU.iRef.txt",stringsAsFactors=FALSE),
               HSAold=read.delim("LR.pairs.HSA.txt",stringsAsFactors=FALSE),
               MMUold=read.delim("LR.pairs.MMU.txt",stringsAsFactors=FALSE))

## Get TF data
TFDAT0 <- list(HSA=read.delim("Homo_sapiens_TF.txt",stringsAsFactors=FALSE),
               MMU=read.delim("Mus_musculus_TF.txt",stringsAsFactors=FALSE))
## Get pathway data
PTHWDAT0 <- list(HSA=strsplit(readLines("c2.cp.v6.2.symbols.gmt"),split='\t'),
                 MMU=strsplit(readLines("c2.cp.v6.2.symbols_mmu_ortho.gmt"),split='\t'))
PTHWGENES0 <- list(HSA=lapply(PTHWDAT0[["HSA"]],function(x) x[seq(3,length(x))]),
                   MMU=lapply(PTHWDAT0[["MMU"]],function(x) x[seq(3,length(x))]))
PTHWNAMES0 <- list(HSA=unlist(lapply(PTHWDAT0[["HSA"]],function(x) x[1])),
                   MMU=unlist(lapply(PTHWDAT0[["MMU"]],function(x) x[1])))

#############################################################
### The step-by-step instruction to running the pipeline:
### 1.  Run enrichment (run.enrichment)
### 2.  Prune the enrichment structure based on the
###     TF-compatibility criteria (TF.compatibility.prune)
### 3a. Match the characteristic/conserved receptors with
###     ligand-signals from elsewhere using mean-field
###     approximation (mean.field.match.lr)
### 3b. Get the characteristic/conserved receptors of the
###     populations (get.conserved.receptors).
#############################################################

## Produce an enrichment structure
run.enrichment <- function(fullExpTable,abs.exp.thr=0.0,
                           org="MMU",LR.db="new",svf.key=NULL,
                           large=FALSE,large.sep='\t'){
    ## Run the pathway enrichment and return an enrichment structure
    ## *** Input:
    ## fullExpTable -- a data matrix (row-genes,column-cells) as prepared by
    ## utils.R/prepare.data() or data file name (if large=T)
    ## abs.exp.thr -- an absolute expression cutoff for booleanization > abs.exp.thr => TRUE
    ## org -- abbreviation for the organism: "MMU" or "HSA"
    ## LR.db -- "old" (Ramilowski et al.) or "new" (iRef index)
    ## svf.key -- if not NULL, will save the enrichment struc to disk with this identifier
    ## large -- TRUE/FALSE, if large data processing is requested (using bigmemory package)
    ## large.sep -- separator for the large data file
    ## *** Output:
    ## enrich.struc -- an enrichment structure needed for the downstream analysis (see also
    ## enrich.R/sc.RTF.pairing())
    #######################################################################################
    ## Environments
    set.envir(org,LR.db=LR.db)
    ## Thresholds: for pathways and for L's
    cat("Abs.exp.thr:",abs.exp.thr,"\n")
    
    if(large){
        suppressPackageStartupMessages({
            require(bigmemory)
        })
        cat("Large data reading ")
        cat("(",fullExpTable,")...")
        binFile <- "expTable.bin"
        descFile <- "expTable.desc"
        fullExpTable <- read.big.matrix(fullExpTable,sep=large.sep,header=TRUE,
                                        has.row.names=TRUE,type="double",backingfile=binFile,
                                        descriptorfile=descFile)
        cat(" Done.\n")
    }
    ## Populations
    pops <- levels(as.factor(colnames(fullExpTable)))
    if(large){
        fullExpTable.b <- make.exp.bool.bm(fullExpTable,exp.thr=abs.exp.thr)
        cat("Original class:",class(fullExpTable),"; bool class:",class(fullExpTable.b),"\n")
    }
    else
        fullExpTable.b <- make.exp.bool(fullExpTable,exp.thr=abs.exp.thr)
    enrich.str <- vector("list",length=length(pops))
    names(enrich.str) <- as.character(pops)
    for(rp in pops){
        popExpTable.b <- select.pop.exp(fullExpTable.b,rp)
        cat("=== Population:",rp,",",dim(popExpTable.b)[2],"cells\n")
        st <- system.time(enrich.str[[rp]] <- sc.RTF.pairing(popExpTable.b))
        ## ---------------------------------
        cat("* Enrichment time:\n")
        print(st)
    }
    if(!is.null(svf.key)){
        tm <- as.character(Sys.time())
        dt <- unlist(strsplit(tm," "))
        fn <- paste("enrich_str",svf.key,dt[1],paste(unlist(strsplit(dt[2],":"))[1:2],collapse="."),
                    sep="_")
        cat("Saving enrich structure to disk (",paste0(fn,".RData"),")...\n",sep="")
        cat("Size:",object.size(enrich.str)/1000000.0,"MB\n")
        ## Save with version=2 for R.ver 3.4 and lower compatibility
        save(enrich.str,file=paste0(fn,".RData"),version=2)
    }
    return(enrich.str)
}

## Filter pathways based on the TF-compatibility criteria
TF.compatibility.prune <- function(fullExpTable,enrich.str,min.q=0.95,abs.exp.thr=0.0,large=FALSE,large.sep='\t'){
    ## Prune the enrichment structure based on the TF compatibility principles
    ## *** Input:
    ## fullExpTable -- data matrix (row-genes,column-cells, as returned by utils.R/prepare.data())
    ## or data file name (if large=T)
    ## enrich.str -- an enrichment struc returned by run.enrichment()
    ## abs.exp.thr -- an absolute expression cutoff for booleanization > abs.exp.thr => TRUE
    ## large -- TRUE/FALSE, if large data processing is requested (using bigmemory package)
    ## large.sep -- separator for the large data file
    ## *** Output:
    ## enrich.str.pruned -- a pruned enrichment str, that is with some pathways removed due to
    ## not satisfied TF compatibility criteria (see enrich.R/compat.prune())
    if(large){
        suppressPackageStartupMessages({
            require(bigmemory)
        })
        cat("Large data reading ")
        cat("(",fullExpTable,")...")
        binFile <- "expTable.bin"
        descFile <- "expTable.desc"
        fullExpTable <- read.big.matrix(fullExpTable,sep=large.sep,header=TRUE,
                                        has.row.names=TRUE,type="double",backingfile=binFile,
                                        descriptorfile=descFile)
        cat(" Done.\n")
        fullExpTable.b <- make.exp.bool.bm(fullExpTable,exp.thr=abs.exp.thr)
        cat("Original class:",class(fullExpTable),"; bool class:",class(fullExpTable.b),"\n")
    }
    else
        fullExpTable.b <- make.exp.bool(fullExpTable,exp.thr=abs.exp.thr)
    enrich.str.pruned <- compat.prune(enrich.str,fullExpTable.b,min.q=min.q,pop="all")
    return(enrich.str.pruned)
}

## Get characteristic receptors from an enrichment structure
get.conserved.receptors <- function(enrich.str,conserv.thr,uniprot=FALSE){
    ## Get lists of conserved receptors for each population
    ## *** Input:
    ## enrich.str -- an enrichment struc
    ## conserv.thr -- a relative (fraction) cutoff for the receptors: if expression(rec) > conserv.thr, accept rec
    ## uniprot -- TRUE/FALSE, whether one checks "membrane" association of the receptor in UniProt DB
    ## *** Output:
    ## r.list -- a named (by pop) list of conserved receptors for each population
    ################################################################################
    r.list <- vector("list",length(enrich.str))
    names(r.list) <- names(enrich.str)
    for(rp in names(enrich.str)){
        R.tab <- get.all.rl.tab(enrich.str[[rp]],rp,mode="r")
        R.tab <- R.tab[R.tab > conserv.thr]
        r.list[[rp]] <- sort(names(R.tab))
    }
    if(uniprot){## Filter more based on the "membrane" annotation in UniProt
        cat("Web-based annotation checking requested: may take some time...\n")
        source("uniprot.R")
        r.list <- lapply(r.list,uniprot.filter.receptors)
    }
    return(r.list)
}

## Get L-R matching using mean-field approximation
mean.field.match.lr <- function(fullExpTable,enrich.str,conserv.thr,pops=NULL){
    ## Match the conserved receptors with the ligands from the niche (mean-field)
    ## *** Input:
    ## fullExpTable -- data matrix (row-genes,column-cells, as returned by utils.R/prepare.data())
    ## or data file name (if large=T)
    ## enrich.str -- an enrichment struc
    ## conserv.thr -- a relative (fraction) cutoff for the receptors: if expression(rec) > conserv.thr, accept rec
    ## pops -- populations to be processed
    ## *** Output:
    ## sig.map -- saturation level for each receptor in a pop from each/all population of the niche
    ############################################################################
    id.r.list <- get.conserved.receptors(enrich.str,conserv.thr)
    if(is.null(pops))
        pops <- levels(as.factor(colnames(fullExpTable)))
    id.signals <- getIdentitySignals(fullExpTable,pops,id.r.list,lrmatchingmethod="none")
    sig.map <- deconvoluteSignalingBySubPop(fullExpTable,pops,id.signals,lrmatchingmethod="fracexp")
    return(sig.map)
}

set.envir <- function(org="MMU",LR.db="new"){
    ## Sets an R environment to work in: DB's and namings
    if( toupper(org) == "HSA" || toupper(org) == "MMU" ){
        cat("We use",toupper(org),"databases...\n")
        if( toupper(org) == "HSA" ){
            assign("ORG","hsapiens",envir = .GlobalEnv)
            ## TF-TF interactions
            load("HSA_TF_TF_interactions_2019-04-01.RData")
            assign("TFTF",HSA.TF.TF.interactions,envir = .GlobalEnv)
            load("HSA_background_interactome_2019-04-01.RData")
            assign("INTERACTOME",HSA.background.interactome,envir = .GlobalEnv)
        }
        if( toupper(org) == "MMU" ){
            assign("ORG","mmusculus",envir = .GlobalEnv)
            ## TF-TF interactions
            load("MMU_TF_TF_interactions_2019-04-01.RData")
            assign("TFTF",MMU.TF.TF.interactions,envir = .GlobalEnv)
            load("MMU_background_interactome_2019-04-01.RData")
            assign("INTERACTOME",MMU.background.interactome,envir = .GlobalEnv)
        }
        if(LR.db == "old")
            assign("LRDAT",LRDAT0[[paste0(toupper(org),"old")]],envir = .GlobalEnv)
        else
            assign("LRDAT",LRDAT0[[toupper(org)]],envir = .GlobalEnv)
        assign("RLIST",unique(LRDAT$Receptor),envir = .GlobalEnv)
        assign("LLIST",unique(LRDAT$Ligand),envir = .GlobalEnv)
        assign("TFDAT",TFDAT0[[toupper(org)]],envir = .GlobalEnv)
        assign("TFLIST",unique(TFDAT$Symbol),envir = .GlobalEnv)
        assign("PTHWDAT",PTHWDAT0[[toupper(org)]],envir = .GlobalEnv)
        assign("PTHWGENES",PTHWGENES0[[toupper(org)]],envir = .GlobalEnv)
        assign("KEGGPTHWGENES",PTHWGENES[grep("KEGG",names(PTHWGENES))],envir = .GlobalEnv)
        assign("KEGGITFLIST",unlist(lapply(KEGGPTHWGENES,function(x) x[x %in% TFLIST])),
               envir = .GlobalEnv)
        assign("PTHWNAMES",PTHWNAMES0[[toupper(org)]],envir = .GlobalEnv)
        eval(expression(names(PTHWGENES) <- PTHWNAMES),envir = .GlobalEnv)
    }
    else
        stop("Organism must be HSA or MMU.")
}

sc.commap <- function(fullExpTable,abs.exp.thr=0.0,conserv.thr=0.0,org="MMU",
                      tbl.svf=NULL,LR.db="new",verb=FALSE,
                      str.svf.suf=NULL,large=FALSE,large.sep='\t'){
    ## =================================
    ## Cell-by-cell L-R pairing function
    ## =================================
    ## Environments
    set.envir(org,LR.db=LR.db)
    ## Thresholds: for pathways and for L's
    cat("Abs.exp.thr:",abs.exp.thr,"\n")
    cat("Conserv.thr:",conserv.thr,"\n")
    ## Open file connection to dump the resulting table line-by-line
    if(!is.null(tbl.svf)){
        f <- file(tbl.svf,"w")
        writeLines(paste("Lig.pop","L.frac","Ligand","Receptor","R.frac","Rec.pop",sep='\t'),f)
    }
    ## Before processing the data check if the data is on disk (very large data were requested).
    if(large){
        suppressPackageStartupMessages({
            require(bigmemory)
        })
        cat("Large data reading ")
        cat("(",fullExpTable,")...")
        binFile <- "expTable.bin"
        descFile <- "expTable.desc"
        fullExpTable <- read.big.matrix(fullExpTable,sep=large.sep,header=TRUE,
                                        has.row.names=TRUE,type="double",backingfile=binFile,
                                        descriptorfile=descFile)
        cat(" Done.\n")
    }
    ## Populations
    pops <- levels(as.factor(colnames(fullExpTable)))
    if(large){
        fullExpTable.b <- make.exp.bool.bm(fullExpTable,exp.thr=abs.exp.thr)
        cat("Original class:",class(fullExpTable),"; bool class:",class(fullExpTable.b),"\n")
    }
    else
        fullExpTable.b <- make.exp.bool(fullExpTable,exp.thr=abs.exp.thr)
    
    ## ==========================
    ## Cell-by-cell R and L lists
    if(large)
        L.list.ext.all <- sc.L.list.ext.bm(fullExpTable.b)
    else
        L.list.ext.all <- sc.L.list.ext(fullExpTable.b)# L list for all pop
    if(is.null(tbl.svf))
        RL.table <- data.frame(L.pop=c(),L.frac=c(),Ligand=c(),
                               Receptor=c(),R.frac=c(),R.pop=c(),pthw.p=c(),id.p=c(),
                               pr=c(),rc=c(),f1=c())
    R.list.ext.pop <- vector("list",length=length(pops))
    names(R.list.ext.pop) <- as.character(pops)
    for(rp in pops){
        popExpTable.b <- select.pop.exp(fullExpTable.b,rp)
        cat("=== Receiver-pop:",rp,",",dim(popExpTable.b)[2],"cells\n")
        st <- system.time(R.list.ext.pop[[rp]] <- sc.RTF.pairing(popExpTable.b))
        ## ---------------------------------
        cat("* RTF pairing time:\n")
        print(st)
        R.tab <- get.all.rl.tab(R.list.ext.pop[[rp]],rp,mode="r")
        R.tab <- R.tab[R.tab > conserv.thr]
        r.list.sorted <- names(R.tab)
        st <- system.time(
        for(lp in pops){
            L.tab <- get.all.rl.tab(L.list.ext.all,lp,mode="l")
            l.list.sorted <- names(L.tab)
            lapply(r.list.sorted,function(x){
                cogn.l <- unique(LRDAT$Ligand[LRDAT$Receptor == x]) ## cognate ligands
                idx <- cogn.l %in% l.list.sorted ## cognate L's in L.tab
                cogn.l <- cogn.l[idx] ## cognate L's in the data
                cogn.l.perc <- L.tab[cogn.l] ## L.tab is named, so use L names
                for(i in seq_along(cogn.l)){
                    if(!is.null(tbl.svf)){
                        txt <- paste(lp,cogn.l.perc[i],cogn.l[i],x,R.tab[x],rp,sep='\t')
                        writeLines(txt,f)
                    }else {
                        RL.table <<- rbind(RL.table,data.frame(Lig.pop=lp,
                                                               L.frac=cogn.l.perc[i],
                                                               Ligand=cogn.l[i],
                                                               Receptor=x,
                                                               R.frac=R.tab[x],
                                                               Rec.pop=rp))
                    }
                }
            })
        })
        cat("* RL pairing time:\n")
        print(st)
    }
    ## Close the connection
    if(!is.null(tbl.svf)){
        close(f)
        cat("Saved to",tbl.svf,"...\n")
    }
    ## Save the RL-structures
    if(!is.null(str.svf.suf)){
        tm <- as.character(Sys.time())
        dt <- unlist(strsplit(tm," "))
        fn <- paste("pthw_data",str.svf.suf,dt[1],paste(unlist(strsplit(dt[2],":"))[1:2],collapse="."),
                    sep="_")
        st <- system.time({
            cat("Saving R-L structures to disk (",paste0(fn,".RData"),")...\n",sep="")
            cat("R.list:",object.size(R.list.ext.pop)/1000000.0,"MB\n")
            cat("L.list:",object.size(L.list.ext.all)/1000000.0,"MB\n")
	    ## Save with version=2 for R.ver 3.4 and lower compatibility
            save(R.list.ext.pop,L.list.ext.all,file=paste0(fn,".RData"),version=2)})
        cat("* Saving time:\n")
        print(st)
    }
    cat("=================\n")
    cat(length(pops),"populations analyzed:",pops,"\n")
    if(large){
        rm(fullExpTable);void <- gc(FALSE) ## recommended in bigmemory
        ## Remove backing files
        unlink(binFile)
        unlink(descFile)
    }
    if(is.null(tbl.svf)){
        if(dim(RL.table)[1] > 0)
            rownames(RL.table) <- seq(1,dim(RL.table)[1])
        return(RL.table)
    }
    else{
        return(NULL)
    }
}

