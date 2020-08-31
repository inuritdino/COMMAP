## ================================
## SYNERGY
## ================================

#####################################################################################
### Statistical "synergy": co-association of two (or more receptors) in a population
### No common iTF is needed, a broader approach
#####################################################################################

association.tbl <- function(enrich.str,conserv.thr){
    ## Get all characteristic receptors
    r.lst <- get.conserved.receptors(enrich.str,conserv.thr,uniprot=FALSE)
    asc.tbls <- vector("list",length(r.lst))
    ## Form all pairs of receptors within each population
    for(ip in 1:length(r.lst)){
        asc.tbls[[ip]] <- association.tbl.r.lst(enrich.str[[ip]],r.lst[[ip]])
    }
    names(asc.tbls) <- names(r.lst)
    return(asc.tbls)
}

association.tbl.r.lst <- function(enrich.str.pop,r.lst){
    ## Process each pop R list here
    combs <- combn(1:length(r.lst),2) ## take combinations by 2
    asc.tbl <- vector("list",dim(combs)[2])
    full.cell.idx <- 1:length(enrich.str.pop) ## all cells in the pop
    cell.id.buf <- vector("list",2)
    for(ic in 1:dim(combs)[2]){
        cell.id.buf[[1]] <- get.cell.id(enrich.str.pop,rec = r.lst[ combs[1,ic] ])[[1]]
        cell.id.buf[[2]] <- get.cell.id(enrich.str.pop,rec = r.lst[ combs[2,ic] ])[[1]]
        asc.tbl[[ic]] <- matrix(c(length(intersect(cell.id.buf[[1]],cell.id.buf[[2]])),## r1 & r2
                           length(intersect(cell.id.buf[[1]],setdiff(full.cell.idx,cell.id.buf[[2]]))),## r1 & not-r2
                           length(intersect(cell.id.buf[[2]],setdiff(full.cell.idx,cell.id.buf[[1]]))),## r2 & not-r1
                           length(intersect(setdiff(full.cell.idx,cell.id.buf[[2]]),
                                            setdiff(full.cell.idx,cell.id.buf[[1]])))## not-r1 & not-r2
                           ),nrow=2,dimnames=list(paste(r.lst[combs[2,ic]],c("ON","OFF"),sep="_"),
                                                  paste(r.lst[combs[1,ic]],c("ON","OFF"),sep="_")))
        names(asc.tbl)[ic] <- paste(r.lst[ combs[,ic] ],collapse="_")
    }
    return(asc.tbl)
}

association.prob <- function(asc.tbl,mode="both",least.asc=FALSE){
    ## min to be optimistic, max to be more conservative
    if(asc.tbl[1,1] == 0)
        return(0.0)
    if(least.asc){
        p0 <- least.association.prob(asc.tbl,mode=mode) ## min n1/p
        ##p1 <- min(asc.tbl[1,1]+asc.tbl[1,2],asc.tbl[1,1]+asc.tbl[2,1])## max n1/p
        p.min <- asc.tbl[1,1]/(asc.tbl[1,1] + min(asc.tbl[1,2],asc.tbl[2,1]))
        p.max <- asc.tbl[1,1]/(asc.tbl[1,1] + max(asc.tbl[1,2],asc.tbl[2,1]))
        p1 <- 1.0
        if(mode == "min")
            p <- p.min
        else if(mode == "max")
            p <- p.max
        else if(mode == "both"){
            p <- p.max + p.min
            p1 <- 2.0
        }
        else{
            cat("Warning: mode is not explicit. Both assumed.\n")
            p <- p.max + p.min
            p1 <- 2.0
        }
        ##cat("p0 =",p0,"p =",p,"p1 =",p1,"\n")
        return((p - p0)/(p1 - p0))
    }
    else{
        p.min <- asc.tbl[1,1]/(asc.tbl[1,1] + min(asc.tbl[1,2],asc.tbl[2,1]))
        p.max <- asc.tbl[1,1]/(asc.tbl[1,1] + max(asc.tbl[1,2],asc.tbl[2,1]))
        if(mode=="min")
            return(p.min)
        else if(mode == "max")
            return(p.max)
        else if(mode == "both")
            return(p.max + p.min)
        else{
            cat("Warning: mode is not explicit. Both assumed.\n")
            return(p.max + p.min)
        }
    }
}

least.association.tbl <- function(asc.tbl){
    N <- sum(asc.tbl)
    n1 <- sum(asc.tbl[1,]) + sum(asc.tbl[,1]) - N
    if(n1 < 0)
        return(matrix(c(0,0,0,0),nrow=2))
    n2 <- sum(asc.tbl[1,]) - n1
    n3 <- sum(asc.tbl[,1]) - n1
    n4 <- 0
    return(matrix(c(n1,n2,n3,n4),nrow=2,byrow=TRUE))
}

least.association.prob <- function(asc.tbl,mode="both"){
    least.asc.tbl <- least.association.tbl(asc.tbl)
    if(sum(least.asc.tbl) == 0)
        return(0.0)
    else
        return(association.prob(least.asc.tbl,mode=mode,least.asc=FALSE))
}

association.phi.max <- function(asc.tbl){
    n11x <- asc.tbl[1,1]
    n00x <- asc.tbl[2,2]
    n10x <- asc.tbl[1,2]
    n01x <- asc.tbl[2,1]
    ## Define max correlation structure possible:
    ## max overlap that may exist for the given configuration
    n11 <- min(n11x+n10x,n11x+n01x)## max n1/p
    n10 <- n10x + n11x - n11
    n01 <- n01x + n11x - n11
    n00 <- sum(asc.tbl) - n11 - n10 - n01
    if( (n11*n00 - n10*n01) == 0 )
        return(0.0)
    phi <- (n11*n00 - n10*n01)/(( as.numeric((n11+n10)*(n01+n00))*as.numeric((n10+n00)*(n11+n01)) )**(0.5))
    return(phi)
}

association.phi.min <- function(asc.tbl){
    return(association.phi(least.association.tbl(asc.tbl)))
}

association.phi <- function(asc.tbl){
    n11 <- asc.tbl[1,1]
    n00 <- asc.tbl[2,2]
    n10 <- asc.tbl[1,2]
    n01 <- asc.tbl[2,1]
    if( (n11*n00 - n10*n01) == 0 )
        return(0.0)
    phi <- (n11*n00 - n10*n01)/(( as.numeric((n11+n10)*(n01+n00))*as.numeric((n10+n00)*(n11+n01)) )**(0.5))
    return(phi)
}

association.phi.corr <- function(asc.tbl){
    ## x <- association.phi.min(asc.tbl)
    y <- association.phi.max(asc.tbl)
    z <- association.phi(asc.tbl)
    if(abs(y) < 1e-04)
        return(z)
    if(z < 0)
        return(z)
    else
        return(z/y)
    ## return((z-x)/(y-x))
}

#####################################################################
### Synergy: combinations of R's with one or more common iTF targets
### Contingency table approach (at least for combinations of two)
#####################################################################

get.ph.tbl.idx <- function(ph.tbls,ct,recs){
    rec.lst <- rec.list(ph.tbls)
    pop.lst <- pop.list(ph.tbls)
    Reduce(c,lapply(recs,function(x) which(rec.lst == x & pop.lst == ct)))
}
get.recs.common.iTF <- function(ph.tbls,ct,recs){
    idx <- get.ph.tbl.idx(ph.tbls,ct,recs)
    comm.iTFs <- Reduce(intersect,lapply(ph.tbls[idx],function(x) x$tbl$if.TF))
    return(comm.iTFs)
}

get.rec.combs <- function(ph.tbls,ct,by=2){
    rec.lst <- rec.list(ph.tbls)
    pop.lst <- pop.list(ph.tbls)
    rec.lst.ct <- rec.lst[pop.lst == ct]
    combs <- combn(1:length(rec.lst.ct),by)
    rec.combs <- matrix(NA,nrow=0,ncol=by+1,dimnames=list(c(),c(paste0("Receptor",1:by),"iTF")))
    for(cmb in 1:dim(combs)[2]){
        ##cat("Rec pair:",rec.lst.ct[combs[,c]],"\n")
        comm.iTFs <- get.recs.common.iTF(ph.tbls,ct,rec.lst.ct[combs[,cmb]])
        ##cat("\t Comm iTFs:",comm.iTFs,"\n")
        ## if(length(comm.iTFs) == 0){
        ##     cat("Zero length:",rec.lst.ct[combs[,cmb]],"\n")
        ## }
        for(itf in seq_along(comm.iTFs))
            rec.combs <- rbind(rec.combs,c(rec.lst.ct[combs[,cmb]],comm.iTFs[itf]))
    }
    return(rec.combs)
}

get.conting.tbl <- function(ph.tbls,ct,rec.lst,iTF,norm=TRUE,rand=FALSE,n.permut=100){
    if(length(rec.lst) > 2)
        cat("Warning: more than 2 receptors requested...\n")
    
    ph.tbl.idx <- get.ph.tbl.idx(ph.tbls,ct,rec.lst) ## which ph.tbls to use
    n.cells <- ph.tbls[[ph.tbl.idx[1]]]$nc ## 1st R's num of cells (nc field) to use, cell type is the same
    full.cell.idx <- seq(1,n.cells,by=1)

    ## Null distribution, permutations
    if(rand){
        rec.n.cells <- lapply(ph.tbls[ph.tbl.idx],function(x)
            length(unlist(num.arr(x$tbl$cell.id[x$tbl$if.TF == iTF]))))
        tbl <- vector("list",n.permut)
        for(p in 1:n.permut){
            cell.ids.perm <- lapply(rec.n.cells,function(n) sample(full.cell.idx,n))
            ## print(lapply(cell.ids.perm,sort))
            ch.arr <- lapply(rec.lst,function(x) rep("out",times=length(full.cell.idx)))
            for(i in 1:length(rec.lst)){
                ch.arr[[i]][ cell.ids.perm[[i]] ] <- "in"
            }
            if(norm)
                tbl[[p]] <- table(ch.arr)/n.cells
            else
                tbl[[p]] <- table(ch.arr)
            tbl[[p]] <- matrix(c(tbl[[p]]["in","in"],tbl[[p]]["out","in"],
                                 tbl[[p]]["in","out"],tbl[[p]]["out","out"]),
                                nrow=2,dimnames=list(c("in","out"),c("in","out")))
        }
    }
    else {## Normal usage
        ch.arr <- lapply(rec.lst,function(x) rep("out",times=length(full.cell.idx)))
        cell.ids <- lapply(ph.tbls[ph.tbl.idx],function(x)
            unlist(num.arr(x$tbl$cell.id[x$tbl$if.TF == iTF])))
        for(i in 1:length(rec.lst)){
            ch.arr[[i]][ cell.ids[[i]] ] <- "in"
        }
        ## print(ch.arr)
        if(norm)
            tbl <- table(ch.arr)/n.cells
        else
            tbl <- table(ch.arr)
        tbl <- matrix(c(tbl["in","in"],tbl["out","in"],
                        tbl["in","out"],tbl["out","out"]),nrow=2,dimnames=list(c("in","out"),c("in","out")))
    }
    return(tbl)
}

### Measure of synergy
mut.info <- function(contig.tbl,norm=FALSE){
    x1 <- rowSums(contig.tbl)
    x2 <- colSums(contig.tbl)
    H1 <- -sum(sapply(x1,function(x) ifelse(x > 0, x*log2(x), 0)))
    H2 <- -sum(sapply(x2,function(x) ifelse(x > 0, x*log2(x), 0)))
    H12 <- -sum(sapply(contig.tbl,function(x) ifelse(x > 0, x*log2(x), 0)))
    ## print(H1+H2-H12)
    ## print(mi.empirical(contig.tbl,unit="log2"))
    if(norm)
        return((H1+H2-H12)/max(H1,H2))
    else
        return(H1+H2-H12)
}

cond.prob <- function(conting.tbl,cond="in"){
    p1 <- conting.tbl[cond,cond]/sum(conting.tbl[cond,])
    p2 <- conting.tbl[cond,cond]/sum(conting.tbl[,cond])
    return(ifelse(p1 > p2, p1, p2))
}

phi.coeff <- function(conting.tbl){
    n11 <- conting.tbl["in","in"]
    n00 <- conting.tbl["out","out"]
    n10 <- conting.tbl["in","out"]
    n01 <- conting.tbl["out","in"]
    phi <- (n11*n00 - n10*n01)/(( as.numeric((n11+n10)*(n01+n00))*as.numeric((n10+n00)*(n11+n01)) )**(0.5))
    return(phi)
}

yule.y <- function(conting.tbl){
    n11 <- conting.tbl["in","in"]
    n00 <- conting.tbl["out","out"]
    n10 <- conting.tbl["in","out"]
    n01 <- conting.tbl["out","in"]
    Y <- (sqrt(n11*n00) - sqrt(n10*n01))/(sqrt(n11*n00) + sqrt(n10*n01))
    return(Y)
}

stat.tests <- function(con.tbl,test="fisher"){
    if(test == "fisher"){
        pval <- fisher.test(con.tbl,alternative="greater")$p.value
    }
    return(pval)
}


### Overall procedures
## Stem plot
make.stem.plot <- function(plt.df,x=x,y=y,p.line=FALSE){
    suppressPackageStartupMessages({
        require(ggplot2)
    })
    f <- function(y)
        c(label=length(y),y=min(y))
    
    g <- ggplot(plt.df,aes(x=plt.df[[x]],y=plt.df[[y]])) +
        ## geom_segment(aes(x=plt.df[[x]],xend=plt.df[[x]],y=0,yend=plt.df[[y]]),color="grey") +
        geom_boxplot() +
        geom_dotplot(binaxis='y', stackdir='center', dotsize=0.2) +
        geom_jitter(shape=8, position=position_jitter(0.15), color='orange') +
        ## geom_point(color="orange",size=4) +
        theme_light() +
        coord_flip() +
        theme(
            panel.grid.major.x = element_blank(),
            panel.border = element_blank(),
            axis.ticks.x = element_blank()
        ) +
        stat_summary(fun.data=f, geom="text", vjust=-0.3, hjust=0.01, col="blue") + 
        xlab("") +
        ylab(y)
    if(p.line)
        g <- g + geom_segment(aes(x=0,xend=length(plt.df[[x]]),y=-log10(0.05),
                                  yend=-log10(0.05)),color="black")
    return(g)
}

multi.stem.plot <- function(ph.tbls,ct,by=2){
    rec.combs <- get.rec.combs(ph.tbls,ct,by) ## receptor combinations with a common TF
    u.rec.combs <- apply(unique(rec.combs[,1:2]),1,paste,collapse="_")
    if(length(u.rec.combs) > 40){
        cat("Warning: more than 30 pairs, reducing...\n")
        rec.combs <- rec.combs[sample(1:dim(rec.combs)[1],30),]
        u.rec.combs <- apply(unique(rec.combs[,1:2]),1,paste,collapse="_")
    }
    ## plt.df <- data.frame(comb=u.rec.combs,phi=rep(NA,length(u.rec.combs)),
    ##                      neg.log.p.fisher=rep(NA,length(u.rec.combs)),
    ##                      max.cond.pr=rep(NA,length(u.rec.combs)),
    ##                      mut.info=rep(NA,length(u.rec.combs)),
    ##                      yule.y=rep(NA,length(u.rec.combs)),
    ##                      n.log.mi.p=rep(NA,length(u.rec.combs)),
    ##                      n.log.phi.p=rep(NA,length(u.rec.combs)),
    ##                      n.TF=rep(0,length(u.rec.combs)))
    all.combs <- apply(rec.combs,1,function(x) paste(x[1:2],collapse="_"))
    plt.df <- data.frame(comb=all.combs,phi=rep(NA,length(all.combs)),
                         neg.log.p.fisher=rep(NA,length(all.combs)),
                         min.cond.pr=rep(NA,length(all.combs)),
                         mut.info=rep(NA,length(all.combs)),
                         yule.y=rep(NA,length(all.combs)))

    for(i in 1:dim(rec.combs)[1]){
        ## +++++++++++++++++++++
        ## Get contingency table
        cnt.tbl <- get.conting.tbl(ph.tbls,ct,rec.combs[i,1:2],rec.combs[i,3],norm=F)
        ## cnt.tbl.perm <- get.conting.tbl(ph.tbls,ct,rec.combs[i,1:2],rec.combs[i,3],norm=F,rand=T,n.permut=200)
        ## ++++++++++++++++++++++++++++++++++++
        ## Index of the combination in question
        idx <- which(plt.df$comb == paste(rec.combs[i,1:2],collapse="_"))
        idx <- idx[ which(is.na(plt.df$phi[idx]))[1] ]

        ## +++++++++++++++++++++
        ## +++ Calculate Phi +++
        phi <- phi.coeff(cnt.tbl)
        plt.df$phi[idx] <- phi
        ## phi.perm <- unlist(lapply(cnt.tbl.perm,phi.coeff))
        ## phi.p <- sum(phi.perm >= phi)/length(phi.perm)
        ## plt.df$n.log.phi.p[idx] <- ifelse(is.na(plt.df$n.log.phi.p[idx]),
        ##                            ifelse(phi.p>0,-log10(phi.p),3),
        ##                            ifelse(phi > plt.df$phi[idx],
        ##                            ifelse(phi.p>0,-log10(phi.p),3),plt.df$n.log.phi.p[idx]))

        ## +++++++++++++++++++++++++++++++++++++
        ## +++ Calculate Fisher test p-value +++
        fisher.p <- -log10(stat.tests(cnt.tbl,test="fisher"))
        plt.df$neg.log.p.fisher[idx] <- fisher.p

        ## +++++++++++++++++++++++++++++++++++++++
        ## +++ Calculate max. conditional prob +++
        min.cond.pr <- cond.prob(cnt.tbl,cond="in")
        plt.df$min.cond.pr[idx] <- min.cond.pr

        ## +++++++++++++++++++++++++++++
        ## +++ Calculate mutual info +++
        mi <- mut.info(cnt.tbl/sum(cnt.tbl))
        plt.df$mut.info[idx] <- mi
        ## mi.perm <- unlist(lapply(cnt.tbl.perm,function(x) mut.info(x/sum(x)))) # mi null distr
        ## mi.p <- sum(mi.perm >= mi)/length(mi.perm) # p-value
        ## How to select an appropriate p-value for the mi selected before?
        ## plt.df$n.log.mi.p[idx] <- ifelse(is.na(plt.df$n.log.mi.p[idx]),ifelse(mi.p>0,-log10(mi.p),3),
        ##                           ifelse(mi > plt.df$mut.info[idx],
        ##                                  ifelse(mi.p>0,-log10(mi.p),3),plt.df$n.log.mi.p[idx]))

        ## ++++++++++++++++++++++++
        ## +++ Calculate Yule Y +++
        yy <- yule.y(cnt.tbl/sum(cnt.tbl))
        ## Change the value in df
        plt.df$yule.y[idx] <- yy
    }
    suppressPackageStartupMessages({
        require(gridExtra)
    })
    p1 <- make.stem.plot(plt.df,"comb","phi")
    ##p7 <- make.stem.plot(plt.df,"comb","n.log.phi.p",p.line=T)
    p2 <- make.stem.plot(plt.df,"comb","neg.log.p.fisher")
    p3 <- make.stem.plot(plt.df,"comb","min.cond.pr")
    p4 <- make.stem.plot(plt.df,"comb","mut.info")
    ##p6 <- make.stem.plot(plt.df,"comb","n.log.mi.p",p.line=T)
    p5 <- make.stem.plot(plt.df,"comb","yule.y")
    p8 <- make.stem.plot(plt.df,"comb","n.TF")
    grid.arrange(p1,p3,p2,p4,p5,nrow=2,ncol=3,top=ct)
    print(plt.df)
}

## Distance
multi.stat.dist <- function(ph.tbls,ct,by=2){
    if(by != 2)
        return(NULL)
    rec.combs <- get.rec.combs(ph.tbls,ct,by) ## receptor combinations with a common TF
    u.recs <- unique(c(rec.combs[,1],rec.combs[,2]))
    dist.mat <- matrix(0.0,nrow=length(u.recs),ncol=length(u.recs),dimnames=list(u.recs,u.recs))
    for(i in 1:dim(rec.combs)[1]){
        contig <- get.conting.tbl(ph.tbls,ct,rec.combs[i,1:2],rec.combs[i,3],norm=F)
        print(1-cond.prob(contig,cond="in"))
        dist.mat[rec.combs[i,2],rec.combs[i,1]] <- max(dist.mat[rec.combs[i,2],rec.combs[i,1]],
                                                       stat.tests(contig),na.rm=TRUE)
        ## dist.mat[rec.combs[i,1],rec.combs[i,2]] <- 1-cond.prob(contig,cond="in")
    }
    ##return(as.dist(dist.mat,diag=T,upper=T))
    return(dist.mat)
}

##############################################################################################
### Synergy: graph-based approach: determine indegree of iTFs to select synergy modules of R's
### Common iTF approach
##############################################################################################

get.connect.table <- function(ph.tbls,ct){
    if.TFs <- lapply(ph.tbls,function(x) x$tbl$if.TF)
    c.fracs <- lapply(ph.tbls,function(x) x$tbl$cell.frac)
    c.types <- pop.list(if.TFs)
    ct.iTFs <- if.TFs[c.types == ct]
    ct.c.fracs <- c.fracs[c.types == ct]
    idx <- sapply(ct.iTFs,length) != 0
    ct.iTFs <- ct.iTFs[idx] ## Remove empty lists (because of phenoTF.table filtering)
    ct.c.fracs <- ct.c.fracs[idx]
    tbl.size <- sum(sapply(ct.iTFs,length))

    connect.table <- data.frame(Source=rep(NA,tbl.size),Target=rep(NA,tbl.size),Cell.frac=rep(NA,tbl.size))
    rec.lst <- rec.list(ct.iTFs)
    start.idx <- 1
    rec.pop.names <- names(ct.iTFs)
    for(i in 1:length(ct.iTFs)){
        old.start.idx <- start.idx
        ## Sorted table: one receptor gets all its targets first, then the second etc.
        for(j in 1:length(ct.iTFs[[i]])){
            connect.table$Target[start.idx] <- ct.iTFs[[i]][j]
            connect.table$Cell.frac[start.idx] <- ct.c.fracs[[i]][j]
            start.idx <- start.idx + 1
        }
        connect.table$Source[old.start.idx:(start.idx-1)] <- rep(strsplit(rec.pop.names[i],"_")[[1]][1],start.idx-old.start.idx)
    }
    return(connect.table)
}

get.cell.id.by.rec.iTF <- function(ph.tbls,ct,rec,iTF){
    ph.tbl <- ph.tbls[paste(rec,ct,sep="_") == names(ph.tbls)][[1]]
    c.idx <- num.arr(ph.tbl$tbl$cell.id[ph.tbl$tbl$if.TF == iTF],",")[[1]]
    n.cells <- ph.tbl$nc
    return(list(c.idx,n.cells))
}

plot.syn.graph <- function(con.tbl){
    suppressPackageStartupMessages({
        require(igraph)
    })
    gr <- graph_from_data_frame(con.tbl)
    E(gr)$width <- E(gr)$Cell.frac*10
    V(gr)$Receptor <- ifelse(V(gr)$name %in% con.tbl$Source,TRUE,FALSE)
    V(gr)$shape <- c("square","circle")[V(gr)$Receptor+1]
    V(gr)$color <- c("steel blue","orange")[V(gr)$Receptor+1]
    plot(gr, vertex.label.color="black",layout=layout_components)
}

get.syn.graph.by.iTF.indegree <- function(con.tbl,cmp.fxn,indeg=2,plot=TRUE){
    ## Count all iTFs
    tbl <- table(con.tbl$Target)
    con.tbl.2 <- subset(con.tbl, Target %in% names(tbl[cmp.fxn(tbl,indeg)]))
    if(dim(con.tbl.2)[1] == 0)
        return(NULL);
    if(plot)
        plot.syn.graph(con.tbl.2)
    
    return(con.tbl.2);
}

get.syn.modules <- function(con.tbl){
    ## Reduced con.tbl, only one indegree is allowed
    tbl <- table(con.tbl$Target)
    if(any(tbl[1] != tbl)){
        cat("Warning: indegree differs, run get.?.by.iTF.indegree first\n")
        return(NULL)
    }
    indeg <- tbl[1]
    iTFs <- names(tbl)
    module.df <- data.frame(rec=c(),iTF=c())
    for(i in 1:length(tbl)){
        recs <- paste(sort(con.tbl$Source[con.tbl$Target == iTFs[i]]),collapse=",")
        idx <- module.df$rec == recs
        if(any(idx)){
            module.df$iTF[idx] <- paste(module.df$iTF[idx],iTFs[i],sep=",")
        }
        else{
            module.df <- rbind(module.df,data.frame(rec=recs,iTF=iTFs[i],stringsAsFactors=FALSE))
        }
    }
    return(module.df)
}

get.cfrac.by.module.iTF <- function(ph.tbls,ct,module.df){
    for(i in 1:dim(module.df)[1]){
        cat("Module:",module.df[i,1],"\n")
        iTFs <- char.arr(module.df[i,2])[[1]]
        recs <- char.arr(module.df[i,1])[[1]]
        for(j in 1:length(iTFs)){
            cat("\t iTF:",iTFs[j],"=> ")
            c.idx.all <- vector("list",length(recs))
            cat("(")
            for(k in 1:length(recs)){
                if(k == 1){
                    x <- get.cell.id.by.rec.iTF(ph.tbls,ct,recs[k],iTFs[j])
                    c.idx.all[[k]] <- x[[1]]
                    n.cells <- x[[2]]
                } else {
                    c.idx.all[[k]] <- get.cell.id.by.rec.iTF(ph.tbls,ct,recs[k],iTFs[j])[[1]]
                }
                cat(recs[k],":",length(c.idx.all[[k]])/n.cells,",",sep="")
            }
            cat("),n.cells:",n.cells,"\n")
            y <- Reduce(intersect,c.idx.all,accumulate=FALSE)
            cat("\t\tTotal synergy strength:",length(y)/n.cells,"\n")
        }
    }
}

synergy.modules <- function(ph.tbls,ct,indegree=2){
    con.tbl <- get.connect.table(ph.tbls,ct)
    con.tbl.2 <- get.syn.graph.by.iTF.indegree(con.tbl,`==`,indegree,plot=FALSE)
    modules <- get.syn.modules(con.tbl.2)
    get.cfrac.by.module.iTF(ph.tbls,ct,modules)
    par(mfrow=c(1,2))
    plot.syn.graph(con.tbl)
    plot.syn.graph(con.tbl.2)
}
