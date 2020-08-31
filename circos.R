suppressPackageStartupMessages(require(circlize))

## Process cMap df
x.cell.margins <- function(cMap,varName,factor){
    x.cell.margin <- list()

    pops <- unique(cMap[factor][,1])
    for (p in pops){
#        cat(p,": \n")
        var.list <- unique(cMap[varName][cMap[factor] == p])
        n.chn.by.var <- vector("integer",length=length(var.list))
        names(n.chn.by.var) <- var.list

        for (r in var.list){
            n.chn.by.var[r] <- sum(cMap[varName] == r & cMap[factor] == p)
        }
#       print(n.chn.by.var)
        x.cell.margin[[p]] <- c(0,cumsum(n.chn.by.var / sum(n.chn.by.var)))
        names(x.cell.margin[[p]]) <- c("NA",as.character(var.list))
    }

    return(x.cell.margin)
}

circMap <- function(cMap,set.colors=T){
    cMap <- lapply(cMap,as.character)
    cMap <- data.frame(Receptor=cMap$Receptor,Rec.pop=cMap$Rec.pop,
                       Ligand=cMap$Ligand,Lig.pop=cMap$Lig.pop,
                       stringsAsFactors=FALSE)

    pop.r <- as.character(cMap$Rec.pop)
    pop.l <- as.character(cMap$Lig.pop)
    pops <- unique(c(pop.r,pop.l))

    circos.par("track.height" = 0.2)
    circos.initialize(factors = pops, xlim = c(0,1))
    x.cell.rec.margin <- x.cell.margins(cMap,"Receptor","Rec.pop")
    x.cell.lig.margin <- x.cell.margins(cMap,"Ligand","Lig.pop")

    col.palette <- c("forestgreen","tomato4","dodgerblue","mediumpurple",
                     "sandybrown","firebrick","plum","limegreen",
                     "coral1","turquoise3","red","blue","green")
    if (set.colors){
        cols <- col.palette[1:length(pops)]
    }
    else {
        cols <- sample(col.palette,length(pops))
    }

    circos.track(ylim=c(0,1),track.height=0.25,panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + uy(6, "mm"), 
                    CELL_META$sector.index,col=cols[CELL_META$sector.index==pops])
        circos.axis(labels.cex = 0.6)
        margins <- x.cell.lig.margin[[CELL_META$sector.index]]
        if(!is.null(margins)){
            for(i in 1:(length(margins)-2)) ## sectors within cells
                circos.lines(rep(margins[i+1],2),c(0,1))
            for(i in 1:(length(margins)-1)){ ## Links and labels
                p1 <- margins[i]+(margins[i+1]-margins[i])/2
                p1 <- c(p1-(margins[i+1]-margins[i])*0.08,p1+(margins[i+1]-margins[i])*0.08)
                ##p1 <- c(p1-0.005,p1+0.005)
                if( (margins[i+1]-margins[i]) > 0.1 )
                    circos.text(p1,0.5,names(margins)[i+1],niceFacing=T)
                else
                    circos.text(p1,0.5,names(margins)[i+1],facing="clockwise",niceFacing=T)
                ## Only unique links, not accounting for the width (#TF,#channels)
                ## Receptor populations == sectors, coupled with the L in question
                rec.pops <- unique(as.character(cMap$Rec.pop[cMap$Ligand == names(margins)[i+1]]))
                for( s in 1:length(rec.pops) ){ ## receiving cells
                    ## List of R's cognate to the L in question
                    rec.list <- unique(cMap$Receptor[cMap$Ligand == names(margins)[i+1] & cMap$Rec.pop == rec.pops[s]])
                    ##                cat(rec.list,"\n")
                    ## R margin array, names = all R's for the Receptor sector/population
                    r.margin <- x.cell.rec.margin[[rec.pops[s]]]
                    ##              cat(names(r.margin),"\n")
                    for ( r in 1:length(rec.list) ){## receiving receptor sectors
                        idx <- which(names(r.margin) == rec.list[r])
                        p2 <- r.margin[idx-1] + (r.margin[idx] - r.margin[idx-1])/2
                        ##                cat(names(r.margin)[idx],p2,"\n")
                        circos.link(CELL_META$sector.index,p1,rec.pops[s],p2,directional=1,
                                    h.ratio=0.5,arr.length=.7,arr.col=cols[rec.pops[s]==pops],
                                    rou1=0.75,rou2=0.5,
                                    border=cols[CELL_META$sector.index==pops],
                                    lwd=2.5,col=cols[CELL_META$sector.index==pops])
                    }
                }
            }
        }
        else {
            circos.lines(c(0,1),c(0,1),lty=2)
            circos.lines(c(0,1),c(1,0),lty=2)

        }
    },bg.col=cols)

    circos.track(ylim=c(0,1),track.height=0.2,panel.fun = function(x, y) {
        margins <- x.cell.rec.margin[[CELL_META$sector.index]]
        if(!is.null(margins)){
            for(i in 1:(length(margins)-2))
                circos.lines(rep(margins[i+1],2),c(0,1))
            for(i in 1:(length(margins)-1)){
                if( (margins[i+1]-margins[i]) < 0.1 )
                    circos.text(margins[i]+(margins[i+1]-margins[i])/2,0.5,
                                names(margins)[i+1],facing="clockwise",niceFacing=T)
                else
                    circos.text(margins[i]+(margins[i+1]-margins[i])/2,0.5,
                                names(margins)[i+1],niceFacing=T)
            }
        }
        else {
            circos.lines(c(0,1),c(0,1),lty=2)
            circos.lines(c(0,1),c(1,0),lty=2)
        }
    },bg.col=cols)



}
