######################################
### Mean-field Ligand pairing routines
######################################

get.id.receptor.list <- function(receptors,rec.pops){
    ## Prepare reported ID receptors for mean-field routines.
    ## ID receptors are those that have passed all selection criteria:
    ## pathway enrichment, TF's criteria etc.
    pops <- unique(rec.pops)
    id.r.list <- vector("list",length = length(pops))
    for(p in seq_along(pops)){
        id.r.list[[p]] <- unique(receptors[rec.pops == pops[p]])
    }
    names(id.r.list) <- pops
    return(id.r.list)
}


#################
# Computes the matching between a receptor and its cognate ligands for all subpopulations and
# their id receptors. Currently only use "fracexp". Does not yet uses the free ligand expression
# Could be tested later depending on the results.
#################
deconvoluteSignalingBySubPop <- function(data,unique_cts,id_signals, id_r_list = NULL, removeCompetingReceptors = FALSE, lrmatchingmethod = "expression"){
    ## =====================================
    ## I add to account for global variables
    data_r <- data[rownames(data) %in% RLIST,]
    data_l <- data[rownames(data) %in% LLIST,]

    influenceMatrices <- list()
    for(to in seq_along(unique_cts)){
        ##print(to)
        cur_rec <- rownames(data_r)[which(rownames(data_r) %in% unique(id_signals[[to]][,2]))]
        cur_InfMat <- matrix(-1,nrow = length(id_signals),ncol = length(cur_rec))
        k <- 1
        for(from in seq_along(unique_cts)){
            ## ==============================
            ## Removed all other lrmatchingmethod's here
            ## ...
            ## lrmatchingmethod == "fracexp"
            ## cells_from <- which(colnames(data) %in% rownames(cts_in_data)[which(cts_in_data$cell.type %in% unique_cts[from])])
            cells_from <- colnames(data) == unique_cts[from]
            ## cells_to <- which(colnames(data) %in% rownames(cts_in_data)[which(cts_in_data$cell.type %in% unique_cts[to])])
            cells_to <- colnames(data) == unique_cts[to]
            fracs <- vector(mode = "numeric", length = length(rownames(data_r)[which(rownames(data_r) %in% unique(id_signals[[to]][,2]))]))
            fracsum <- 0
            i_fracs <- 1
            for(i in seq_along(cur_rec)){
                r <- cur_rec[i]
                ligands_r <- unique(id_signals[[to]][which(id_signals[[to]][,2] == r),1])
                ligands_r <- ligands_r[which(ligands_r %in% rownames(data_l))]
                ##print(head(ligands_r))
                ligands_r_sum <- apply(data_l[ligands_r,cells_from,drop=FALSE],2,sum)
                ##print(head(rownames(data_l)[ligands_r]))
                ## rownames(data_l)[ligands_r] are not valid entities
                ## What is meant here???
                ## Is it just all ligands that are present in the data???
                ## Then LRDAT$Ligand %in% rownames(data_l) should suffice
                ## BUT competingReceptors are never used really
                ## competingReceptors <- setdiff(unique(LRDAT$Receptor[which(LRDAT$Ligand %in% rownames(data_l)[ligands_r])]),r)
                ## competingReceptorLigands <- setdiff(unique(LRDAT$Ligand[which(LRDAT$Receptor %in% competingReceptors)]),rownames(data_l)[ligands_r])
                ## freecompligand_expr <- max(sum(unname(apply(data_r[competingReceptors,],1,sum)))-sum(unname(apply(data_l[competingReceptorLigands,],1,sum))),0)
                ## freecompligand_expr currently unused. Takes into account proportions, so is problematic...
                ##freecompligand_expr <- freecompligand_expr/ncol(data_r)
                frac <- ((sum(as.numeric(unname(ligands_r_sum))))/sum(as.numeric(unname(data_r[r,cells_to]))))*(length(cells_to)/(1.0*length(cells_from)))
                fracsum <- fracsum + frac
                fracs[i_fracs] <- frac
                i_fracs <- i_fracs + 1
            }
            ##print(fracs)
            ##influenceMat[k,to] <- fracsum/length(rownames(data_r)[which(rownames(data_r) %in% unique(id_signals[[to]][,2]))])
            cur_InfMat[from,] <- fracs
        }
        colnames(cur_InfMat) <- cur_rec
        rownames(cur_InfMat) <- unique_cts
        influenceMatrices[[to]] <- cur_InfMat
        k <- k + 1
    }
    names(influenceMatrices) <- unique_cts
    return(influenceMatrices)
}

###########
# Get all ligands for the identity receptors.
# Multiple methods have been implemented, but so far only a workaround should be used here.
# The workaround is lrmatchingmethod = "none", to not filter any ligand before.
# The other lrmatchingmethods are problematic, since "expression" takes into account proportions
                                        # and "distribution" does not work as expected.
# id_r_list is a list of identity/conserved receptors per population
###########
getIdentitySignals <- function(data,unique_cts,id_r_list, removeCompetingReceptors = FALSE, lrmatchingmethod = "expression"){
    ## =====================================
    ## I add to account for global variables
    data_r <- data[rownames(data) %in% RLIST,]
    data_l <- data[rownames(data) %in% LLIST,]
    ## =====================================
    lr_pairs <- list()
    for(i in seq_along(unique_cts)){
        cur_r <- c()
        cur_l <- c()
        curSubPop <- unique_cts[i]
        for(r in id_r_list[[i]]){
            ## cells <- which(colnames(data) %in% rownames(cts_in_data)[which(cts_in_data$cell.type == curSubPop)])
            cells <- colnames(data) == curSubPop
            receptor_expr_r <- sum(data_r[r,cells]) ## not used at all for method "none"
            ligands_r <-LRDAT$Ligand[which(LRDAT$Receptor == r)]
            ligands_r <- ligands_r[which(ligands_r %in% rownames(data_l))]
            if(lrmatchingmethod == "none"){
                idx <- ligands_r
            }
### Removed other lrmatchingmethod's here
### ....
            cur_r <- c(cur_r,rep(r,length(idx)))
            cur_l <- c(cur_l,idx)
        }
        lr_pairs[[i]] <- data.frame(Ligand = cur_l, Receptor = cur_r, stringsAsFactors = FALSE)
    }
    return(lr_pairs)
}
