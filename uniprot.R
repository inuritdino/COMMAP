## Add UniProt DB interface to select further the receptors

uniprot.filter.receptors <- function(rec.list){
    suppressPackageStartupMessages({
        if(ORG == "hsapiens")
            require(org.Hs.eg.db)
        if(ORG == "mmusculus")
            require(org.Mm.eg.db)
        require(UniProt.ws)
    })
    if(ORG == "hsapiens"){
        up <- UniProt.ws(taxId=9606)
        ens.sym.df <- select(org.Hs.eg.db,rec.list,"ENSEMBL","SYMBOL")
    }
    if(ORG == "mmusculus"){
        up <- UniProt.ws(taxId=10090)
        ens.sym.df <- select(org.Mm.eg.db,rec.list,"ENSEMBL","SYMBOL")
    }
    locs <- select(up,ens.sym.df$ENSEMBL,"SUBCELLULAR-LOCATIONS","ENSEMBL")
    ens.sym.loc.df <- merge(locs,ens.sym.df)
    idx <- grep("membrane",ens.sym.loc.df$`SUBCELLULAR-LOCATIONS`)
    rec.list.new <- unique(ens.sym.loc.df$SYMBOL[idx])
    return(rec.list.new)
}
