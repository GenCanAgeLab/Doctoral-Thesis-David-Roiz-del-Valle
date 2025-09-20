
require(GenomicRanges)
require(rtracklayer)
require(data.table)


circs <- fread('circRNA_analysis/Results/All_circRNAs_annotated.tsv')
circs <- makeGRangesFromDataFrame(circs,keep.extra.columns = T)

getlinjnc <- function(jncfile, circgrs){
  # Read lineal junctions
  jnc <- fread(jncfile,
               col.names = c('seqnames','start','end','strandN','InMot','Intype','Nuni','Nmul','Over'))
  jnc <- makeGRangesFromDataFrame(jnc,keep.extra.columns = T,starts.in.df.are.0based = F)
  start(jnc) <- start(jnc)-1 # because are intron coords
  end(jnc) <- end(jnc) + 1 # because are intron coords
  # Find overlaps (upstream and downstream)
  lin5 <- findOverlaps(resize(unstrand(circs),1,'end'),resize(jnc,1,'start'),ignore.strand=T,maxgap = 2)
  lin3 <- findOverlaps(resize(unstrand(circs),1,'start'),resize(jnc,1,'end'),ignore.strand=T,maxgap = 2)
  # Extract the N unique reads mapping the junction (field 7 in Sj.out.tab)
  lin5 <- sum(splitAsList(jnc$Nuni[to(lin5)],from(lin5)))
  lin3 <- sum(splitAsList(jnc$Nuni[to(lin3)],from(lin3)))
  # Arrange results
  data.frame(
    circID = circs$circID,
    'Njnc5' = lin5[as.character(1:length(circs))],
    'Njnc3' = lin3[as.character(1:length(circs))]
  ) 
}

myfiles <- list.files('circRNA_analysis/STAR_chim_aligns',recursive = T,pattern = 'SJ.out.tab',full.names=T)

lapply(myfiles, function(x){
  res <- getlinjnc(x,circs)
  outn <- paste0(basename(dirname(x)), "_lineal_junctions.tsv")
  write.table(res,
              file.path('circRNA_analysis/Lineal_splicing/',outn),
              quote=F, row.names = F, sep='\t')
})



myfiles <- list.files('circRNA_analysis/Lineal_splicing/', full.names=T)
mydata <- lapply(myfiles,fread)
mydata <- sapply(mydata, function(x) rowSums(x[,2:3],na.rm=T))
circs2 <- cbind(as.data.frame(circs),Njunc=rowSums(mydata))
