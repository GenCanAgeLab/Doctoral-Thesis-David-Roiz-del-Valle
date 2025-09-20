
require(Rsamtools)
require(GenomicRanges)
require(rtracklayer)


.find_inforf <- function(seq){
  # Find frames without stop codons
  seq <- rep(seq,3)
  e <- gregexpr('TAA|TAG|TGA',as.character(seq))[[1]]
  inf_frames <- setdiff(0:2,e %% 3)
  if (length(inf_frames) > 0){
    s <- matchPattern('ATG',seq)
    s <- s[start(s) %% 3 %in% inf_frames]
    s <- s[start(s) < length(seq)/3]
    ir <- Views(seq,IRanges(start(s),start(s) + (length(seq)/3) -1))
    res <- as(ir,'DNAStringSet')
    res <- translate(res)
    return (data.frame(
      'type' = 'Inf_ORF',
      'name'=paste0('ORF_',start(ir)),
      'frame'= (start(ir) %% 3) + 1,
      'start'=start(ir),
      'length'=width(ir)/3,
      'seq'= as.character(res))
    )
  } else {
    return(NULL)
  }
}
  

.find_finorfs <- function(seq,f=1,orfmin=25){
  seq3 <- rep(seq,3)
  suppressWarnings(b<- translate(subseq(seq3,start = f),no.init.codon = T))
  orfs <- matchLRPatterns('M','*',Inf,b)
  orfs <- orfs[start(orfs) < (length(b) / 3)]
  orfs <- unlist(endoapply(split(ranges(orfs), start(orfs)), function(x) x[which.min(width(x))]))
  orfs <- unlist(endoapply(split(ranges(orfs), end(orfs)), function(x) x[which.max(width(x))]))
  orfs <- orfs[width(orfs) >= orfmin]
  s <- (start(orfs)*3)+f-3
  if (length(orfs)==0){return(NULL)}
  data.frame(
    'type'= 'Def_ORF',
    'name' = paste0('ORF_',s),
    'frame' = f,
    'start' = s,
    'length'= width(orfs)-1,
    'seq' = as.character(AAStringSet(Views(b,orfs))),
    row.names = NULL
  )
  }

find_orfs <- function(seq,orfmin=25){
  finorfs <- do.call('rbind',lapply(1:3, .find_finorfs, seq=seq,orfmin=orfmin))
  inforfs <- .find_inforf(seq)
  rbind(finorfs,inforfs)
}


pcr_temp <- function(dnastr){
  i <- floor(length(dnastr)/2)
  h5 <- subseq(dnastr,start=2,end=i)
  mid <- c(dnastr[length(dnastr)],dnastr[1])
  h3 <- subseq(dnastr,start=i+1, end=-2)
  paste0(
    as.character(h3)
    ,'[',
    as.character(mid),
    ']',
    as.character(h5))
}


getcircseq <- function(circgr, exongrs,fafile, returnexns=F,fixcoords=F){
  # Try to find one exon
  exin <- sort(unique(subsetByOverlaps(exongrs,circgr,type='equal',maxgap=4)))
  # Try to find several exons
  if (length(exin) == 0){
    sex <- subsetByOverlaps(exongrs,circgr,type='start',maxgap=2)
    eex <- subsetByOverlaps(exongrs,circgr,type='end',maxgap=2)
    txs <- intersect(sex$transcript_id,eex$transcript_id)
    exin <- subsetByOverlaps(exongrs,resize(circgr,width(circgr)+8,fix = 'center'),type='within')
    exin <- exin[exin$transcript_id %in% txs]
    exin <- sort(unique(exin))
  }
  if (length(exin) == 0){
    stop("circRNA does not match exon junctions")
  }
  edif <- end(exin[length(exin)]) - end(circgr)
  sdif <- start(exin[1]) - start(circgr)
  if (sdif != 0 | edif !=0){
    warning(sprintf("circ start and edn do not match an exon start. sdif = %s, edif = %s", sdif,edif))
    if (fixcoords == T) {
      start(circgr) <- start(circgr) + sdif
      end(circgr) <- end(circgr) + edif
      warning(sprintf("Exon coords are %s",as.character(circgr)))
      warning(sprintf("New circ coord are %s",as.character(exin)))
    }
  }
  exin <- subsetByOverlaps(exin, circgr, type='within')
  # endexns <- c(exin[1],exin[length(exin)])
  # extraex <- exin[countOverlaps(exin,endexns) == 0]
  # exin <- sort(c(endexns,extraex))
  if (!isDisjoint(exin)){
    warning("exons are not disjoint")
  }
  s <- unlist(getSeq(fafile, unstrand(exin)))
  if (unique(strand(circgr)) == '-'){
    s <- reverseComplement(s)
  }
  if (returnexns) {return(exin)}
  return(s)
  }


myex <- import('/data/genomes/GRCm39_GENCODE_All/gencode.vM32.chr_patch_hapl_scaff.annotation.gtf')
myex <- subset(myex, type == 'exon')

mycircs <- read.delim('/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/All_circRNAs_annotated_cross_with_circDBs.tsv')
mycircs <- makeGRangesFromDataFrame(mycircs,keep.extra.columns = T)

fa <- FaFile('/data/genomes/GRCm39_GENCODE_All/GRCm39.genome.fa')

(s <- getcircseq(mycircs[mycircs$circID == 'cRNA_1690'],myex,fa))
(o <- find_orfs(s))

# Ankrd17 seq does not match with circatlas


cRNA_0887

flag <- 'DYKDDDDK'
flag <- 'EQKLISEEDL' # HA
flag <- unlist(strsplit(flag,split=''))
flagd <- lapply(flag, function(x) names(GENETIC_CODE[GENETIC_CODE == x]))
seqs <- DNAStringSet(apply(expand.grid(flagd),1,paste0,collapse=''))
s <- seqs[vcountPattern('AGTG',seqs)>0]

# Original myc
# GAA CAA AAA CTC ATC TCA GAA GAG GAT CTG
# Splicing version (change ser TCA to AGT) They have similar codon usage
# GAA CAA AAA CTC ATC AGT GAA GAG GAT CTG
