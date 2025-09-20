
# Plot circRNAs

# Old function to plot reads supporting a circRNA
circ_readsPlot <- function(circ, bamfiles, exongtf, maxreadstoplot=50, seed=NULL){
  set.seed(seed)
  
  # Check and process exon info
  stopifnot('GRanges' %in% class(exongtf))
  exns <- subset(exongtf,type == 'exon')
  if (length(exns) <1) {
    stop("exongtf does not contain exon information")
  }
  
  # Get chimeric reads
  reads <- getchimreads(circ,bamfiles)
  # Subsample to avoid overplotting
  keepreads <- sample(
    x=unique(names(reads)),
    size= min(maxreadstoplot, length(unique(names(reads)))),
    replace=F)
  reads <- reads[names(reads) %in% keepreads]
  
  # Tranform to a datatable for plotting
  df <- data.table(as.data.frame(reads,row.names = NULL,use.outer.mcols=T))
  
  # Invert coordinates for - strand 
  df[,start_d := start]
  df[,end_d := end]
  df[strand == '-',end_d := start]
  df[strand == '-',start_d := end]
  
  # Assign an index to each read pair (ordered by mate1 start)
  idx <- df[,.(st=min(start[mate == '1'])),by=qname]
  idx[order(st),idx:=1:.N]
  df <- merge(df,idx,by='qname')
  
  # Exon scheme
  dfex <- exonTrack(circ,exongtf)
  # Plot 
  ggplot(NULL)+
    geom_segment(data=df,
                 aes(x=start_d, xend=end_d, y=idx+2,yend=idx+2,col=as.factor(mate)),
                 arrow=arrow(length = unit(3,'points'),type = 'closed'))+
    geom_vline(xintercept=unique(c(dfex$xmin,dfex$xmax)),linetype=2,linewidth=0.25)+
    geom_hline(yintercept=0.9,col='gray20')+
    geom_rect(data=dfex,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),col='gray20',fill="gray20",linewidth=0.25)+
    scale_color_manual(name = "Mate",values = c("#4DACD6","#C02D45"))+
    scale_y_continuous(breaks=(1:max(df$idx))+2)+
    coord_cartesian(xlim=c(start(circ)-10,end(circ)+10))+
    theme_classic()+ xlab('Genomic position') + ylab(NULL) +
    theme(
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank())+
    ggtitle(sprintf("%s - %s (%s)",circ$gene_id,circ$circID, as.character(circ)))
}



# Low level function to extract chimeric reads spanning a circRNA backsplice junction
.getchimreads <- function(circ, bamfile){
  # Check objects are as expected
  stopifnot('GRanges' %in% class(circ) & length(circ) == 1)
  
  # Check bam file is sorted
  stopifnot(file.exists(bamfile))
  h <- scanBamHeader(bamfile)
  if (!"SO:coordinate" %in% h[[1]][['text']][['@HD']]){
    stop("bamfile does not seem to be sorted")
  }
  
  # Load alignments
  reads <- readGAlignments(
    bamfile,
    param = ScanBamParam(
      which=circ,
      what=c('qname','flag'),
      tag='ch'
    )
  )
  names(reads) <- mcols(reads)$qname
  
  # Keep pairs where there are chimeric alignments
  reads <- reads[names(reads) %in% unique(names(reads[!is.na(mcols(reads)$ch)]))]
  # Note that this would be equivalent, since the whole pair is tag as ch 
  # reads <- reads[!is.na(mcols(reads)$ch)]
  
  # Add mate info
  mcols(reads)$mate <- bamFlagTest(mcols(reads)$flag, 'isFirstMateRead')+1
  
  # Extract cigar matches in reference space coordinates
  ivs <- grglist(reads, use.names=T,use.mcols=T)
  
  # Filters 
  # Keep reads that span the junction 
  keepreads <- union(
    findOverlaps(ivs,resize(circ,1,'start'),maxgap=1,ignore.strand=T),
    findOverlaps(ivs,resize(circ,1,'end'),maxgap=1,ignore.strand=T))
  ivs <- ivs[names(ivs) %in% unique(names(ivs)[from(keepreads)])]
  
  # Keep pairs that are completely within the circRNA
  ivsl <- unlist(ivs)
  ivsl <-  unlist(range(split(ivsl, names(ivsl))))
  keepreads <- findOverlaps(ivsl,circ,type='within',ignore.strand=T)
  ivs <- ivs[names(ivs) %in% unique(names(ivsl)[from(keepreads)])]
  
  # keep only complete pairs
  keepreads <- mcols(ivs)
  keepreads <- with(keepreads, aggregate(flag,list(qname),function(x) length(unique(x)) > 2))
  keepreads <- with(keepreads, Group.1[x])  
  
  # Apply filters
  ivs[names(ivs) %in% keepreads]
}

# Function to extract chimeric reads spanning circRNA BSJ from one or multiple bam files
getchimreads <- function(circ,bamfiles){
  do.call(
    'c',
    lapply(bamfiles, .getchimreads,circ=circ)
  ) 
}

# Low level function to generate coverage of chimeric reads spanning a circRNA BSJ
.genBSJCov <- function(circ,bamfile, norm_fact=1, label='label'){
  # Get chimeric reads
  reads <- getchimreads(circ,bamfile)
  cov <- coverage(reads)[circ][[1]]
  cov <- data.frame(
    x=(1:length(cov))+start(circ),
    y=decode(cov)*norm_fact,
    label=label,
    row.names = NULL
  )
  cov
}


# Function to generate a plottable track with the coverage of chimeric reads around a BSJ
circ_covTrack <- function(circ, bamfiles,norm_facts=1,labels='trackname',makeplot=F){
  df <- data.frame(
    'f' = bamfiles,
    'n' = norm_facts,
    'label'=labels
  )
  covs <- apply(df,1, function(x) .genBSJCov(circ,x[[1]],as.numeric(x[[2]]),x[3]))
  cov <- do.call('rbind',covs)
  if (makeplot){
    ggplot(cov, aes(x=x,y=y))+
      stat_summary(fun='sum',geom='area',aes(fill=label),show.legend = F)+
      facet_wrap(~label,strip.position = 'right',ncol=1)+
      theme_bw()+ xlab(NULL)+ylab('Coverage')+
      theme(panel.grid=element_blank()
      )
  } else {
    cov
  }
}

# Function to generate a plottable track with the exon structure of a circRNA
exonTrack <- function(circ, exongtf, makeplot=FALSE){
  # Check objects are as expected
  stopifnot('GRanges' %in% class(circ) & length(circ) == 1)
  stopifnot('GRanges' %in% class(exongtf))
  # Get exons (unique)
  exns <- subset(exongtf,type == 'exon')
  if (length(exns) <1) {
    stop("exongtf does not contain exon information")
  }
  # Exon scheme
  dfex <- unique(subsetByOverlaps(exns, circ,ignore.strand=T))
  dfex <- data.frame(
    'xmin'= start(dfex),
    'xmax'= end(dfex),
    'ymin'= 0,
    'ymax'= 2,
    'exon_id'= dfex$exon_id
  )
  if (makeplot){
    return(
      ggplot(dfex, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax))+
        geom_hline(yintercept=1)+
        geom_rect(fill='gray80',col='gray20',linewidth=0.25)+
        theme_classic()+ 
        theme(
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          axis.text= element_blank())
    )}
  else {
    return(dfex)
  }
}


# Load exons
mygtf <- import('/data/genomes/GRCm39_GENCODE_All/gencode.vM32.chr_patch_hapl_scaff.annotation.gtf',format='gtf')
mygtf <- subset(mygtf, type=='exon')


# Load circs
mycircs <- fread('Data/circRNA_LAKI/All_circRNAs_annotated.tsv')
mycircs <- makeGRangesFromDataFrame(mycircs,keep.extra.columns = T)
circ <- mycircs[mycircs$circID == 'cRNA_0854'] # circInsr

# Load samples info
sampinfo <- fread('Data/Info/Sample_info_complete.tsv')
sampinfo[,files:=file.path('/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/',SampleID,'Aligned.out.sorted.bam')]


# Previous function to plot reads
circ_readsPlot(circ,'/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/Tube_106/Aligned.out.sorted.bam',mygtf)
# Now it accepts multiple files
circ_readsPlot(circ, sampinfo[Tissue == 'Heart' & Genotype == 'KO'][['files']],mygtf,seed=234)
ggsave("Out/LAKI/circRNA/circInsr_reads.pdf", height = 3, width = 5)

# New function to plot coverage of chimeric reads
# requires a vector of files
# An optional vector of normalization factors can be provided (to normalized by library reads, for example)
# An optional vector of the labels of each bamfile
# Example:
# p1 <- circ_covTrack(circ = circ,
#                     bamfiles = sampinfo[Tissue == 'Heart'][['files']],
#                     norm_facts = 1,
#                     labels = sampinfo[Tissue == 'Heart'][['Genotype']],
#                     makeplot = T)
# p1  
# 
# # Example introducing normalization factors (simulated to reduce coverage to 0.5 fold)
# sampinfo[,norm_fact:=rnorm(nrow(sampinfo),0.5,0.1)]
# p1 <- circ_covTrack(circ = circ,
#                     bamfiles = sampinfo[Tissue == 'Heart'][['files']],
#                     norm_facts = sampinfo[Tissue == 'Heart'][['norm_fact']],
#                     labels = sampinfo[Tissue == 'Heart'][['Genotype']],
#                     makeplot = T)
# p1  
# 
# 
# # New function to generate a exon structure track
# p2 <- exonTrack(circ,mygtf,makeplot=T)
# p2
# 
# # Now, we can combine the tracks 
# p<- wrap_plots(p1,p2,heights = c(10,1),ncol=1)
# ggsave('circRNA_analysis/cov_example_plot.pdf',p)

# # Load exons
# mygtf <- import('/data/genomes/GRCm39_GENCODE_All/gencode.vM32.chr_patch_hapl_scaff.annotation.gtf',format='gtf')
# mygtf <- subset(mygtf, type=='exon')
# 
# 
# mybam <- '/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/Tube_106/Aligned.out.sorted.bam'
# mybamWT <- '/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/Tube_112/Aligned.out.sorted.bam'
# 
# pdf("Out/LAKI/circRNA/plot_BSJ.pdf", width = 8, height = 4)
# 
# # circBfar
# circ_readsPlot(mycircs[mycircs$circID == 'cRNA_0234'],
#               sampinfo[Tissue == 'Heart' & Genotype == 'KO'][['files']],
#               mygtf)
# 
# # circArhgap5
# circ_readsPlot(mycircs[mycircs$circID == 'cRNA_2382'],
#               sampinfo[Tissue == 'Heart' & Genotype == 'KO'][['files']],
#               mygtf)
# 
# # circInsr
# circ_readsPlot(mycircs[mycircs$circID == 'cRNA_0854'],
#               sampinfo[Tissue == 'Heart' & Genotype == 'KO'][['files']],
#               mygtf)
# 
# # circR3hcc1l
# circ_readsPlot(mycircs[mycircs$circID == 'cRNA_1421'],
#               sampinfo[Tissue == 'Heart' & Genotype == 'WT'][['files']],
#               mygtf)
# 
# dev.off()
