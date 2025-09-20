
# Import Circ annotations
  circ_annotation <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/circRNA_ID_across_tools.tsv")
  circmeta <- read_tsv("/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Results/All_circRNAs_annotated_cross_with_circDBs.tsv")


# Load BSJ data
  ciri2Res <- read_tsv("Data/circRNA_LAKI/Ciri2_circRNA_counts.tsv") |> 
    mutate(
      Tool = 'ciri2',
      ID = circ_annotation$circID[match(ID,circ_annotation$CIRI2)]) |>
    pivot_longer(cols=starts_with('Tube'))
  
  circExp2Res <- read_tsv("Data/circRNA_LAKI/circExplorer2_circRNA_counts.tsv") |> 
    rename_all(function(x) sub('_circRNA_known.txt','',x)) |>
    mutate(
      Tool = 'circExp',
      ID = circ_annotation$circID[match(ID,circ_annotation$circExplorer2)])|>
    pivot_longer(cols=starts_with('Tube'))
  
  findCircRes <- read_tsv("Data/circRNA_LAKI/findcirc_circRNA_counts.tsv") |> 
    mutate(
      Tool = 'find_circ',
      ID = circ_annotation$circID[match(ID,circ_annotation$find_circ)])|>
    pivot_longer(cols=starts_with('Tube'))
  
  starchipRes <- read_tsv("Data/circRNA_LAKI/starchip_circRNA_counts.tsv") |> 
    mutate(
      Tool = 'starchip',
      ID = circ_annotation$circID[match(ID,circ_annotation$starchip)])|>
    pivot_longer(cols=starts_with('Tube'))


# Read lineal junctions
  myfiles <- list.files('/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/Lineal_splicing/',
                        full.names=T)
  
  linjunc <- read_tsv(myfiles,id = 'name') |>
    mutate(name=gsub('_lineal_junctions.tsv','',basename(name)))



# Read sample info
  sampdf <- read.table(
    "Data/Info/Sample_info_complete_genotiped.tsv",
    header = T,
    sep = "\t") |> 
    filter(Notes != "Ratones de Face-NR") |>
    # Fix wrong genotypes and remove HT animals (see Bash/Regenotyping)
    dplyr::rename(Regenotype = Infered_Genotype) |>
    filter(Regenotype != "HT") |> 
    mutate(Regenotype = if_else(Regenotype == "WT","WT","LAKI")) |>
    dplyr::mutate(Genotype = Regenotype) |> 
    #Tube_100 and Tube_20 are incorrectly classified as males,
    #and Tube_91 and Tube_11 incorrectly as females. So we change them
    mutate(Gender = case_when(
      SampleID %in% c("Tube_100","Tube_20") ~ "Female",
      SampleID %in% c("Tube_91","Tube_11") ~ "Male",
      .default = Gender
    ))
  
  
# Calculate LSJ, filter by more than 10 BSJ, calculate BSJ ratio (just for circRNA exons) 
  plotdt <- rbind(ciri2Res,circExp2Res,findCircRes,starchipRes) |> 
      left_join(linjunc,join_by(ID == circID,name == name)) |>
      inner_join(circmeta, join_by(ID == circID)) |> 
      filter(value >= 10) |>
      filter(type %in% c('JNC-JNC','FULLEX')) |>
      rowwise() |>
      mutate(TSJ=value+sum(Njnc5,Njnc3,na.rm=T)/2) |>
      mutate(BSJratio=value/TSJ*100) |>
      inner_join(sampdf, join_by( name == SampleID)) |>
      group_by(Tissue,Genotype,Gender,Tool,name) |>
      summarise(MeanBSJratio = mean(BSJratio,na.rm=T)) |>
    mutate(Genotype = factor(Genotype, levels=c('WT','LAKI')))
  
  
  # Plot BSJratio
  p <- ggplot(plotdt, aes(x=Tissue,y=MeanBSJratio,fill=Genotype)) +
    geom_boxplot(outliers = F,position=position_dodge(0.75)) +
    geom_jitter(position=position_jitterdodge(0.4,0,0.75),
                size=0.5,
                shape = 16) +
    facet_wrap(~Tool, ncol = 2) +
    ylab('Mean BSJ ratio (%)') +
    theme_classic() +
    # theme(
    #   axis.text.x=element_text(angle=45,hjust=1)
    # ) +
    scale_fill_manual(values = c("#0072B2", "#ff3b30")) +
    ggpubr::stat_compare_means(
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")),
      aes(label = after_stat(p.signif)),
      hide.ns = TRUE,
      size = 8,
      vjust = 1
    )


ggsave('Out/LAKI/circRNA/BSJ2LSJratio_of_circRNA_exons.pdf',p,width=7,height = 5)

  
  
######### NUMBER of BSJ normalized and unnormalized #############
  
# Number of BSJs (Raw values without normalizing)
  plotdt <- rbind(ciri2Res,circExp2Res,findCircRes,starchipRes) |> 
    inner_join(sampdf, join_by( name == SampleID)) |>
    mutate(Genotype = factor(Genotype, levels=c('WT','LAKI'))) |>
    group_by(Tissue,Genotype,Gender,Tool,name)|>
    summarize(BSJ=sum(value))
    
  ggplot(plotdt, aes(x=Tissue,y=BSJ,fill=Genotype))+
    # stat_summary(fun=mean, geom='bar',position=position_dodge(0.8),width=0.8)+
    geom_boxplot(outliers=F,position=position_dodge(0.8),width=0.8)+
    geom_jitter(aes(group=Genotype),
                position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),
                size=0.7)+
    facet_grid(~Tool)+
    theme_bw()
  
  # Read ALL linear splicing for normalization (It doesn't make much sense)
  # lsjfiles <- list.files('/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/',
  #        pattern='SJ.out.tab',
  #        full.names=T,
  #        recursive = T)
  # lsj <- sapply(lsjfiles, function(x) sum(read_tsv(x,col_names = F,show_col_types = F)[['X7']]))  
  # lsj <- data.frame('name'=basename(dirname(names(lsj))),LSJ=unname(lsj))
  # 
  # plotdt <- plotdt |>
  #   inner_join(lsj,join_by(name))
  
  # ggplot(plotdt, aes(x=Tissue,y=BSJ/(BSJ+LSJ),fill=Genotype))+
  #   # stat_summary(fun=mean, geom='bar',position=position_dodge(0.8),width=0.8)+
  #   geom_boxplot(outliers=F,position=position_dodge(0.8),width=0.8)+
  #   geom_jitter(aes(group=Genotype),
  #               position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),
  #               size=0.7)+
  #   facet_grid(~Tool)+
  #   theme_bw()

  # Normalize by total number of uniquely mapped reads
  genefiles <- list.files('/home/jandrulas/20231017_RNASeq_LAKI/circRNA_analysis/STAR_chim_aligns/',
                         pattern='ReadsPerGene.out.tab',
                         full.names=T,
                         recursive = T)
  treads <- sapply(genefiles, function(x) sum(read_tsv(x,col_names = F,show_col_types = F,skip=4)[['X4']]))  
  treads <- data.frame('name'=basename(dirname(names(treads))),Treads=unname(treads))
  
  plotdt <- plotdt |>
    inner_join(treads,join_by(name))
  
  p <- ggplot(plotdt, aes(x=Tissue,y=BSJ/Treads*100,fill=Genotype))+
    # stat_summary(fun=mean, geom='bar',position=position_dodge(0.8),width=0.8)+
    geom_boxplot(outliers=F,position=position_dodge(0.8),width=0.8)+
    geom_jitter(aes(group=Genotype),
                position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),
                size=0.5,
                shape = 16)+
    facet_wrap(~Tool, ncol = 2)+ ylab('% of BSJ reads')+
    theme_classic()+
    # theme(
    #   axis.text.x=element_text(angle=45,hjust=1)
    # ) +
    scale_fill_manual(values = c("#0072B2", "#ff3b30")) +
    ggpubr::stat_compare_means(
      symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")),
      aes(label = after_stat(p.signif)),
      hide.ns = TRUE,
      size = 8,
      vjust = 1
    )
  ggsave('Out/LAKI/circRNA/BSJ_number_total_reads_normalized.pdf',p,width=7,height = 5)
  

##### Table of BSJ ratio to use for filtering ##########  
  
  
  
# Join (summarizing BSJ by max betweTreads# Join (summarizing BSJ by max between tools)
allres <- rbind(ciri2Res,circExp2Res,findCircRes,starchipRes) |> 
  group_by(ID,name) |> 
  summarize(BSJ=max(value,na.rm=T)) |>
  left_join(linjunc,join_by(ID == circID,name == name)) |>
  rowwise() |>
  mutate(LSJ = sum(Njnc5,Njnc3,na.rm = T)) |>
  mutate(BSJratio = BSJ/(BSJ+LSJ)) |>
  inner_join(sampdf, join_by(name == SampleID)) |>
  group_by(ID,Tissue,Gender, Genotype) |>
  summarize(
    BSJratio=mean(BSJratio,na.rm=T),
    BSJ = mean(BSJ,na.rm=T),
    LSJ = mean(LSJ,na.rm=T)) |>
  inner_join(circmeta,join_by(ID == circID))

b<-allres |>
  filter(BSJ >= 10) |>
  group_by(Tissue,ID,Gender) |>
  mutate(BSJrFold=BSJratio/mean(BSJratio[Genotype=='WT'],na.rm=T))


# Generate a bed track for UCSC
bedtrack <- rbind(ciri2Res,circExp2Res,findCircRes,starchipRes) |> 
  group_by(ID,name) |> 
  summarize(BSJ=max(value,na.rm=T)) |>
  group_by(ID) |>
  summarize(score=sum(BSJ))|>
  inner_join(circmeta,join_by(ID == circID)) |>
  dplyr::select(seqnames,start,end, ID,score,strand) |>
  dplyr::rename(name = ID) |>
  mutate(score = round(score*1000/max(score,na.rm=T)))
bedtrack <- makeGRangesFromDataFrame(bedtrack,keep.extra.columns = T)
rtracklayer::export(bedtrack,'Out/LAKI/circRNA/circRNA_UCSC_bedtrack.bed',format='bed',
                    trackLine=as("track type=bed name=LAKI_circRNAs useScore=1",'TrackLine'))

                    