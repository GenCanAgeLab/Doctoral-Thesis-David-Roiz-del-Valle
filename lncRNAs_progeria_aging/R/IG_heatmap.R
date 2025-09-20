


# Load gencode

gcode <- readRDS('Data/geneMetadata.rds')
IG_code <- subset(gcode, grepl('^IG',gene_type))
IG_genes <- gcode$gene_id2[grepl('^IG',gcode$gene_type)]                
H_genes <- gcode$gene_id2[grepl('^H2-',gcode$gene_name)]
Other_genes <- gcode$gene_id2[gcode$gene_name %in% c('Jchain','C4a','C4b')]
goi <- c(IG_genes,H_genes,Other_genes)

## With Fold Change

tms <- readRDS('Out/TMS/DESeq_results_TMS_oldvsyoung.rds')
# goi <- unique(unlist(lapply(tms, function(x) rownames(subset(x,padj < 0.05 & log2FoldChange > 0)))))
tms <- lapply(tms,function(x) {
  x <- x[rownames(x) %in% goi,]
  x$gene_id <- rownames(x)
  x[,c('baseMean','log2FoldChange','padj','gene_id')]})
tmsnames <- rep(names(tms),sapply(tms,nrow))
tms <- do.call('rbind',tms)
tms <- cbind(tms, matrix(unlist(strsplit(tmsnames,'\\.')),ncol=2,byrow = T,dimnames = list(NULL,c('DS','Tissue'))))
# tms <- left_join(as.data.frame(tms),gcode,join_by(gene_id == gene_id2))
tms$DS <- paste0('TMS_',tms$DS)
# tms <- subset(tms, sex=='male')



zhou <- readRDS('Out/Zhou/DESeq_results_Zhou_oldvsyoung.rds')
zhou <- lapply(zhou,function(x) {
  x <- x[rownames(x) %in% goi,]
  x$gene_id <- rownames(x)
  x[,c('baseMean','log2FoldChange','padj','gene_id')]})
zhounames <- rep(names(zhou),sapply(zhou,nrow))
zhou <- do.call('rbind',zhou)
zhou$DS <- 'Zhou'
zhou$Tissue <- zhounames


laki <- readRDS('Out/LAKI/DESeq_results_LAKI.rds')
laki <- lapply(laki,function(x) {
  rownames(x) <- gsub('\\.[0-9]+','',rownames(x))
  x <- x[rownames(x) %in% goi,]
  x$gene_id <- rownames(x)
  x[,c('baseMean','log2FoldChange','padj','gene_id')]})
lakinames <- rep(names(laki),sapply(laki,nrow))
laki <- do.call('rbind',laki)
laki$DS <- 'laki'
laki$Tissue <- lakinames

df <- as.data.frame(rbind(laki,zhou,tms))
signgoi <- unique(subset(df,padj < 0.05 & abs(log2FoldChange) > 1)$gene_id)
df <- subset(df,gene_id %in% signgoi)

# Inject missing values (not needed, it is easier to set a gray background or a color scale without middle white)
# df <- tidyr::expand_grid(unique(df[,c('DS','Tissue')]),'gene_id'=signgoi)|>
#   left_join(df,join_by(DS,Tissue,gene_id))


df <- left_join(as.data.frame(df),gcode,join_by(gene_id == gene_id2))

df$gene_type <- gsub('_gene','',df$gene_type)
df$gene_type[!grepl('^IG',df$gene_type)] <- 'Others'


# Remove IGkv
df <- subset(df, gene_name != 'Igkv4-77')
p <- ggplot(df,aes(x=Tissue,y=gene_name,fill=log2FoldChange))+
  geom_tile()+theme_bw()+
  facet_grid(rows=vars(gene_type),cols=vars(DS),scale='free',space='free')+
  scale_fill_gradient2('Log2FC',low='blue',mid='white',high='red',midpoint = 0)+
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5),
        panel.spacing.y = unit(0.5,'point'),
        axis.text.y=element_text(size = 6))
p
ggsave('Out/IG_genes_heatmap.pdf',p,width=7,height = 12)

p <- ggplot(df,aes(y=Tissue,x=gene_name,fill=log2FoldChange))+
  geom_tile()+theme_bw()+
  facet_grid(cols=vars(gene_type),rows=vars(DS),scale='free',space='free')+
  scale_fill_gradient2('Log2FC',low='blue',mid='white',high='red',midpoint = 0)+
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=0.5),
        panel.spacing.y = unit(0.5,'point'),
        axis.text.y=element_text(size=8))
  p
ggsave('Out/IG_genes_heatmap2.pdf',p,width=18,height = 6)

lims <- range(df$log2FoldChange)
p <- ggplot(df,aes(y=Tissue,x=gene_name,fill=log2FoldChange))+
  geom_tile(na.rm=F)+
  facet_grid(cols=vars(gene_type),rows=vars(DS),scale='free',space='free')+
  scale_fill_gradientn(values=scales::rescale(c(lims[1],seq(-3,3,length.out=99),lims[2])),
                       colours = c('green','white','red'))+
  xlab(NULL)+ylab(NULL)+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=6,angle=90,hjust=1,vjust=0.5),
        panel.spacing.y = unit(5,'point'),
        axis.text.y=element_text(size=8))
p
ggsave('Out/IG_genes_heatmap3.pdf',p,width=18,height = 6)



# Plot baseMean of IGKC

ggplot(subset(df,gene_name == 'Igkc'), aes(x=Tissue,y=baseMean))+
  geom_bar(stat='identity')+
  facet_grid(cols=vars(DS),scale='free_x',space='free_x')+
  xlab(NULL)+ylab('Igkc mean levels (counts per million)')+theme_bw()+
  theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
