# R 4.0.4-sdc
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)



tmp=read.table("~/cyt/out/inferCNV13/infercnv.references.txt", header=T)
cnv_table <- read.table("~/cyt/out/inferCNV13/infercnv.observations.txt", header=T)

all_cnv_table <- cbind(tmp,cnv_table)

down=mean(rowMeans(tmp)) - 2 * mean( apply(tmp, 1, sd))  #0.948
up=mean(rowMeans(tmp)) + 2 * mean( apply(tmp, 1, sd))  #1.0569
oneCopy=up-down
oneCopy  #0.1085806
a1= down- 2*oneCopy
a2= down- 1*oneCopy
down;up
a3= up +  1*oneCopy
a4= up + 2*oneCopy 

# Replicate the table 
all_cnv_score_table <- as.matrix(all_cnv_table)

all_cnv_score_mat <- as.matrix(all_cnv_table)
# Scoring
all_cnv_score_table[all_cnv_score_mat > 0 & all_cnv_score_mat < a2] <- "A" #complete loss. 2pts
all_cnv_score_table[all_cnv_score_mat >= a2 & all_cnv_score_mat < down] <- "B" #loss of one copy. 1pts
all_cnv_score_table[all_cnv_score_mat >= down & all_cnv_score_mat <  up ] <- "C" #Neutral. 0pts
all_cnv_score_table[all_cnv_score_mat >= up  & all_cnv_score_mat <= a3] <- "D" #addition of one copy. 1pts
all_cnv_score_table[all_cnv_score_mat > a3  & all_cnv_score_mat <= a4 ] <- "E" #addition of two copies. 2pts
all_cnv_score_table[all_cnv_score_mat > a4] <- "F" #addition of more than two copies. 2pts

# Check
table(all_cnv_score_table[,2])
# Replace with score 
all_cnv_score_table_pts <- all_cnv_table
rm(cnv_score_mat)
# 
all_cnv_score_table_pts[all_cnv_score_table == "A"] <- 2
all_cnv_score_table_pts[all_cnv_score_table == "B"] <- 1
all_cnv_score_table_pts[all_cnv_score_table == "C"] <- 0
all_cnv_score_table_pts[all_cnv_score_table == "D"] <- 1
all_cnv_score_table_pts[all_cnv_score_table == "E"] <- 2
all_cnv_score_table_pts[all_cnv_score_table == "F"] <- 2

# Scores are stored in “cnv_score_table_pts”. Use colSums to add up scores for each cell and store as vector 
all_cell_scores_CNV <- as.data.frame(colSums(all_cnv_score_table_pts))
colnames(all_cell_scores_CNV) <- "cnv_score"
head(all_cell_scores_CNV)
write.csv(x = all_cell_scores_CNV, file = "~/cyt/rscript/epi-one-by-one-inferCNV/all_cnv_scores.csv")

pbmc_harmony2 <- Hepatocyte
phe=pbmc_harmony2@meta.data
head(rownames(phe))
head(rownames(all_cell_scores_CNV)) 

rownames(all_cell_scores_CNV)=gsub('^X','',rownames(all_cell_scores_CNV))
rownames(all_cell_scores_CNV)=gsub('[.]','-',rownames(all_cell_scores_CNV))
head(rownames(all_cell_scores_CNV))

phe=phe[rownames(phe) %in% rownames(all_cell_scores_CNV),]  
phe$cnv_scores  =  all_cell_scores_CNV[rownames(phe),]

head(rownames(phe))
dim(phe) 

colnames(phe)
library(ggpubr)
p2=ggboxplot(phe,'group3','cnv_scores', fill = "group3")
p2 


group_name <- c('N1','N2','T1-1','T1-2','T2-1','T2-2')
phe$group3 <- factor(phe$group3,levels = group_name)

p2 <- ggplot(phe,aes(x = group3, y = cnv_scores,fill = group3)) + geom_violin()
p2
p2=p2+theme(legend.position="none",
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text( size=rel(1)),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size=1)
)
p2

p2=p2+geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")
p2

p2=p2+theme_bw()
p2 = p2+NoLegend()
p2

my_comparisons2 <- list( c('N1',"N2"), c("T1-1", "T1-2"), c("T1-2", "T2-1"),c("T2-1","T2-2") )
o <-  ggplot(phe,aes(x=group3,y=cnv_scores,fill=group3))+geom_violin()+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size=rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1)
  )+
 # scale_fill_manual(values = color)+
  stat_compare_means( comparisons = my_comparisons2)+
  stat_compare_means(label.y = 15)+
  ggtitle('CNV score between groups')
o




##
phe$group5[phe$group3 == 'N1'] <- 'Normal hepatocytes'
phe$group5[phe$group3 == 'N2'] <- 'Normal hepatocytes'
#phe$group5[phe$group3 == 'CSC-like cells'] <- 'CSC-like cells'
phe$group5[phe$group3 == 'T1-1'] <- 'Malignant cells'
phe$group5[phe$group3 == 'T1-2'] <- 'Malignant cells'
phe$group5[phe$group3 == 'T2-1'] <- 'Malignant cells'
phe$group5[phe$group3 == 'T2-2'] <- 'Malignant cells'
table(phe$group3,phe$group5)

group_name <- c('Normal hepatocytes','Malignant cells')
phe$group5 <- factor(phe$group5,levels = group_name)


p1=ggboxplot(phe,'group5','cnv_scores', fill = "group5")
p1 

p <- ggplot(phe,aes(x = group5, y = cnv_scores,fill = group5)) + geom_violin()
p
p=p+theme(legend.position="none",
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text( size=rel(1)),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size=1)
)
p

p=p+geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")
p

p=p+theme_bw()
p


library(ggplot2)
library(ggpubr)
color = c("#33B44A","#EE3536")
#color = c("#33B44A","#EE3536","#3A429B")
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

my_comparisons <- list( c("Normal hepatocytes", "Malignant cells") )
y <-  ggplot(phe,aes(x=group5,y=cnv_scores,fill=group5))+geom_violin()+
  geom_boxplot(width=0.2,position=position_dodge(0.9),outlier.colour = NA,fill="white")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5),
        axis.text.y = element_text(hjust = 0.5, vjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text( size=rel(1)),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black",size=1)
  )+
  scale_fill_manual(values = color)+
  stat_compare_means( comparisons = my_comparisons,hide.ns = T,label = "p.signif",method = 'wilcox.test',label.x.npc=5,)+
  stat_compare_means(label.y = 15)+
  ggtitle('CNV score between groups')
y




phe2 <- phe
phe2$cellID <- rownames(phe2)
phe2 <- phe2[,c('cellID','cnv_scores')]
pbmc_harmony3 <- pbmc_harmony2  #肝细胞
#pbmc_harmony3@meta.data <- merge(pbmc_harmony3@meta.data ,phe2,by = 'cellID',all.x = T )
pbmc_harmony3 <- AddMetaData(pbmc_harmony3,phe2)
head(pbmc_harmony3)
x1 <- FeaturePlot(pbmc_harmony3,features = 'cnv_scores',order = T,cols = c('#440256','#2A768E','#48C16E','#F8E621'))
x1
x2 <- FeaturePlot(pbmc_harmony3,features = 'cnv_scores',split.by = 'group3',order = T,cols = c('#440256','#2A768E','#48C16E','#F8E621'),ncol = 3) / DimPlot(pbmc_harmony3,split.by = 'group3',ncol = 3)
x2
#pbmc_harmony <- AddMetaData(pbmc_harmony,pbmc_harmony3@meta.data)
pbmc_harmony4 <- pbmc_harmony
pbmc_harmony4 <- AddMetaData(pbmc_harmony4,pbmc_harmony3@meta.data)
head(pbmc_harmony4)
x3 <- FeaturePlot(pbmc_harmony4,features = 'cnv_scores',order = T,cols = c('#440256','#2A768E','#48C16E','#F8E621'))
x3
DimPlot(pbmc_harmony4,split.by = 'group3')
x1|x2
#saveRDS(pbmc_harmony3,"~/cyt/rscript/inferCNV/Hepatocyte_CSC_like_cells_cnv_score.rds")

z <- DimPlot(pbmc_harmony3)+ggtitle('CellType')
z
i <- DimPlot(pbmc_harmony3,group.by = 'group3')
i
c <- DimPlot(pbmc_harmony)
c
(x|z|i)/(c|y)
dev.off()

pdf('~/cyt/out/cnv_score/Hepatocyte_CSC_like_cells_cnv_score.pdf',width = 11,height = 8)
(x|z|i)/(c|y)
dev.off()











value <- pbmc_harmony3@meta.data$cnv_scores
thresthold <- median(value) - 1* sd(value) #941 这里是一个标准差
thresthold <- quantile(value,probs = 0.75)- 1* sd(value)
thresthold <- quantile(value,probs = 0.75)- 1* sd(value)
thresthold <- max(pbmc_harmony3@meta.data$cnv_scores[pbmc_harmony3$group == 'Normal']) #3191
thresthold <- median(pbmc_harmony3@meta.data$cnv_scores[pbmc_harmony3$group == 'Normal']) #705
thresthold <- mean(pbmc_harmony3@meta.data$cnv_scores[pbmc_harmony3$group == 'Normal']) #746.7
thresthold <- 1100

pbmc_harmony3@meta.data %>% mutate(condition=if_else(.$cnv_scores <= thresthold,'non-malignant','malignant')) -> pbmc_harmony3@meta.data # 分出了恶性和非恶性
head(pbmc_harmony3)
table(pbmc_harmony3@meta.data$condition)
table(pbmc_harmony3@meta.data$condition,pbmc_harmony3@meta.data$group3)
p1 <- DimPlot(object = pbmc_harmony3,label=F, pt.size=0.5,group.by='condition',cols = c('#E41A1C','#00C4BF')) # 画图
p1
p2 <- DimPlot(object = pbmc_harmony3,label=F, pt.size=0.5,split.by='condition') # 画图
p2
p3 <- DimPlot(pbmc_harmony3,split.by = 'group3')
(p1|x3)/ p3/(y|p2)



#cyt:左边这个图可以做成右边这种形式吗？就是在原来的细胞注释umap上显示出恶性细胞和非恶性细胞。其它细胞细胞类型标灰。
{
  pbmc_harmony5 <- pbmc_harmony
  pbmc_harmony5 <- AddMetaData(pbmc_harmony5,pbmc_harmony3@meta.data)
  head(pbmc_harmony5)
  xxx <- DimPlot(pbmc_harmony5,group.by = 'condition')
  xxx
  ggsave("/share/home/lx/cyt/out/cnv_score/malignant_cells_umap.pdf", xxx, width = 6, height = 5)
}


#开始绘图
library(Seurat)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(scales)
pbmc <- pbmc_harmony3
colourCount = length(unique(pbmc@meta.data$condition))
library(colorspace)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

celltype_colors <- getPalette(colourCount)

#提取样本和细胞数据，并且进行长宽数据转换
pbmc2@meta.data <- pbmc2@meta.data[pbmc2@meta.data$group3 %in% selected_group,]
table(pbmc2@meta.data$group3)
plotC <- table(pbmc2@meta.data$group3, pbmc2@meta.data$condition) %>% melt()

colnames(plotC) <- c("Sample", "condition","Number")

#绘制每个组织中细胞数目柱状图

pC1 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = condition)) +
  
  geom_bar(stat = "identity", width=0.8,aes(group=condition),position="stack")+
  
  scale_fill_manual(values=celltype_colors) +
  
  theme_bw()+
  
  theme(panel.grid =element_blank()) +
  
  labs(x="",y="Average number")+
  
  theme(axis.text = element_text(size=12, colour = "black"))+
  
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))

#绘制每个组织中细胞比例柱状图

pC2 <- ggplot(data = plotC, aes(x = Sample, y = Number, fill = condition)) +
  
  geom_bar(stat = "identity", width=0.8,aes(group=condition),position="fill")+
  
  scale_fill_manual(values=celltype_colors) +
  
  theme_bw()+
  
  theme(panel.grid =element_blank()) +
  
  labs(x="",y="Cell proportion")+
  
  #scale_y_continuous(labels = percent)+ ####用来将y轴移动位置
  
  theme(axis.text = element_text(size=12, colour = "black"))+
  
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8, vjust = 0.6))#让横轴上的标签倾斜45度

#两个图片进行拼图

pC <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')
pC
pC1+pC2+p
pp <- pC1 + pC2 + plot_layout(ncol = 2, widths = c(1,1),guides = 'collect')+DimPlot(pbmc_harmony)+p +DimPlot(pbmc_harmony,group.by = 'group')
pp


ppp <- (x3|xxx|p1)/ p3/(y|p2)/pC +plot_annotation(tag_levels = "A")
ppp
ggsave("/share/home/lx/cyt/out/cnv_score/malignant_cells_by_cnv_scores.pdf", ppp, width = 13, height = 18)
ggsave("/share/home/lx/cyt/out/cnv_score/malignant_cells_by_cnv_scores.png", ppp,width = 13, height = 15)
#ggsave("/share/home/lx/cyt/out/copykat/pred_mallignant2.pdf", pp, width = 10, height = 11)

mallignant <- read.delim("../rscript/copykat/Hepatocytes_sub_copykat_prediction.txt", row.names = 1)
pbmc_harmony3 <- AddMetaData(pbmc_harmony3, metadata = mallignant)
head(pbmc_harmony3)

yyy <- DimPlot(pbmc_harmony3, group.by = "copykat.pred") + scale_color_manual(values = c("red", "gray",'blue'))

ppp <- (x3|xxx|p1|yyy)/ p3/(y|p2)/pC +plot_annotation(tag_levels = "A")
ppp
ggsave("/share/home/lx/cyt/out/cnv_score/malignant_cells_by_cnv_scores.pdf", ppp, width = 16, height = 21)
ggsave("/share/home/lx/cyt/out/cnv_score/malignant_cells_by_cnv_scores.png", ppp,width = 16, height = 22)











