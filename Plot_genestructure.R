# =======================
# 1. 加载包
# =======================
library(ggplot2)
library(dplyr)
library(gridExtra)

# =======================
# 2. 读取数据
# =======================
gene_info <- read.table("gene_info.txt", header=TRUE, sep="\t")
gene_family <- read.table("gene_family.txt", header=FALSE, sep="\t")
fam_class <- read.table("family.freq_class.txt", header=FALSE, sep="\t")

colnames(gene_family) <- c("Family","geneID")
colnames(fam_class) <- c("Family","Size","Class")

# 避免 factor 问题
gene_info$geneID <- as.character(gene_info$geneID)
gene_family$geneID <- as.character(gene_family$geneID)
gene_family$Family <- as.character(gene_family$Family)
fam_class$Family <- as.character(fam_class$Family)

# =======================
# 3. 合并
# =======================
df <- gene_family %>%
  inner_join(fam_class, by="Family") %>%
  inner_join(gene_info, by="geneID")

df$Class <- factor(df$Class, levels=c("Core","Softcore","Shell","Cloud"))

# =======================
# 4. log10 转换
# =======================
df <- df %>%
  mutate(
    logGene   = log10(GeneLen + 1),
    logCDS    = log10(CDSLen + 1),
    logExon   = log10(ExonLen + 1),
    logIntron = log10(IntronLen + 1)
  )

# =======================
# 5. 颜色
# =======================
mycol <- c(
  "Core"="#3d5887",
  "Softcore"="#a1bcd9",
  "Shell"="#86984a",
  "Cloud"="#d69c26"
)

# =======================
# 6. Wilcoxon函数
# =======================
get_p_label <- function(data, var, g1, g2){
  x <- data[data$Class==g1, var]
  y <- data[data$Class==g2, var]
  p <- wilcox.test(x, y)$p.value
  
  if(p < 0.001) return("***")
  else if(p < 0.01) return("**")
  else if(p < 0.05) return("*")
  else return("ns")
}

# =======================
# 7. 绘图函数
# =======================
plot_func <- function(var, ylab){

  labels <- c(
    get_p_label(df, var, "Core","Softcore"),
    get_p_label(df, var, "Core","Shell"),
    get_p_label(df, var, "Core","Cloud")
  )

  p <- ggplot(df, aes(x=Class, y=.data[[var]], color=Class, fill=Class)) +

    geom_boxplot(width=0.5, outlier.shape=NA, alpha=0.3) +

    geom_jitter(width=0.25, size=1.5, alpha=0.7) +

    scale_color_manual(values=mycol) +
    scale_fill_manual(values=mycol) +

    theme_bw() +
    labs(x=NULL, y=ylab) +

    theme(
      legend.position="none",
      axis.text=element_text(size=11),
      axis.title=element_text(size=12)
    )

  # 显著性
  y_max <- max(df[[var]], na.rm=TRUE)

  for(i in 1:3){
    p <- p +
      annotate("segment",
               x=1, xend=i+1,
               y=y_max + i*0.3,
               yend=y_max + i*0.3) +
      annotate("text",
               x=(1+i+1)/2,
               y=y_max + i*0.3 + 0.1,
               label=labels[i],
               size=4)
  }

  return(p)
}

# =======================
# 8. 四张图
# =======================
p1 <- plot_func("logGene",   "log10(Gene Length)")
p2 <- plot_func("logCDS",    "log10(CDS Length)")
p3 <- plot_func("logExon",   "log10(Exon Length)")
p4 <- plot_func("logIntron", "log10(Intron Length)")

# =======================
# 9. 输出（2x2排版）
# =======================
pdf("Comparison_4panels.pdf", width=10, height=8)

grid.arrange(
  p1, p2, p3, p4,
  ncol=2
)

dev.off()
