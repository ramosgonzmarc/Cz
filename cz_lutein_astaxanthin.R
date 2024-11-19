## Gene expression
normalized.gene.expression <- read.table(file="gene_expression.tsv",
                              header = T,as.is=T,sep="\t")
head(normalized.gene.expression)

gene.ids <- normalized.gene.expression$geneID
normalized.gene.expression <- normalized.gene.expression[,2:13]
rownames(normalized.gene.expression) <- gene.ids
head(normalized.gene.expression)

norm.exp.red <- data.frame(numeric(length = nrow(normalized.gene.expression)))
names.red <- c()

for (i in seq(1, ncol(normalized.gene.expression), by = 3))
{
  new_vec <- rowMeans(normalized.gene.expression[,c(i, i+1, i+2)])
  norm.exp.red <- cbind(norm.exp.red, new_vec)
  new_name <- colnames(normalized.gene.expression)[i]
  names.red <- c(names.red, new_name)
}

norm.exp.red <- norm.exp.red[,-1]
names.red <- sapply(strsplit(names.red, split = "3_1"), function(x) x[1])
names.red <- paste0(names.red, "3")
colnames(norm.exp.red) <- names.red
head(norm.exp.red)

# Carotenoids content
area_carotenoids_percentage <- read.table(file="carotenoids_area.tsv",sep="\t",header = T,as.is = T)
head(area_carotenoids_percentage)

car.red <- data.frame(numeric(length = 19))
names.car <- c()

for (i in seq(1, ncol(area_carotenoids_percentage), by = 4))
{
  new_vec <- rowMeans(area_carotenoids_percentage[,c(i, i+1, i+2, i+3)])
  car.red <- cbind(car.red, new_vec)
  new_name <- colnames(area_carotenoids_percentage)[i]
  names.car <- c(names.car, new_name)
}

car.red <- car.red[,-1]
names.car <- sapply(strsplit(names.car, split = "3_1"), function(x) x[1])
names.car <- paste0(names.car, "3")
colnames(car.red) <- names.car
head(car.red)

# Correlation between important carotenoid-related enzymes and carotenoids levels
enzymes.names <- c("CYP97A1","CYP97A2","CYP97C","CHYb", "CYP4F", "BKT1", "BKT2")
enzymes.ids <- c("Cz13g16110", "Cz09g14130", "Cz09g07100","Cz12g16080", "Cz12g27180", "Cz13g13100", "Cz04g11250")

enzymes.red <- norm.exp.red[enzymes.ids,]
merged_whole <- t(rbind(enzymes.red, car.red))

asoc_whole <- cor(merged_whole)[1:7,-c(1,2,3,4,5,6,7)]
rownames(asoc_whole) <- enzymes.names
asoc_whole

# Plots
# Carotenoids
barplot.heatmap.carotenoids <- function(car.name,normalized.car.content, max.fc)
{
  expression.hl.1mM <- as.numeric(normalized.car.content[car.name, paste("HL_1mMNO3", 1:4,sep="_")])
  expression.hl.10mM <- as.numeric(normalized.car.content[car.name, paste("HL_10mMNO3", 1:4,sep="_")])
  expression.ll.1mM <- as.numeric(normalized.car.content[car.name, paste("LL_1mMNO3", 1:4,sep="_")])
  expression.ll.10mM <- as.numeric(normalized.car.content[car.name, paste("LL_10mMNO3", 1:4,sep="_")])
  
  means <- c(mean(expression.hl.1mM),
             mean(expression.hl.10mM),
             mean(expression.ll.1mM),
             mean(expression.ll.10mM))
  
  sds <- c(sd(expression.hl.1mM),
           sd(expression.hl.10mM),
           sd(expression.ll.1mM),
           sd(expression.ll.10mM))
  
  png(filename = paste0("figs_carotenoids/barplot_", car.name,".png"),
      width = 400)
  par(lwd=3,mar=c(7,6,2,2))
  xpos <- barplot(means,col=c("coral4","coral2",
                              "darkseagreen3","darkseagreen4"),
                  names.arg = c("HL 1mM", "HL 10 mM",
                                "LL 1mM", "LL 10 mM"),
                  las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=car.name,
                  cex.main=2)
  
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  dev.off()
  
  # color palette blue to red
  rwb <- colorRampPalette(colors = c("blue", "white", "red"))
  range.colors <- rwb(max.fc*200)
  
  log2.fcs <- 100*log2(c(means[1]/means[4],
                         means[2]/means[4],
                         means[3]/means[4]))
  
  log2.fcs <- log2.fcs + max.fc*100
  log2.fcs <- round(log2.fcs)
  log2.fcs[log2.fcs > max.fc*200] <- max.fc*200
  log2.fcs[log2.fcs < 1] <- 1
  
  png(filename = paste0("figs_carotenoids/heatmap_", car.name, ".png"),
      width = 1000,height=300)
  plot(1, 1, col = "white", axes=F,xlab="",ylab="",xlim=c(-1,31),ylim=c(-1,11))  
  polygon(x = c(0, 0, 10, 10),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][1],lwd=6)
  polygon(x = c(10, 10, 20, 20),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][2],lwd=6)
  polygon(x = c(20, 20, 30, 30),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][3],lwd=6)
  dev.off()
  
  
  return(list(HL_1=means[1], HL_10=means[2],
              LL_1=means[3], LL_10=means[4],
              fc_HL1=means[1]/means[4],
              
              fc_HL10=means[2]/means[4],
              
              fc_LL5=means[3]/means[4]))
}

# Before running this, create a folder named figs_carotenoids in this directory 
for (x in rownames(area_carotenoids_percentage))
{
  barplot.heatmap.carotenoids(car.name = x, 
                              normalized.car.content = area_carotenoids_percentage, 
                              max.fc = 6)
}

# Carotenoids-related enzymes
library(data.table)

other.enzymes <- fread("carotenoids_gene_table.csv", fill = T, select = c(1,2), header = T,
                       sep="\t")
other.enzymes.data <-  as.data.frame(other.enzymes)[-c(53,54,55,56),]

barplot.heatmap.gene.other <- function(gene.id,gene.name,normalized.gene.expression,max.fc)
{
  expression.hl.1mM <- 2^unlist(normalized.gene.expression[gene.id,
                                                           paste("HL_1mMNO3", 1:3,sep="_")]) - 1
  expression.hl.10mM <- 2^unlist(normalized.gene.expression[gene.id,
                                                            paste("HL_10mMNO3", 1:3,sep="_")]) - 1
  
  expression.ll.1mM <- 2^unlist(normalized.gene.expression[gene.id,
                                                           paste("LL_1mMNO3", 1:3,sep="_")]) - 1
  expression.ll.10mM <- 2^unlist(normalized.gene.expression[gene.id,
                                                            paste("LL_10mMNO3", 1:3,sep="_")]) - 1
  
  means <- c(mean(expression.hl.1mM),
             
             mean(expression.hl.10mM),
             mean(expression.ll.1mM),
             
             mean(expression.ll.10mM))
  
  sds <- c(sd(expression.hl.1mM),
           
           sd(expression.hl.10mM),
           sd(expression.ll.1mM),
           
           sd(expression.ll.10mM))
  
  png(filename = paste(paste0("figs_other_enzymes/barplot_", gene.name),paste0(gene.id,".png"),sep="_"),
      width = 400)
  par(lwd=3,mar=c(7,6,2,2))
  xpos <- barplot(means,col=c("coral4","coral2",
                              "darkseagreen3","darkseagreen4"),
                  names.arg = c("HL 1mM", "HL 10 mM",
                                "LL 1mM", "LL 10 mM"),
                  las=2,cex.names = 1.5,
                  ylim=c(0,max(means+sds)*1.2),cex.axis = 1.5,lwd=3,
                  main=paste(c(gene.name, "-", gene.id),collapse=" "),
                  cex.main=2)
  mtext(side = 2,text = "FPKM",line = 4,cex = 2)
  arrows(x0 = xpos,y0 = means + sds, x1 = xpos, y1 = means - sds,length = 0.1,
         code = 3,angle=90,lwd=2)
  dev.off()
  
  # color palette blue to red
  rwb <- colorRampPalette(colors = c("blue", "white", "red"))
  range.colors <- rwb(max.fc*200)
  
  log2.fcs <- 100*log2(c(means[1]/means[4],
                         means[2]/means[4],
                         means[3]/means[4]))
  
  log2.fcs <- log2.fcs + max.fc*100
  log2.fcs <- round(log2.fcs)
  log2.fcs[log2.fcs > max.fc*200] <- max.fc*200
  log2.fcs[log2.fcs < 1] <- 1
  
  png(filename = paste(paste0("figs_other_enzymes/heatmap_", gene.name),paste0(gene.id,".png"),sep="_"),
      width = 1000,height=300)
  plot(1, 1, col = "white", axes=F,xlab="",ylab="",xlim=c(-1,31),ylim=c(-1,11))  
  polygon(x = c(0, 0, 10, 10),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][1],lwd=6)
  polygon(x = c(10, 10, 20, 20),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][2],lwd=6)
  polygon(x = c(20, 20, 30, 30),                           
          y = c(0, 10, 10, 0),                           
          col = range.colors[log2.fcs][3],lwd=6)
  dev.off()
  
  return(list(HL_1=means[1], HL_10=means[2],
              LL_1=means[3], LL_10=means[4],
              fc_HL1=means[1]/means[4],
              
              fc_HL10=means[2]/means[4],
              
              fc_LL5=means[3]/means[4]))
}

# Before running this, create a folder named figs_other_enzymes in this directory
for (i in 1:nrow(other.enzymes.data))
{
  barplot.heatmap.gene.other(gene.id = other.enzymes.data$ID[i], gene.name = other.enzymes.data$`Gene name`[i],
                       normalized.gene.expression = normalized.gene.expression,
                       max.fc = 3)
}

# Two way Anova for all genes
# Function for computing the test for a given gene
two.way.anova <- function(gene.id)
{
  
  new.table <- subset(anova.data, gene == gene.id)[,-1]
  mod <- aov(FPKM ~ nitrate * light,
             data = new.table)
  
  # Check assumptions
  if (shapiro.test(mod$residuals)$p.value < 0.05 || car::leveneTest(mod)$`Pr(>F)`[1] < 0.05)
  {
    assumptions <- F
  }
  else
  {
    assumptions <- T
  }
  
  return(data.frame(nitrate = summary(mod)[[1]][["Pr(>F)"]][1],
                    light = summary(mod)[[1]][["Pr(>F)"]][2],
                    interaction = summary(mod)[[1]][["Pr(>F)"]][3],
                    assumptions = assumptions))
  
}

# Create data table
total.ids <- c("Cz13g13100", "Cz04g11250", "Cz13g16110", "Cz09g14130", "Cz09g07100","Cz12g16080")
total.names <- c("BKT1", "BKT2","CYP97A1","CYP97A2","CYP97C","CHYb")

library(data.table)

other.enzymes <- fread("carotenoids_gene_table.csv", fill = T, select = c(1,2), header = T,
                       sep="\t")
other.enzymes.data <- as.data.frame(other.enzymes)[-c(53,54,55,56),]

complete.genes <- data.frame(id = c(total.ids,other.enzymes.data$ID), name = c(total.names,other.enzymes.data$`Gene name`))
complete.genes <- complete.genes[!duplicated(complete.genes),]
dim(complete.genes)                             

total.red.table <- normalized.gene.expression[complete.genes$id,]
head(total.red.table)

anova.data <- data.frame(gene=rep(complete.genes$id, each = 12), light = rep(c("high", "low"),each = 6), 
                         nitrate = rep(c("1","1","1","10","10","10"),2), FPKM=as.numeric(t(total.red.table)))
head(anova.data)

# Apply function
anova.res <- sapply(unique(anova.data$gene), function(x) two.way.anova(x))
anova.res.t <- t(anova.res)
anova.res <- cbind(anova.res.t, name = complete.genes$name)
anova.res <- as.data.frame(anova.res)

# For those without significant interaction, repeat test
rep.ids <- which(anova.res$interaction > 0.05)

for (x in names(rep.ids))
{
  new.table <- subset(anova.data, gene == x)[,-1]
  mod <- aov(FPKM ~ nitrate + light,
             data = new.table)
  anova.res[x, "nitrate"] <- summary(mod)[[1]][["Pr(>F)"]][1]
  anova.res[x, "light"] <- summary(mod)[[1]][["Pr(>F)"]][2]
    
  
}

library(data.table)

anova.dt <- as.data.table(anova.res)
anova.sig <- anova.dt[, nit_sig := ifelse(nitrate < 0.05, "*", "-")] 
anova.sig <- anova.sig[, light_sig := ifelse(light < 0.05, "*", "-")]
anova.sig <- anova.sig[, int_sig := ifelse(interaction < 0.05, "*", "-")]
anova.sig <- anova.sig[, ids := rownames(anova.res)]
fwrite(anova.sig, "anova_res.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


# Two way Anova for carotenoids
anova.data.car <- data.frame(gene=rep(rownames(area_carotenoids_percentage), each = 16), light = rep(c("high", "low"),each = 8), 
                         nitrate = rep(c("1","1","1","1","10","10","10","10"),2), content=as.numeric(t(area_carotenoids_percentage)))
head(anova.data.car)

# Function for carotenoids
two.way.anova.car <- function(gene.id, anova.data)
{
  
  new.table <- subset(anova.data, gene == gene.id)[,-1]
  mod <- aov(content ~ nitrate * light,
             data = new.table)
  
  # Check assumptions
  if (ks.test(mod$residuals, "pnorm", mean = mean(mod$residuals), sd = sd(mod$residuals))$p.value < 0.05 || car::leveneTest(mod)$`Pr(>F)`[1] < 0.05)
  {
    assumptions <- F
  }
  else
  {
    assumptions <- T
  }
  
  return(data.frame(nitrate = summary(mod)[[1]][["Pr(>F)"]][1],
                    light = summary(mod)[[1]][["Pr(>F)"]][2],
                    interaction = summary(mod)[[1]][["Pr(>F)"]][3],
                    assumptions = assumptions))
  
}


# Apply function
anova.res.car <- sapply(unique(anova.data.car$gene), function(x) two.way.anova.car(x, anova.data = anova.data.car))
anova.res.car.t <- t(anova.res.car)
anova.res.car <- as.data.frame(cbind(anova.res.car.t, name = rownames(anova.res.car.t)))

# For those without significant interaction, repeat test
rep.ids <- which(anova.res.car$interaction > 0.05)

for (x in names(rep.ids))
{
  new.table <- subset(anova.data.car, gene == x)[,-1]
  mod <- aov(content ~ nitrate + light,
             data = new.table)
  anova.res.car[x, "nitrate"] <- summary(mod)[[1]][["Pr(>F)"]][1]
  anova.res.car[x, "light"] <- summary(mod)[[1]][["Pr(>F)"]][2]
  
}

library(data.table)

anova.dt.car <- as.data.table(anova.res.car)
anova.sig <- anova.dt.car[, nit_sig := ifelse(nitrate < 0.05, "*", "-")] 
anova.sig <- anova.sig[, light_sig := ifelse(light < 0.05, "*", "-")]
anova.sig <- anova.sig[, int_sig := ifelse(interaction < 0.05, "*", "-")]
fwrite(anova.sig, "anova_res_car.tsv", sep = "\t", quote = F, col.names = T, row.names = F)
