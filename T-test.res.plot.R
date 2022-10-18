library(tidyverse)
library(ggrepel)

setwd("./")
dat <- read.table("T-test.res",header = FALSE)
names(dat) <- c("Genes","Fst_mean_gene","Fst_mean_back","P_value")

dat_tmp <- dat %>% gather(key,value,-Genes,-P_value)

#### significant gene
SigGene <- subset(dat_tmp,P_value < 0.05)

p <- ggplot(dat_tmp,aes(x=key,y=value,color=key)) + geom_point(size=1) +
  geom_jitter(alpha=0.4) +
  #geom_text(data = subset(dat_tmp,P_value < 3.063396e-13 & key == "Fst_mean_gene"),mapping = aes(label=Genes)) +
  geom_text_repel(data = subset(dat_tmp,P_value < 8.060141e-06 & key == "Fst_mean_gene"),mapping = aes(x=key,y=value,label=Genes)) +
  labs(x="",y="Fst")

ggsave(p,filename = "wanna.pdf")
write.table(SigGene,file = "wanna.txt",quote = FALSE,row.names = FALSE)
