#!/bin/R

args <- commandArgs(T)
back <- paste0(args[1],".back")
gene <- paste0(args[1],".gene")
setwd("/data/00/user/user196/01.liangxl/wanna/Olfactory_Receptor_Analysis/File")
backfst <- read.table(back)
genefst <- read.table(gene)
# t-test
T_test <- t.test(backfst,genefst)
P_value <- T_test$p.value

Fst_mean_gene <- mean(genefst$V1)
Fst_mean_back <- mean(backfst$V1)

# HEADER
# cat("Genes","Fst_mean_gene","Fst_mean_back","P_value")
cat(args[1],Fst_mean_gene,Fst_mean_back,P_value,"\n")
