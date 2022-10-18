#!/bin/python3
import re
import sys,os
import random

def Genes_generate():
    genes = []
    [genes.append("OR" + str(i)) for i in range(1,23)]
    [genes.append("Tas2r" + str(i)) for i in range(1,21)]
    [genes.append("NC" + str(i)) for i in range(1,19)]
    return genes

def Gene_Match_InGff(gff):
    genes = Genes_generate()
    Match_gene = {}
    Match_gene_one = {}
    with open(gff,"r") as f:
        for line in f:
            cons = line.strip().split("\t")
            if cons[2] != "gene":continue
            Gene_Symbol = cons[8].split("=")[2]
            #print(Gene_Symbol)
            for i in genes:
                if re.match(i,Gene_Symbol):
                    gene = Gene_Symbol
                    if gene not in Match_gene:
                         Match_gene[gene] = [[cons[0],cons[3],cons[4]]]
                    else:
                        Match_gene[gene].append([cons[0],cons[3],cons[4]])
                else:
                    continue
    #print(Match_gene)
    #### 若基因在gff中多次出现，则取第一个
    for i in Match_gene:
        if len(Match_gene[i]) > 1:
            Match_gene_one[i] = Match_gene[i][0]
        else:
            Match_gene_one[i] = [Match_gene[i][0][0],Match_gene[i][0][1],Match_gene[i][0][2]]
    #print(Match_gene_one)
    return Match_gene_one

def Coding_fst():
    Match_gene = Gene_Match_InGff(gff)
    for i in Match_gene:
        CHR = re.search(r'[0-9]+',Match_gene[i][0])[0]
        regions = CHR + ":" +  Match_gene[i][1] + "-" + Match_gene[i][2]
        out = i + "_" + Match_gene[i][0] + "_" + Match_gene[i][1] + "-" + Match_gene[i][2]
        
        start = Match_gene[i][1]
        end = Match_gene[i][2]
        out = "./fst/" + i
        out_ = out + ".recode.vcf"
        fst = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {} && vcftools --vcf {} --fst-window-size 50 --fst-window-step 20 --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
        print(fst.format(CHR,start,end,out,out_,out))
        
        #### calculate Fst
        #step1 = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {}"
        #print(step1.format(CHR,start,end,out))
        #step2 = "vcftools {} --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
        #print(step2.format(out_,out))

        # pick up the region of genes
        #python gene.py > genes_regions.sh
        #print("bcftools filter M.vcf.gz --regions {} > ./genes/{}".format(regions,out))
        #python gene.py > genes_regions.sh
        #python gene.py > genes_regions.fst.sh
        #fst = "vcftools --vcf ./genes/{} --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out ./genes/{}" + ".fst"
        #print(fst.format(out,out))

def Noncoding_fst():
    Match_gene = Gene_Match_InGff(gff)
    Non_code_res = Non_code(gff)
    for i in Match_gene:
        CHR = Match_gene[i][0]
        start = Match_gene[i][1]
        end = Match_gene[i][2]
        length = int(end) - int(start)

        tmp = []
        for j in Non_code_res:
            if int(j[1]) > int(start) + 200000 and int(j[2]) - int(j[1]) > length:
                #print(i,j[1],j[2])
                tmp.append((j[0],j[1],j[2]))
        rand = random.randint(0,len(tmp))
        Non = tmp[rand]
        #print(i,Non)
        #### check!!!!right!!!!
        #out = "GeneName:" + i + "_" + "GeneCHR:" + CHR + "_"  + "NoncodeCHR:" + Non[0]  +  "_" + "GeneLength: " + str(length) + " " + "START: " + start + " " + "END: " + end + "_" + "NoncodeLength:" + str(int(Non[1]) + length - int(Non[1]))
        #print(out)
        out = "/fst/" + i + "." + "Backgrond"
        out_ = out + ".recode.vcf"
        fst = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {} && vcftools --vcf {} --fst-window-size 50 --fst-window-step 20 --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
        print(fst.format(Non[0],Non[1],str(int(Non[1]) + length),out,out_,out))

        #### calculate Fst
        #step1 = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {}" 
        #print(step1.format(Non[0],Non[1],str(int(Non[1]) + length),out))
        #step2 = "vcftools {} --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
        #print(step2.format(out_,out))

def Non_code(gff):
    None_codes = []
    Non_code_res = []
    with open(gff,"r") as f:
        for line in f:
            cons = line.strip().split("\t")
            if cons[2] != "gene":continue
            tmp = (cons[3],cons[4],cons[0])
            None_codes.append(tmp)
    #print(None_codes)
    #i = 0
    start = "0"
    end = "0"
    for j in range(len(None_codes)):
        if j == 0: 
            #print("The None code regions from %s to %s" %(0,None_codes[j][0]))
            #Non_code_res.append((0,None_codes[j][0]))
            continue
        start = None_codes[j-1][1]
        end = None_codes[j][0]
        CHRID = re.search(r'[0-9]+',None_codes[j][2])[0]
        #print("The None code regions from %s to %s" %(start,end))
        #print("\t".join([str(CHRID),start,end]))
        Non_code_res.append((str(CHRID),start,end))
    #print(Non_code_res)
    return Non_code_res

if __name__ == "__main__":
    gff = "03.addnozhushi.final.genesymbol.gff"
    #Gene_Match_InGff(gff)
    Coding_fst()
    #### python genes.py > None_code.txt
    #Non_code(gff)
    Noncoding_fst()   #python genes.py > Gene_Gackground.fst.sh
