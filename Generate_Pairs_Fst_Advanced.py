import re
import sys,os
import random
import time

#### Input
#workdir = "C:/Users/xllia/Desktop"
#os.chdir(workdir)  #修改路径
#gff = "addnozhushi.final.genesymbol.gff"
gff = sys.argv[1]

#### create gene
def Generate_gene():
    genes = []
    [genes.append("OR" + str(i)) for i in range(1, 23)]
    [genes.append("Tas2r" + str(i)) for i in range(1, 21)]
    [genes.append("NC" + str(i)) for i in range(1, 19)]
    return genes

genes = Generate_gene()
genes = tuple(genes)

#### gene mapping in gff file
def Gene_Match_InGff(gff):
    #genes = Genes_generate()
    Match_gene = {}
    Match_gene_one = {}
    with open(gff,"r") as f:
        for line in f:
            cons = line.strip().split("\t")
            if cons[2] != "gene":continue
            Gene_Symbol = cons[8].split("=")[2]
            for i in genes:
                if re.match(i,Gene_Symbol):
                    gene = Gene_Symbol
                    if gene not in Match_gene:
                         Match_gene[gene] = [[cons[0],cons[3],cons[4]]]
                    else:
                        Match_gene[gene].append([cons[0],cons[3],cons[4]])
                else:
                    continue
    #### 若基因在gff中多次出现，则取第一个
    for i in Match_gene:
        if len(Match_gene[i]) > 1:
            Match_gene_one[i] = Match_gene[i][0]
        else:
            Match_gene_one[i] = [Match_gene[i][0][0],Match_gene[i][0][1],Match_gene[i][0][2]]
    return Match_gene_one

Match_genes = Gene_Match_InGff(gff)
#print(Match_genes)

#### generate Non-coding region
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
    return Non_code_res
Non_code = Non_code(gff)

class Fst:
    def __init__(self,gene,chr,start,end,length,out):
        self.gene = gene
        self.chr = chr
        self.start = start
        self.end = end
        self.length = length
        self.out = out

    def fst_gene(self):
        fst = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {} && vcftools --vcf {} --fst-window-size 50 --fst-window-step 20 --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
        gene_fst = fst.format(self.chr,self.start,self.end,self.out,self.out + ".recode.vcf",self.out)
        return gene_fst

    def fst_NonCoding(self):
        for j in Non_code:
            if int(j[1]) > int(self.start) + 200000 and int(j[2]) - int(j[1]) > self.length:
                chr = j[0]
                start = j[1]
                end = j[2]
                out = self.out + "_NonCoding"
                fst = "vcftools --vcf M.vcf --chr {} --from-bp {} --to-bp {} --recode --recode-INFO-all --out {} && vcftools --vcf {} --fst-window-size 50 --fst-window-step 20 --weir-fst-pop M.aspalax --weir-fst-pop M.psilurus --out {}"
                NonCoding_fst = fst.format(chr, start, end, out, out + ".recode.vcf", out)
                break
        return NonCoding_fst


if __name__ == "__main__":
    """
    How to run:
    python Generate_Pairs_Fst.py gff
    """
    startss = time.time()
    for gene,value in Match_genes.items():
        gene = gene
        chr = value[0]
        chr = re.search(r"[0-9]+",chr)[0]
        start = int(value[1])
        end = int(value[2])
        length = end - start
        f = Fst(gene=gene,chr=chr,start=start,end=end,length=length,out=gene)
        print(f.fst_gene())
        print(f.fst_NonCoding())
    endss = time.time()
    print(f"the running time is {endss - startss}s")  ####字符串格式化新方法