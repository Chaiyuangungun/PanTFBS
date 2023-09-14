from scipy.stats import levene, ttest_ind
import argparse
import numpy as np
import math
def get_gene_id(gene_ids_file):#获得要查找的gene_id
    gene_ids = []
    with open(gene_ids_file,"r") as f:
        for line in f:
            gene_ids.append(line.strip())
    return gene_ids

def write_TFnumbers(samples_file,jaspar_database_pathway,gene_ids,out_file)  :
    TFnums = {}
    TFsums = {}
    TFnames = []
    samples = []
    family = {}
    fams = []
    with open(samples_file,"r") as f1:
        for line in f1:
            lines = line.strip().split()
            fam = lines[0]
            family[fam] = []
            fams.append(fam)
        fams    = list(set(fams))    
    with open(samples_file,"r") as f1:
        for line in f1:
            lines = line.strip().split()
            sample = lines[1]
            fam = lines[0]
            family[fam].append(sample)
            samples.append(sample)
    with open(jaspar_database_pathway,"r") as f2:
        for line in f2:
            lines = line.strip().split()
            if lines[0] == "ID":
                continue
            TFname = lines[1]
            TFnames.append(TFname)
        TFnames = list(set(TFnames))
    for sample in samples:
        TFnums[sample] = {}
        TFsums[sample] = {}
        for geneid in gene_ids:
            TFnums[sample][geneid] = {}
            TFsums[sample][geneid] = {}
            for TFname in TFnames:
                TFnums[sample][geneid][TFname] = 0
                TFsums[sample][geneid][TFname] = 0
    for sample in samples:
        with open(sample+".TFname.promoter.sites","r") as f3:
            for line in f3:
                lines = line.strip().split()
                if lines[0] == "geneid" :
                    continue
                geneid = lines[0]
                TFname = lines[1]
                TFnums[sample][geneid][TFname] = 1 
                TFsums[sample][geneid][TFname] += 1
    nums = {}
    Pvalue = {}    
    for geneid in gene_ids:
        nums[geneid] = {}
        Pvalue[geneid] = {}
        for TFname in TFnames:
            nums[geneid][TFname] = {}
            Pvalue[geneid][TFname] = {}
            for fam in fams:
                nums[geneid][TFname][fam] = []
                for sample in family[fam]:
                    nums[geneid][TFname][fam].append(TFnums[sample][geneid][TFname])
            ref_list = nums[geneid][TFname][fams[0]]
            alt_list = nums[geneid][TFname][fams[1]]
            levene_val = levene(ref_list, alt_list).pvalue
            if levene_val > 0.05:
                equal_var = True
            else:
                equal_var = False
            Pvalue[geneid][TFname] = ttest_ind(ref_list, alt_list, equal_var=equal_var).pvalue
    with open(out_file+".Pvalue","w") as f:
        f.write("geneid\tTFname\t")
        for sample in samples:
            f.write(sample+"\t")
        f.write("Pvalue\n")
        for geneid in gene_ids:
            for TFname in TFnames:
                if Pvalue[geneid][TFname] < 0.05 :
                    f.write(geneid+"\t")
                    f.write(TFname+"\t")
                    for sample in samples:
                        f.write(TFsums[sample][geneid][TFname]+"\t")
                    f.write(str(Pvalue[geneid][TFname])+"\t")
                  

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("-d","--database", type=str,default="/share/home/stu_chaikun/data/Script/python_Script/jaspar_database.csv")#jaspar数据库
parser.add_argument("-i","--id", type=str)#geneid文件
parser.add_argument("-s","--samples", type=str)#samples
parser.add_argument("-o","--out", type=str, default="out_file")#输出文件

args = parser.parse_args()
jaspar_database_pathway = args.database
gene_ids_file = args.id 
samples_file = args.samples
out_file = args.out

gene_ids = get_gene_id(gene_ids_file)
write_TFnumbers(samples_file,jaspar_database_pathway,gene_ids,out_file)
