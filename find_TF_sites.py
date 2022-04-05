import re
import argparse
from functools import partial
from multiprocessing.pool import Pool
import time

#######################################################
########################编写函数########################
######################################################
def get_PFM(jaspar_database_pathway):#获得TF位置频率矩阵
    base_A = {}
    base_C = {}
    base_G = {}
    base_T = {}
    base_types = ["A","C","G","T"]
    gene_logo = {}
    IDs = []   
    with open(jaspar_database_pathway,"r") as f1:
        for line in f1:
            lines = line.strip().split("\t")
            if lines[0] == "ID":
                continue
            IDs.append(lines[0])
            base_A[lines[0]] = lines[-4]
            base_C[lines[0]] = lines[-3]
            base_G[lines[0]] = lines[-2]
            base_T[lines[0]] = lines[-1]
    for id in IDs:
        gene_logo[id] = {}
        bases = {}
        As = base_A[id].strip().split(";")[:-1]
        Cs = base_C[id].strip().split(";")[:-1]
        Gs = base_G[id].strip().split(";")[:-1]
        Ts = base_T[id].strip().split(";")[:-1]
        for num in range(len(As)):
            gene_logo[id][num] = ["",""]
            bases[num] = {}
            A = int(As[num])
            C = int(Cs[num])
            G = int(Gs[num])
            T = int(Ts[num])
            sum = A+C+G+T
            bases[num]["A"] = A
            bases[num]["C"] = C
            bases[num]["G"] = G
            bases[num]["T"] = T
            for base_type in base_types:              
                gene_logo[id][num][0] += base_type+"|"
                gene_logo[id][num][1] += str(format(float(bases[num][base_type]/sum),".3f"))+"|"
            gene_logo[id][num][0] = gene_logo[id][num][0].strip("|")
            gene_logo[id][num][1] = gene_logo[id][num][1].strip("|")
    return gene_logo

def get_threshold_80(PFMs):#将临界值设为80%，来确定每个位点的是否可以错配，及其可以错配的碱基类型
    TF_IDs = []
    threshold_80_result = {}
    for id in PFMs:
        TF_IDs.append(id)
    for id in TF_IDs:
        threshold_80_result[id] = []
        for sites in PFMs[id]:
            result = ""
            site_type = {}
            base_type = PFMs[id][sites][0].split("|")
            frequencys = PFMs[id][sites][1].split("|")
            for num in range(len(base_type)):
                site_type[base_type[num]] = frequencys[num]
            frequency_threshold = format(float(max(frequencys))*0.8,".3f")
            for base in ["A","C","G","T"]:
                if site_type[base] < frequency_threshold:
                    continue
                result += base+"|"
            threshold_80_result[id].append(result.strip("|"))    
    return threshold_80_result

def get_gene_promoter(genome_file,gff_file) :#输入基因组文件和gff文件，获得每个基因的启动子区（上游2500bp）序列字典
    promoter_fasta = {}    
    global gene_gff
    with open(genome_file,"r") as f1:
        genome = {}
        for line in f1:
            if ">" in line:
                chr_num = line.strip().strip(">")
                genome[chr_num] = ""
                continue
            genome[chr_num] += line.strip()
    with open(gff_file,"r") as f2:
        gene_gff = {}
        for line in f2:
            if "gene" in line :
                lines = line.strip().split()
                chr = lines[0]
                start = lines[3]
                end = lines[4]
                direction = lines[6]
                id = lines[-1].split(";")[0].replace("ID=","")
                gene_gff[id] = [chr,start,end,direction]
    for id in gene_gff:
        chr = gene_gff[id][0]
        start = gene_gff[id][1]
        end = gene_gff[id][2]
        direction = gene_gff[id][3]
        if direction == "+":
            promoter_fasta[id] = genome[chr][int(start)-2500:int(start)]
        else :
            promoter_fasta[id] = genome[chr][int(end):int(end)+2500]
            promoter_fasta[id] = promoter_fasta[id][::-1]
    return promoter_fasta

def get_TF_fasta(threshold_80_logo):#得到TF阈值为80%以上的具体序列
    threshold_80_fasta = {}
    for TF_id in threshold_80_logo:
        threshold_80_fasta[TF_id] = ""
        for base in threshold_80_logo[TF_id]:
            if "|" in base :
                threshold_80_fasta[TF_id] += "("+base+")"
            else:
                threshold_80_fasta[TF_id] += base
    return threshold_80_fasta

def get_gene_id(gene_ids_gile):#获得要查找的gene_id
    gene_ids = []
    with open(gene_ids_gile,"r") as f:
        for line in f:
            gene_ids.append(line.strip())
    return gene_ids

def get_promoter_motif_sites(threshold_80_fasta,gene_ids):#启动子区的具体TF结合位点
    TF_longs = {}
    for TF_id in threshold_80_fasta:
        motif_long = 0
        for TF_long in threshold_80_fasta[TF_id]:
            if TF_long == "(" or TF_long == ")":
                continue
            elif TF_long == "|":
                motif_long = motif_long -1
            else:
                motif_long += 1
        TF_longs[TF_id] = motif_long
    TF_site[gene_ids] = {}
    fasta = promoter_fastas[gene_ids]      
    promoter_long = len(fasta)
    for TF_id in threshold_80_fasta:
        TF_site[gene_ids][TF_id] = []
        for num in range(promoter_long-TF_longs[TF_id]+1):  
            if gene_gff[gene_ids][3] == "+":
                if re.search(threshold_80_fasta[TF_id],promoter_fastas[gene_ids][num:num+TF_longs[TF_id]]):
                    start  = int(gene_gff[gene_ids][1])-2500+num
                    end = int(gene_gff[gene_ids][1])-2500+num+TF_longs[TF_id]
                    TF_site[gene_ids][TF_id].append(str(start+1)+"\t"+str(end)+"\t"+promoter_fastas[gene_ids][num:num+TF_longs[TF_id]]+"\t"+"+")
                promoter_fasta_opposite = promoter_fastas[gene_ids][::-1]             
                if re.search(threshold_80_fasta[TF_id],promoter_fasta_opposite[num:num+TF_longs[TF_id]]):
                    start  = int(gene_gff[gene_ids][2])-num
                    end = int(gene_gff[gene_ids][2])-num-TF_longs[TF_id]
                    TF_site[gene_ids][TF_id].append(str(start+1)+"\t"+str(end)+"\t"+promoter_fasta_opposite[num:num+TF_longs[TF_id]]+"\t"+"-")
            if gene_gff[gene_ids][3] == "-":
                if re.search(threshold_80_fasta[TF_id],promoter_fastas[gene_ids][num:num+TF_longs[TF_id]]):
                    start  = int(gene_gff[gene_ids][2])+num
                    end = int(gene_gff[gene_ids][2])+num+TF_longs[TF_id]
                    TF_site[gene_ids][TF_id].append(str(start+1)+"\t"+str(end)+"\t"+promoter_fastas[gene_ids][num:num+TF_longs[TF_id]]+"\t"+"+")
                promoter_fasta_opposite = promoter_fastas[gene_ids][::-1]             
                if re.search(threshold_80_fasta[TF_id],promoter_fasta_opposite[num:num+TF_longs[TF_id]]):
                    end  = int(gene_gff[gene_ids][2])+2500-num
                    start = int(gene_gff[gene_ids][2])+2500-num-TF_longs[TF_id]
                    TF_site[gene_ids][TF_id].append(str(start+1)+"\t"+str(end)+"\t"+promoter_fasta_opposite[num:num+TF_longs[TF_id]]+"\t"+"-")
        if not TF_site[gene_ids][TF_id] :
            del TF_site[gene_ids][TF_id]
        
    return TF_site

def write_TF_sites(promoter_TF_sites):#书写结果文件  
    new_promoter_TF_sites = {}
    for dict in promoter_TF_sites:
        for gene_id in dict:
            new_promoter_TF_sites[gene_id] = {}
            for TF_id in dict[gene_id]: 
                new_promoter_TF_sites[gene_id][TF_id] = dict[gene_id][TF_id]
            new_promoter_TF_sites[gene_id][TF_id] = list(set(new_promoter_TF_sites[gene_id][TF_id]))
    with open(out_file,"w") as f1:
        for  gene_id in new_promoter_TF_sites:
            f1.write(gene_id+"\t")
            for TF_id in new_promoter_TF_sites[gene_id]:
                f1.write(TF_id+"\t")
            f1.write("\n")
            for TF_id in new_promoter_TF_sites[gene_id]:   
                for site in new_promoter_TF_sites[gene_id][TF_id]:
                    f1.write(TF_id+" : "+site+"\n")

#######################################################
########################设置参数########################
######################################################
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("-database", type=str,default="/share/home/stu_chaikun/data/Script/python_Script/jaspar_database.csv")#jaspar数据库
parser.add_argument("-genome", type=str)#基因组文件
parser.add_argument("-gff", type=str)#gff文件
parser.add_argument("-id", type=str)#geneid文件
parser.add_argument("-thread", type=int, default=10)#并行进程数
parser.add_argument("-out", type=str, default="out_file")#输出文件
args = parser.parse_args()

#######################################################
########################调用函数########################
######################################################
jaspar_database_pathway = args.database
genome_file = args.genome 
gff_file = args.gff 
gene_ids_file = args.id 
thread_num = args.thread
out_file = args.out
start_time = time.time() 
PFMs = get_PFM(jaspar_database_pathway)#jaspar数据库位置
threshold_80_logo = get_threshold_80(PFMs)
promoter_fasta = get_gene_promoter(genome_file,gff_file)
global promoter_fastas 
global sum_num
global TF_site
gene_ids = get_gene_id(gene_ids_file)
TF_site = {}
promoter_fastas = promoter_fasta 
#基因组序列文件和基因组gff文件
threshold_80_fasta = get_TF_fasta(threshold_80_logo)
pool = Pool(processes=thread_num)
pfunc = partial(get_promoter_motif_sites,threshold_80_fasta)
promoter_TF_sites = pool.map(pfunc,gene_ids)
write_TF_sites(promoter_TF_sites)
end_time = time.time()
print(end_time-start_time)
