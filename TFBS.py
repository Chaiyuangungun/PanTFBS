import argparse
from functools import partial
from multiprocessing.pool import Pool
import time
import math
from scipy.stats import levene, ttest_ind

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
                bases[num][base_type] = str(bases[num][base_type])        
                gene_logo[id][num][0] += base_type+"|"
                if bases[num][base_type] == "0":
                    gene_logo[id][num][1] += "-4.6|"                
                if  (bases[num][base_type] != "0" and base_type == "A") or (bases[num][base_type] != "0" and base_type == "T") :
                    gene_logo[id][num][1] += str(format(math.log(float(bases[num][base_type])/sum*5,math.e),".3f"))+"|"
                if (bases[num][base_type] != "0" and base_type == "G") or (bases[num][base_type] != "0" and base_type == "C"):
                    gene_logo[id][num][1] += str(format(math.log(float(bases[num][base_type])/sum/3*10,math.e),".3f"))+"|"
            gene_logo[id][num][0] = gene_logo[id][num][0].strip("|")
            gene_logo[id][num][1] = gene_logo[id][num][1].strip("|")
    return gene_logo

def get_threshold_80(PFMs):#
    TF_IDs = []
    threshold_80_result = {}
    for id in PFMs:
        TF_IDs.append(id)
    for id in TF_IDs:
        sum = 0
        for sites in PFMs[id]:
            frequencys = PFMs[id][sites][1].split("|")
            max_frequencys = max(frequencys)
            sum += float(max_frequencys)     
        threshold_80_result[id] = format(sum,".3f") 
    return threshold_80_result

def get_gene_promoter(genome_file,gff_file) :#输入基因组文件和gff文件，获得每个基因的启动子区（上游2000bp）序列字典
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
                if "path" in id :
                    id = id.split(".")[0]
                gene_gff[id] = [chr,start,end,direction]
    for id in gene_gff:
        chr = gene_gff[id][0]
        start = gene_gff[id][1]
        end = gene_gff[id][2]
        direction = gene_gff[id][3]
        if direction == "+":
            promoter_fasta[id] = genome[chr][int(start)-2000:int(start)]
        if direction == "-":
            promoter_fasta[id] = genome[chr][int(end):int(end)+2000]
            promoter_fasta[id] = promoter_fasta[id][::-1]
    return promoter_fasta

def get_gene_id(gene_ids_file):#获得要查找的gene_id
    gene_ids = []
    with open(gene_ids_file,"r") as f:
        for line in f:
            gene_ids.append(line.strip())
    return gene_ids

def get_promoter_motif_sites(PFMs,gene_id):#启动子区的具体TF结合位点
    promoter_motif_sites = {}
    promoter_motif_sites[gene_id] = {}
    for TF_id in PFMs:  
        promoter_motif_sites[gene_id][TF_id] = []
        longs = len(PFMs[TF_id])
        try:
            fasta = promoter_fastas[gene_id]#正向
            F_fasta = promoter_fastas[gene_id][::-1]#反向
        except:
            continue
        for num in range(len(fasta)-longs):#正向
            sum = 0
            R_sum = 0
            search_fasta = fasta[num:num+longs]
            R_search_fasta = search_fasta.replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper()
            for base in range(len(search_fasta)):
                A_type = PFMs[TF_id][base][1].split("|")[0]
                C_type = PFMs[TF_id][base][1].split("|")[1]
                G_type = PFMs[TF_id][base][1].split("|")[2]
                T_type = PFMs[TF_id][base][1].split("|")[3]
                if search_fasta[base] == "A":
                    sum += float(A_type)
                    R_sum += float(T_type)#互补
                if search_fasta[base] == "G":
                    sum += float(G_type)
                    R_sum += float(C_type)
                if search_fasta[base] == "C":
                    sum += float(C_type)
                    R_sum += float(G_type)
                if search_fasta[base] == "T":
                    sum += float(T_type)
                    R_sum += float(A_type)
            if sum >= float(threshold_80_logo[TF_id])*threshold:
                score = format(sum/float(threshold_80_logo[TF_id]),".3f")
                if gene_gff[gene_id][3] == "+":#
                    start = int(gene_gff[gene_id][1])
                    chr = gene_gff[gene_id][0]
                    promoter_motif_sites[gene_id][TF_id].append(chr+"\t"+str(start-2000+num+1)+"\t"+str(start-2000+longs+num+1)+"\t+\t"+search_fasta+"\t"+str(score))

            if R_sum  >= float(threshold_80_logo[TF_id])*threshold:
                R_score = format(R_sum/float(threshold_80_logo[TF_id]),".3f")
                
                if gene_gff[gene_id][3] == "-":#
                    end = int(gene_gff[gene_id][2])
                    chr = gene_gff[gene_id][0]
                    promoter_motif_sites[gene_id][TF_id].append(chr+"\t"+str(end+num+1)+"\t"+str(end+longs+num+1)+"\t+\t"+R_search_fasta+"\t"+str(R_score))
        for num in range(len(fasta)-longs):#反向
            sum = 0
            R_sum = 0
            F_search_fasta = F_fasta[num:num+longs]
            RF_search_fasta = F_search_fasta.replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper()
            for base in range(len(F_search_fasta)):
                A_type = PFMs[TF_id][base][1].split("|")[0]
                C_type = PFMs[TF_id][base][1].split("|")[1]
                G_type = PFMs[TF_id][base][1].split("|")[2]
                T_type = PFMs[TF_id][base][1].split("|")[3]
                if F_search_fasta[base] == "A":
                    sum += float(A_type)
                    R_sum += float(T_type)#互补
                if F_search_fasta[base] == "G":
                    sum += float(G_type)
                    R_sum += float(C_type)#互补
                if F_search_fasta[base] == "C":
                    sum += float(C_type)
                    R_sum += float(G_type)#互补
                if F_search_fasta[base] == "T":
                    sum += float(T_type)
                    R_sum += float(A_type)#互补
            if sum >= float(threshold_80_logo[TF_id])*threshold:
                score = format(sum/float(threshold_80_logo[TF_id]),".3f")
                if gene_gff[gene_id][3] == "-":#
                    end = int(gene_gff[gene_id][2])
                    chr = gene_gff[gene_id][0]
                    promoter_motif_sites[gene_id][TF_id].append(chr+"\t"+str(end+2000-num-1)+"\t"+str(end+2000-longs-num-1)+"\t-\t"+F_search_fasta[::-1].replace("A","t").replace("T","a").replace("G","c").replace("C","g").upper()+"\t"+str(score))
            if R_sum >= float(threshold_80_logo[TF_id])*threshold:
                R_score = format(R_sum/float(threshold_80_logo[TF_id]),".3f")
                if gene_gff[gene_id][3] == "+":#
                    start = int(gene_gff[gene_id][1])
                    chr = gene_gff[gene_id][0]
                    promoter_motif_sites[gene_id][TF_id].append(chr+"\t"+str(start-num-1-longs)+"\t"+str(start-num-1)+"\t-\t"+F_search_fasta[::-1]+"\t"+str(R_score))
   
        if not promoter_motif_sites[gene_id][TF_id] :
            del promoter_motif_sites[gene_id][TF_id]
    return promoter_motif_sites

def write_TF_sites(promoter_TF_sites,out_file):#书写结果文件  
    new_promoter_TF_sites = {}
    TF_ids = []
    for dict in promoter_TF_sites:
        for gene_id in dict:
            new_promoter_TF_sites[gene_id] = {}
            for TF_id in dict[gene_id]: 
                new_promoter_TF_sites[gene_id][TF_id] = dict[gene_id][TF_id]
                TF_ids.append(TF_id)
                new_promoter_TF_sites[gene_id][TF_id] = list(set(new_promoter_TF_sites[gene_id][TF_id]))
    TF_ids = list(set(TF_ids))
    with open(out_file+".TFids.genome.sites","w") as f1:
        f1.write("geneid\tTFid\tchr\tstart\tend\tdirection\tTFseq\tscore\n")
        for  gene_id in new_promoter_TF_sites:
            for TF_id in new_promoter_TF_sites[gene_id]:   
                for site in new_promoter_TF_sites[gene_id][TF_id]:
                    f1.write(gene_id+"\t"+TF_id+"\t"+site+"\n")    
    return new_promoter_TF_sites,TF_ids

####鉴定完转录因子结合位点肯定要找出差异，比如我有一组经过胁迫处理的品种，有些差异表达基因，我需要搞清楚某个转录因子调控那个基因，而且我又可以在差异表达
####的基因中找到这个转录因子，那么就可以绘制一个调控网络出来，比如某个转录因子可能在这个胁迫条件下调控某些基因以应对这种胁迫环境
def get_TF_Gene(new_promoter_TF_sites,TF_ids):##new_promoter_TF_sites[gene_id][TF_id]
    TF_Gene_ids = {}
    for id in TF_ids:
        TF_Gene_ids[id] = []
        for gene_id in new_promoter_TF_sites:
            for TF_id in new_promoter_TF_sites[gene_id] :  
                if TF_id == id:
                    TF_Gene_ids[id].append(gene_id)
    return TF_Gene_ids

####得到一个文件，以TFgenes为后缀，每行第一列记录转录因子id，第二列记录该转录因子在提供的基因id列表中，多少个基因有该转录因子的结合位点；后边列数记录geneid

#####TF_sites更新为转录因子名字对应该转录因子在启动子区的起始终止位置
def draw(new_promoter_TF_sites,gene_gff,jaspar_database_pathway,out_file) :
    gene_ids = []
    TF_ids = []
    TF_id_name = {}
    TF_sites = {}
    with open(jaspar_database_pathway,"r") as f1:
        for line in f1:
            lines = line.strip().split()
            if lines[0] == "ID":
                continue
            TF_id = lines[0]
            TF_name = lines[1]
            TF_id_name[TF_id] = TF_name    
    for gene_id in new_promoter_TF_sites:
        gene_ids.append(gene_id)
        for TF_id in new_promoter_TF_sites[gene_id]:
            TF_ids.append(TF_id)
    for gene_id in gene_ids:
        try:
            start = gene_gff[gene_id][1]
            end = gene_gff[gene_id][2]
            direction = gene_gff[gene_id][3]
        except:
            continue
        TF_sites[gene_id] = {}
        for TF_id in new_promoter_TF_sites[gene_id]:
            TF_sites[gene_id][TF_id_name[TF_id]] = []
            for site in new_promoter_TF_sites[gene_id][TF_id]:
                sites = site.strip().split()
                if direction == "+":
                    new_site_start = str(int(sites[1]) - int(start))
                    new_site_end = str(int(sites[2]) - int(start))
                    new_site_direction = sites[3]
                    new_site_seq = sites[4]
                    new_site_score = sites[5]
                if direction == "-":
                    new_site_start = str(int(end) - int(sites[2]))
                    new_site_end = str(int(end) - int(sites[1]))
                    new_site_direction = sites[3]
                    new_site_seq = sites[4]
                    new_site_score = sites[5]
                new_site = new_site_start+"\t"+new_site_end+"\t"+new_site_direction+"\t"+new_site_seq+"\t"+new_site_score
                TF_sites[gene_id][TF_id_name[TF_id]].append(new_site)  ###TF_sites改TF_id_name为family
    with open(out_file+".TFname.promoter.sites","w") as f1:
        f1.write("geneid\tTFname\tstart\tend\tdirection\tTFseq\tscore\n")
        for geneid in TF_sites:
            for TFname in TF_sites[geneid]:
                for newsite in TF_sites[geneid][TFname]:
                    f1.write(geneid+"\t"+TFname+"\t"+"\t"+newsite+"\n")

#######################################################
########################设置参数########################
######################################################
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument("-d","--database", type=str,default="/share/home/stu_chaikun/data/Script/python_Script/jaspar_database.csv")#jaspar数据库
parser.add_argument("-g","--genome", type=str)#基因组文件
parser.add_argument("-a","--gff", type=str)#gff文件
parser.add_argument("-i","--id", type=str)#geneid文件
parser.add_argument("-k","--threshold", type=float, default=0.8)#阈值
parser.add_argument("-t","--thread", type=int, default=10)#并行进程数
parser.add_argument("-o","--out", type=str, default="out_file")#输出文件
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
global threshold
threshold = args.threshold
start_time = time.time() 
PFMs = get_PFM(jaspar_database_pathway)#jaspar数据库位置
global  threshold_80_logo
threshold_80_logo = get_threshold_80(PFMs)
promoter_fasta = get_gene_promoter(genome_file,gff_file)
global promoter_fastas 
global TF_site
gene_ids = get_gene_id(gene_ids_file)
TF_site = {}
promoter_fastas = promoter_fasta 
#基因组序列文件和基因组gff文件
pool = Pool(processes=thread_num)
pfunc = partial(get_promoter_motif_sites,PFMs)
promoter_TF_sites = pool.map(pfunc,gene_ids)
new_promoter_TF_sites,TF_ids = write_TF_sites(promoter_TF_sites,out_file)
TF_Gene_ids = get_TF_Gene(new_promoter_TF_sites,TF_ids)
draw(new_promoter_TF_sites,gene_gff,jaspar_database_pathway,out_file)
end_time = time.time()
print(end_time-start_time)
