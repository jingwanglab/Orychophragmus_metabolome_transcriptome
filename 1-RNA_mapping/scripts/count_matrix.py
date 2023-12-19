# -*- coding: utf-8 -*-
"""
Changfu Jia : 1020160171@qq.com

This is a scripts for merging htseq count file 
"""
import os
import re
import pybedtools 
import sys
from multiprocessing import Pool
pybedtools.set_bedtools_path('/usr_storage/software/bedtools2/bin/')


def multi_bedtools(gene, bed, leng):
    #print(gene+"sub process active")
    exon_len=0;
    bed_inter=pybedtools.bedtool.BedTool(bed,from_string=True).merge()
    for i in range(0,len(bed_inter)):
        exon_len+=(bed_inter[i]['end']-bed_inter[i]['start']+1)
    #return(exon_len)
    leng[gene]=exon_len
    #print(gene+"sub process finished")
   
    
   
def calculate_merged_exon_length(gtf):
#lenth = end - start + 1  for gtf or gff format 
    s={}
    for line in open(gtf):
        data=line.strip().replace(".TAIR10","").replace(".v1.1","").replace(".v2.1","").split()
        if data[2]=='exon':
            id_patt='\ttranscript_id\s+"([^"]*)";\s+gene_id\s"([^"]*)";'
            #id_patt='\ttranscript_(id)\s+'
            id=re.findall(id_patt, line.strip().replace(".TAIR10","").replace(".v1.1","").replace(".v2.1",""), re.S)
            #print(id[0][0])
            if id[0][1] in s:
                s[id[0][1]]+=data[0]+'\t'+data[3]+'\t'+data[4]+'\t'+id[0][1]+'\n'
            else:
                s[id[0][1]]=data[0]+'\t'+data[3]+'\t'+data[4]+'\t'+id[0][1]+'\n'
                    #s[l[0]]=l[0]
        else:
            continue

#test pybedtools
#c=pybedtools.bedtool.BedTool(s['Ovio00003'],from_string=True)
#cc=c.merge()
#cc[0]
#cc[4]
#cc[4]['end']-cc[4]['start']+1
#len(c)


#foreach 
    leng={}
    pool=Pool(128)
    for gene in s.keys():
        pool.apply_async(multi_bedtools(gene, s[gene], leng))
    pool.close()
    pool.join()
        #bed_inter=pybedtools.bedtool.BedTool(bed,from_string=True).merge()
        #for i in range(0,len(bed_inter)):
            #if gene in leng:
                #leng[gene]+=(bed_inter[i]['end']-bed_inter[i]['start']+1)
            #else:
                #leng[gene]=(bed_inter[i]['end']-bed_inter[i]['start']+1)
                
    return leng

   

    

def generate_count_matrix_for_fpkm(gtf, start, path, output):
    #get sum of merged exon length of each gene
    
    gtfDB = {}
    gtfDB = calculate_merged_exon_length(gtf)
    #gtfDB = {}
    #gtfDB = leng.copy()
    
    #stats each counts file and print header
    countsDB = {}
    g=''
    for file in os.listdir(path):
        if bool(re.search(start,file)) == True:
            for line in open(path+'/'+file):
                data=line.strip().replace(".TAIR10","").replace(".v1.1","").replace(".v2.1","").split()
                if bool(re.match("_",data[0])) == True:
                    continue
                else:
                    if data[0] in countsDB:
                    #countsDB[data[0]]={}
                        countsDB[data[0]][file]=data[1]
                    else:
                        countsDB[data[0]]={}
                        countsDB[data[0]][file]=data[1]
                g=data[0]
        else:
            continue

    with open(output, 'w') as fOut:
        fOut.write('geneid')
        for i in sorted(countsDB[g].keys()):
            fOut.write('\t'+i)
        fOut.write("\tgene_length\n")
        
 
        for gene in sorted(countsDB.keys()):
            fOut.write(gene)
            for file in sorted(countsDB[gene].keys()):
                fOut.write("\t"+countsDB[gene][file])
            fOut.write("\t"+str(gtfDB[gene])+"\n")
  
        

    
    
###test debug

#path='/usr_storage/jcf/2.zhugecai/00.diploid/19.sub_noe_functionalization/4.counts_matrix/'
#start='P'
#gtf=''
#generate_count_matrix_for_fpkm(gtf, start, path)





if __name__ == "__main__":
    generate_count_matrix_for_fpkm(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    

#def generate_count_matrix_for_fpkm(gtf, start, path, output):






