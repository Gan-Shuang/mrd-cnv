import sys
import os
import argparse
from subprocess import call
from subprocess import check_output
from subprocess import check_call
import pymysql
import yaml
import pandas as pd
import numpy as np

def load_yml(dir):
    file=open(dir+"/parameter.yml","r")
    file_data=file.read()
    file.close()
    config=yaml.load(file_data,Loader=yaml.FullLoader)
    return(config)

def cnv_info(user,passwd):
    config = {
          'host':'',
          'port':,
          'user':user,
          'password':passwd,
          'database':'',
          'charset':'utf8mb4',
          'cursorclass':pymysql.cursors.Cursor
          }
    db = pymysql.connect(**config)
    cursor = db.cursor()
    cursor.execute("SELECT * FROM mutation_region")
    cnv_info = cursor.fetchall()
    gene_list=[]
    for mutations in cnv_info:
        gene=mutations[1]
        gene_group=mutations[6]
        if gene_group=='CNV':
            gene_list.append(gene)
    return(gene_list)

def run_readCounter(inputbam,sample_name,outputdir):
    run_readCounter_commond=("readCounter --window 1000000 --quality 20 --chromosome '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y' "+
                         inputbam+" > "+outputdir+"/"+sample_name+".wig")
    check_call(run_readCounter_commond,shell=True)

def cor_wig(outputdir,sample_name,correct_ratio):
    ratio=correct_ratio
    file=outputdir+"/"+sample_name+".wig"
    out_file=outputdir+"/"+sample_name+".cor.wig"
    wig_read=open(file,"r")
    wig_write=open(out_file,"w")
    wig_chr_read=pd.DataFrame()
    for i in wig_read.readlines():
        i=i.rstrip()
        try:
            chr=i.split(" ")[1].split("=")[1]
        except IndexError:
            pass
        try:
            wig_chr_read=wig_chr_read.append(pd.DataFrame({"chrom":[chr],"depth":[int(i)]}))
        except ValueError:
            wig_chr_read=wig_chr_read.append(pd.DataFrame({"chrom":[chr],"depth":[i]}))
    wig_chr_read=wig_chr_read.reset_index(drop=True)
    wig_chr_read.index=wig_chr_read.index+1
    out_wig=pd.DataFrame()
    out_put=[]
    for i in wig_chr_read.index:
        if wig_chr_read.loc[i].values[0]=='19' and isinstance(wig_chr_read.loc[i].values[1],int)==True:
            num_read=round(int(wig_chr_read.loc[i].values[1])*ratio,3)
            out_wig=out_wig.append(pd.DataFrame({"wig":[num_read]}))
            out_put.append(num_read)
        else:
            out_wig=out_wig.append(pd.DataFrame({"wig":[wig_chr_read.loc[i].values[1]]}))
            out_put.append(wig_chr_read.loc[i].values[1])
    for i in out_put:
        wig_write.write(str(i)+"\n")
    wig_write.close()

def run_inchorCNA(sample_name,outputdir,work_dir,sensitivity,maxCN):
    run_ichorcna_commod=('Rscript '+str(work_dir)+'/scripts/runIchorCNA.R --id '+str(sample_name)+' --WIG '+str(outputdir)+"/"+str(sample_name)+
                         '.cor.wig --ploidy "c(2,3)" --normal "'+str(sensitivity)+'" --maxCN '+str(maxCN)+
                         ' --gcWig '+str(work_dir)+'/ref/gc_hg19_1000kb.wig '+
                         ' --mapWig '+str(work_dir)+'/ref/map_hg19_1000kb.wig '+
                         ' --centromere '+str(work_dir)+'/ref/GRCh37.p13_centromere_UCSC-gapTable.txt '+
                         ' --normalPanel '+str(work_dir)+'/ref/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds '+
                         ' --includeHOMD False --chrs "c(1:22)" --chrTrain "c(1:22)" '+
                         ' --estimateNormal True --estimatePloidy True --estimateScPrevalence FALSE '+
                         ' --scStates "c()" --txnE 0.9999 --txnStrength 10000 --plotFileType=png --outDir '+str(outputdir))
    print(run_ichorcna_commod)
    check_call(run_ichorcna_commod,shell=True)

def tumor_fraction_ploidy(sample_name,outputdir):
    tumor_fraction_commond="sed -n '2p' "+outputdir+"/"+sample_name+".params.txt|awk -F'\t' '{print $2}'"
    tumor_ploidy_commond="sed -n '2p' "+outputdir+"/"+sample_name+".params.txt|awk -F'\t' '{print $3}'"
    tumor_fraction=check_output(tumor_fraction_commond,shell=True).decode().rstrip()
    tumor_ploidy=check_output(tumor_ploidy_commond,shell=True).decode().rstrip()
    return(tumor_fraction,tumor_ploidy)

def get_gene_info(outputdir,sample_name):
    get_bed_commond="awk -F'\t' -v OFS='\t' '{if($11!=2)print $2,$3,$4}' "+outputdir+"/"+sample_name+".seg.txt"+" > "+outputdir+"/"+sample_name+"_cna.bed"
    call(get_bed_commond,shell=True)
    rm_header_commond="sed -i '1d' "+outputdir+"/"+sample_name+"_cna.bed"
    call(rm_header_commond,shell=True)
    anno_bed_commond=("bedtools intersect -a ref/annotation_genesymbol.bed -b "+outputdir+"/"+
                      sample_name+"_cna.bed"+"|awk -F '\t' '{print $4}' |sort|uniq > "+
                      outputdir+"/"+sample_name+"_gene_list.txt")
    call(anno_bed_commond,shell=True)
    get_common_genes_commond=("cat "+outputdir+"/"+sample_name+"_gene_list.txt "+outputdir+
                              "/anno_gene.txt|sort|uniq -d > "+outputdir+"/reported_genes.txt")
    call(get_common_genes_commond,shell=True)
    get_geneinfo_commond=("for i in `cat "+outputdir+
                          "/reported_genes.txt`;do grep -P '\t'${i}'$' ref/annotation_genesymbol.bed;done > "+
                          outputdir+"/gene_info.txt")
    call(get_geneinfo_commond,shell=True)

def cnv_gene(outputdir,sample_name):
    gene_list=[]
    for i in open(outputdir+"/tmp/gene_info.txt","r").readlines():
        gene=i.split("\t")[3].rstrip()
        chr1=i.split("\t")[0].rstrip()
        start=i.split("\t")[1].rstrip()
        end=i.split("\t")[2].rstrip()
        for line in open(outputdir+"/tmp/"+sample_name+".seg.txt","r").readlines():
            chr2=line.split("\t")[1].rstrip()
            start1=line.split("\t")[2].rstrip()
            end1=line.split("\t")[3].rstrip()
            num=line.split("\t")[10].rstrip()
            event=line.split("\t")[11].rstrip()
            if chr1==chr2 and int(start)>int(start1) and int(end)<int(end1) :
#                 print(gene+"\t"+start1+"\t"+end1+"\t"+event+"\t"+num+"\t"+str(int(end)-int(start)))
                gene_list.append(gene+"\t"+event+"\t"+num+"\t"+str(int(end)-int(start)))
            elif chr1==chr2 and int(start)<int(start1) and int(end)>int(start1) :
#                 print(gene+"\t"+start1+"\t"+end1+"\t"+event+"\t"+num+"\t"+str(int(end)-int(start1)))
                gene_list.append(gene+"\t"+event+"\t"+num+"\t"+str(int(end)-int(start1)))
            elif chr1==chr2 and int(start)<int(end1) and int(end)>int(end1) :
#                 print(gene+"\t"+start1+"\t"+end1+"\t"+event+"\t"+num+"\t"+str(int(end1)-int(start)))
                gene_list.append(gene+"\t"+event+"\t"+num+"\t"+str(int(end1)-int(start)))
    gene_list=sorted(list(set(sorted(gene_list))))
    #############################################################################
    last_gene=""
    last_len=""
    set_gene_list=[]
    for i in gene_list:
        gene=i.split("\t")[0]
        num=i.split("\t")[2]
        length=i.split("\t")[3]
        if gene!=last_gene:
            set_gene_list.append(i)
        if gene==last_gene:
            if length>last_len:
                set_gene_list.pop()
                set_gene_list.append(i)
            if length<last_gene:
                pass
        last_gene=gene
        last_len=length
    ###############################################################################
    file=open(outputdir+"/"+sample_name+".gene_fraction","w")
    file.write("ApprovedSymbol"+"\t"+"CopyNumberVariant"+"\t"+"CopyNumber"+"\n")
    for i in set_gene_list:
        gene=i.split("\t")[0]
        num=i.split("\t")[2]
        if str(num)!='NA':
            if int(num)<2:
                file.write(gene+"\tloss\t"+num+"\n")
            if int(num)>2:
                file.write(gene+"\tgain\t"+num+"\n")
    file.close()

def chr_region(outputdir,sample_name):
    file=outputdir+"/tmp/"+sample_name+".seg.txt"
    result=open(outputdir+"/"+sample_name+".chr_fraction","w")
    result.write("Chr\tStart\tEnd\tCopyNumber\tAnnotaionGene\n")
    for line in open(file,"r").readlines():
        chrom=line.rstrip().split("\t")[1]
        start=line.rstrip().split("\t")[2]
        end=line.rstrip().split("\t")[3]
        copy_event=line.rstrip().split("\t")[11]
        copy_number=line.rstrip().split("\t")[10]
        if copy_number!='Corrected_Copy_Number' and int(copy_number)!=2:
            info=chrom+"\t"+start+"\t"+end
            f1=open(outputdir+"/tmp/info.bed","w")
            f1.write(info+"\n")
            f1.close()
            commond=("bedtools intersect -b "+outputdir+"/tmp/info.bed -a ref/annotation_genesymbol.bed > "+
                    outputdir+"/tmp/gene.bed")
            call(commond,shell=True)
            gene_list=[]
            for info in open(outputdir+"/tmp/gene.bed","r").readlines():
                gene=info.rstrip().split("\t")[3]
                gene_list.append(gene)
            chr_gene=str(gene_list).replace("[","").replace("]","").replace(" '","").replace("'","")
            result.write(chrom+"\t"+start+"\t"+end+"\t"+copy_number+"\t"+chr_gene+"\n")
    result.close()

def plot_chr(seg_file,gene_file,ploidy,outputdir):
    if len(open(gene_file,"r").readlines())==0:
        plot_commond=("Rscript "+"scripts/plot_nogene.r "+
                  " --input "+seg_file+" --genes "+gene_file+" --ploidy "+ploidy+" --outdir "+outputdir)
        check_call(plot_commond,shell=True)
    if len(open(gene_file,"r").readlines())!=0:
        plot_commond=("Rscript "+"scripts/plot.r "+
                  " --input "+seg_file+" --genes "+gene_file+" --ploidy "+ploidy+" --outdir "+outputdir)
        check_call(plot_commond,shell=True)

parser=argparse.ArgumentParser(description="Analyze CNA wit ichorCNA")
parser.add_argument("-i","--inputbam",dest="bam_file",help="input bam file",required=True)
parser.add_argument("-o","-outputdir",dest="output_dir",help="output dir",required=True)
args = parser.parse_args()
sys.stdout.write("Bam file: "+args.bam_file+"\n")
sys.stdout.write("Output dir: "+args.output_dir+"\n")
input_bam=args.bam_file
output_dir=args.output_dir

if __name__=="__main__":
    image_run="singularity exec -B /mnt:/mnt,/fastzone:/fastzone ichorcna.sif"
    sample_name=input_bam.split("/")[-1].split(".")[0].split("-")[0]
    ##mkdir for output
    call(["mkdir",output_dir])
    call(["mkdir",output_dir+"/tmp"])
    tmpdir=output_dir+"/tmp"
    ##load_config
    work_dir=sys.path[0]
    config=load_yml(work_dir)
    mysql_user=config['mysql']['user']
    mysql_passwd=config['mysql']['passwd']
    correct_ratio=config['CHR19_correct_ratio'][0]
    sensitivity=config['ichorcna']['sensitivity']
    maxCN=config['ichorcna']['maxCN']
    #get anno genes
    anno_genes=cnv_info(mysql_user,mysql_passwd)
    get_anno_gene=open(tmpdir+"/anno_gene.txt","w")
    for gene in anno_genes:
        get_anno_gene.write(gene+"\n")
    get_anno_gene.close()
    ##run inchorCNA
    run_readCounter(input_bam,sample_name,tmpdir)
    cor_wig(tmpdir,sample_name,correct_ratio)
    attempts1 = 0
    success = False
    while attempts1 < 10 and not success:
        try:
            run_inchorCNA(sample_name,tmpdir,work_dir,sensitivity,maxCN)
            success = True
        except:
            attempts1 += 1
            if attempts1 == 10:
                break
    call(["cp",tmpdir+"/"+sample_name+"/"+sample_name+"_genomeWide.png",output_dir])
    ##get tumor fraction and ploidy
    tumor_fraction,tumor_ploidy=tumor_fraction_ploidy(sample_name,tmpdir)
    fraction_cnv=open(output_dir+"/"+sample_name+".fraction_cnv","w")
    fraction_cnv.write("Tumor_fraction_cnv\n"+tumor_fraction+"\n")
    fraction_cnv.close()
    ##get gene info
    get_gene_info(tmpdir,sample_name)
    ##annotation cnv genes
    cnv_gene(output_dir,sample_name)
    ##chr_region
    chr_region(output_dir,sample_name)
    ##plot chr-----##try 10 times
    attempts2 = 0
    success = False
    while attempts2 < 10 and not success:
        try:
            plot_chr(tmpdir+"/"+sample_name+".cna.seg",tmpdir+"/reported_genes.txt",tumor_ploidy,output_dir)
            success = True
        except:
            attempts2 += 1
            if attempts2 == 10:
                break
