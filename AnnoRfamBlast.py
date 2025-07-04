#!/usr/bin/env python
# -*-coding:utf-8-*-
import argparse
import sys
import time
import os
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import random
import math

def parse_args():
    parser = argparse.ArgumentParser(description="""用途 ： 
    为miRNA过滤后的reads与Rfam比对结果添加RNA类型""", formatter_class=argparse.RawDescriptionHelpFormatter, epilog="""使用举例： 
    python AnnoRfamBlast.py -infile /work2/users/chengjiangnan/Test/Whole_trans_test/SmallRNA_Result/1_CleanData/SRR19278766/SRR19278766.rfam.blastout \\
    -info /work1/DataBase/cjndatabase/Rfam/RF_seq.relation.txt \\
    -outfile /work2/users/chengjiangnan/Test/Whole_trans_test/SmallRNA_Result/1_CleanData/SRR19278766/SRR19278766.blast.result.txt
    """)
    parser.add_argument('-infile', help='输入样本与Rfam库的比对结果文件', dest='infile', required=True)
    parser.add_argument('-info', help='输入Rfam库中各个小RNA家族对应分类文件', dest='info', default='/work1/DataBase/cjndatabase/Rfam/RF_seq.relation.txt')
    parser.add_argument('-outfile', help='输出样本汇总结果文件', dest='outfile', required=True)
    argv = parser.parse_args()
    return argv

def read_info_table(info):
    dict_rfam = {}
    with open(info,'r') as f:
        for line in f:
            line = line.strip().split('\t')
            value_rfam = '#'.join(line[1:])
            dict_rfam[line[0]] = value_rfam
    return dict_rfam

def even_solve(tmp_sizes,tmp_labels,half_n):
    tmp_sizes.insert(0, 1)
    tmp_labels.insert(0, 'one')
    sizes = []
    labels = []
    for i in range(1,half_n + 1):
        sizes.append(tmp_sizes[i])
        sizes.append(tmp_sizes[-i])
        labels.append(tmp_labels[i])
        labels.append(tmp_labels[-i])
    return sizes,labels

def odd_solve(tmp_sizes,tmp_labels,half_n):
    tmp_sizes.insert(0, 1)
    tmp_labels.insert(0, 'one')
    sizes = []
    labels = []
    for i in range(1,half_n + 1):
        if i != half_n:
            sizes.append(tmp_sizes[i])
            sizes.append(tmp_sizes[-i])
            labels.append(tmp_labels[i])
            labels.append(tmp_labels[-i])
        else:
            sizes.append(tmp_sizes[i])
            labels.append(tmp_labels[i])
    return sizes,labels

def sort_list(old_sizes,old_labels):
    indexed_sizes = [(index, value) for index, value in enumerate(old_sizes)]
    indexed_sizes.sort(key=lambda x: x[1], reverse=True)
    sorted_index = [index for index, value in indexed_sizes]
    tmp_sizes = [value for index, value in indexed_sizes]
    tmp_labels = [old_labels[i] for i in sorted_index]
    if len(old_sizes) % 2 == 0:
        half_n = int(len(old_sizes)/2)
        sizes,labels = even_solve(tmp_sizes,tmp_labels,half_n)
    else:
        half_n = int(math.ceil(len(old_sizes)/2))
        sizes, labels = odd_solve(tmp_sizes,tmp_labels,half_n)
    return sizes, labels

def mkout(infile,info,outfile):
    indir,file = os.path.split(infile)
    prefix = file.split('.')[0]
    header = ['Query_Identifier','Subject_Identifier','Percentage_of_identity','Alignment_length','Number_of_mismatches',
              'gaps','Query_sequence_start','Query_sequence_end','Subject_sequence_start','Subject_sequence_end',
              'Expect_value','bit_score','Rfam_id','small_RNA_type']
    dict_rfam = read_info_table(info)
    with open(infile,'r') as fn,open(outfile,'w') as f:
        f.write('\t'.join(header) + '\n')
        for line in fn:
            line = line.strip().split('\t')
            rna_list = dict_rfam[line[1]].split('#')
            f.write('\t'.join(line + rna_list) + '\n')
    miRNA_file = os.path.join(indir,prefix + '.miRNA.id')
    n = 0
    list_total = []
    with open(outfile,'r') as fn,open(miRNA_file,'w') as f:
        for line in fn:
            if line.startswith('Query'):
                continue
            line = line.strip().split('\t')
            if 'Gene' not in line[13] and 'RNA' not in line[13]:
                continue
            ncRNA_list = [i.strip() for i in line[13].split(';') if i and i!='Gene' and 'RNA' in i]
            list_total.extend(ncRNA_list)
            if not len(ncRNA_list):
                continue
            n += 1
            if 'miRNA' in ncRNA_list:
                f.write(line[0] + '\n')
    old_labels = []
    old_sizes = []
    with open(f'{indir}/{prefix}.RNA.type.txt','w') as f:
        f.write('\t'.join(['RNA_type',f'{prefix} _total',f'{prefix}_total(%)']) + '\n')
        for name in Counter(list_total):
            rna_percent = str(round(int(Counter(list_total)[name])/float(n)*100,2)) + '%'
            content = [name,str(Counter(list_total)[name]),rna_percent]
            old_labels.append(name)
            old_sizes.append(int(Counter(list_total)[name]))
            f.write('\t'.join(content) + '\n')
    sizes,labels = sort_list(old_sizes,old_labels)
    total_colors = ['#ff9999','#66b3ff','#99ff99','#ffcc99','#c2c2f0',"#8b86bf", "#18aeda",
                    "#ba9333", "#f07570","#4575B4","#91BFDB","#E0F3F8","#FFFFBF","#FEE090","#FC8D59","#D73027"]
    colors = random.sample(total_colors, len(sizes))
    piechart_png = f'{indir}/{prefix}.RNA_type.piechart.png'
    piechart_pdf = f'{indir}/{prefix}.RNA_type.piechart.pdf'
    fig, ax = plt.subplots(figsize=(16, 12))
    fig.set_facecolor('white')
    wedges, texts, autotexts = ax.pie(
        sizes,
        colors=colors,
        startangle=90,
        wedgeprops={'edgecolor': 'black', 'linewidth': 1},
        autopct='',#autopct='%3.2f%%',
        pctdistance=0.85  # 百分比标签离中心的距离
    )
    plt.setp(autotexts, size=10, weight='bold', color='black')
    for i, (wedge, label) in enumerate(zip(wedges, labels)):
        ang = (wedge.theta2 - wedge.theta1) / 2. + wedge.theta1  # 扇形中心角度
        x = np.cos(np.deg2rad(ang))  # 角度转弧度
        y = np.sin(np.deg2rad(ang))
        horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
        connectionstyle = f"angle,angleA=0,angleB={ang}"
        label_distance = 1.3  # 标签距离中心的倍数
        label_x = label_distance * x
        label_y = label_distance * y
        ax.annotate(
            f"{label}: {sizes[i]} ({sizes[i] / sum(sizes):.2%})",  # 标签文本
            xy=(x * 0.8, y * 0.8),  # 引导线起点（扇形边缘内）
            xytext=(label_x, label_y),  # 标签位置
            horizontalalignment=horizontalalignment,
            verticalalignment='center',
            arrowprops=dict(
                arrowstyle="-",  # 直线箭头
                color="black",
                connectionstyle=connectionstyle,
                linewidth=1
            ),
            bbox=dict(
                boxstyle="round,pad=0.3",
                facecolor=colors[i],
                alpha=0.7,
                edgecolor='black'
            )
        )
    plt.title('Distribution of RNA types', fontsize=14, pad=30)
    plt.tight_layout()
    plt.savefig(piechart_png)
    plt.show()
    #plt.pie(rate_list, labels=kinds_list, autopct="%3.2f%%", startangle=150,pctdistance=0.9,labeldistance = 0.5, radius = 1.5)
    #plt.savefig(piechart_png)
    #plt.show()
    with PdfPages(piechart_pdf) as pdf:
        pdf.savefig()
    plt.show()
    plt.close()

def main():
    sys.stdout.write('%s 程序开始运行： -----------\n\n' % time.strftime('%Y-%m-%d %H:%M:%S'))
    argvs = parse_args()
    mkout(argvs.infile,argvs.info,argvs.outfile)
    sys.stdout.write('%s 程序顺利运行结束 !!!\n' % time.strftime('%Y-%m-%d %H:%M:%S'))

if __name__ == "__main__":
    main()