# coding: utf-8
#This is a clustering script for clustering minimap2's reads' overlaps results.
#Author: Yilei Fu

from collections import defaultdict
import json
import heapq
import time
import os
import sys
import multiprocessing
#import progressbar

PASSED_LIST = [[]]
WHOLE_PASSED_LIST = []
ALL_PASSED = []
PASSED = [[]]
COUNT = 0
NAME_CLUSTER = []
OVERLAP_LIST = multiprocessing.Manager().list()

#FILE_LIST = []


def execute_minimap_command(input_file_name, output_file_name, threads):
    print("No paf file found. Executing minimap2 overlap detecting commands...")
    os.system("minimap2 -x ava-pb "+input_file_name+" " +input_file_name+ " > "+output_file_name+" -t "+ str(threads))
    print("Minimap2 command executed, paf file generated.")

def read_paf_file(file_name):
    #Open overlap file
    f = open(file_name)
    print("Reading paf files...")
    file_list = f.readlines()
    f.close()
    print("Finished")
    return file_list

def delete_useless_info(file_list, process_id):
    #Store overlap information, delete unless information
    length = len(file_list)
    print("Process "+ str(process_id) +" is deleting usless information in the list...")
    for i in range(0, length):
        temp = file_list[i].split('\t')
        if float(temp[9])/float(temp[10]) >= 0.9:
            OVERLAP_LIST.append((temp[0], temp[5]))
    file_list = []
    print("Process "+  str(process_id) +" is finished.")

'''
def key_position(dic, key):
    return dic.keys().index(key)
'''

def generate_overlap_dic():
    #Use dictionary to store overlap information
    print("Generating overlap dictionary...")
    ovlpdic = {}
    cnt = 0
    for i in OVERLAP_LIST:
        if (i[0] not in ALL_PASSED) and (i[1] not in ALL_PASSED):
            ovlpdic.update({i[0]:[i[1]]})
            ALL_PASSED.append(i[0])
            ALL_PASSED.append(i[1])
            print("A new cluster "+ str(cnt) +" appeared.")
            cnt = cnt+1
        elif i[0] in ovlpdic:
            ovlpdic[i[0]].append(i[1])
            ALL_PASSED.append(i[1])
        elif i[1] in ovlpdic:
            ovlpdic[i[1]].append(i[0])
            ALL_PASSED.append(i[0])
        elif i[0] in ALL_PASSED:
            for key in ovlpdic:
                if i[0] in ovlpdic[key]:
                    ovlpdic[key].append(i[1])
                    ALL_PASSED.append(i[1])
        else:
            for key in ovlpdic:
                if i[1] in ovlpdic[key]:
                    ovlpdic[key].append(i[0])
                    ALL_PASSED.append(i[0])
    print("Finished")
    return ovlpdic

def generate_passed_list(info_dic, dict_name):
    f = open(dict_name, 'w+')
    temp = []
    for key in info_dic:
        temp.append(key)
        temp = temp+info_dic[key]
        #print(temp)
        f.write(str(temp)+"\n")
        temp = []
'''
def generate_clustered_file(overlap_dic, dict_name):
    #Generate Clustered file
    global COUNT
    f = open(dict_name, 'w+')
    for i in overlap_dic:
        if (i not in WHOLE_PASSED_LIST):
            print("Clustering reads based on read " + str(i))
            gene_search(overlap_dic, i)
            f.write(str(PASSED_LIST[COUNT])+"\n")
            #print(PASSED_LIST[COUNT],file = f)
            print("Cluster " + str(COUNT) +" generated")
            COUNT = COUNT+1
            PASSED_LIST.append([])
    f.close()

    def gene_search(dic, node):
        #recursively clustering overlap files to get cluster info
        global COUNT
        if node not in WHOLE_PASSED_LIST:
            PASSED_LIST[COUNT].append(node)
            WHOLE_PASSED_LIST.append(node)
            if node in dic:
                for i in dic[node]:
                    gene_search(dic, i)
'''

def preprocess_clustered_file(dict_name, file_type):
    #Dividing the whole reads file to divide whole reads into different clustered reads
    clusterfile = open(dict_name)
    clusterlist = clusterfile.readlines()
    clusterfile.close()
    cluster = []
    global NAME_CLUSTER
    for i in range(len(clusterlist)):
        templist = []
        splited_cluster_list = clusterlist[i].split("'")
        for j in range(len(splited_cluster_list)):
            if j%2 == 1:
                templist.append(splited_cluster_list[j])
        cluster.append(templist)
    if file_type == "fastq" or file_type == "fq":
        for i in cluster:
            if len(i) >= 10:
                temp = []
                for j in i:
                    temp.append("@" + j + '\n')                                                                 
                NAME_CLUSTER.append(temp)
    else:
        for i in cluster:
            if len(i) >= 10:
                temp = []
                for j in i:
                    temp.append(">" + j + '\n')                                        
                NAME_CLUSTER.append(temp)



def file_select(file_name, cluster_list, cluster_num, file_type):
    genome_size = 0
    if file_type == "fastq" or file_type == "fq":
        with open(file_name, 'r') as f: 
            while True:
                line_0 = f.readline()
                #print(line_0[0:5])
                if line_0 == '':
                    break
                elif line_0 in cluster_list:
                    #print(line_0)
                    fname = "clustered_reads" + str(cluster_num) + "." + file_type                                     
                    line_1 = f.readline()
                    genome_size = genome_size + len(line_1)
                    line_2 = f.readline()
                    line_3 = f.readline()
                    with open(fname, 'a') as clusterfile:
                        clusterfile.write(line_0)
                        clusterfile.write(line_1)
                        clusterfile.write(line_2)
                        clusterfile.write(line_3)
    else:
        with open(file_name, 'r') as f: 
            while True:
                line_0 = f.readline()
                if line_0 == '':
                    break
                elif line_0 in cluster_list:
                    fname = "clustered_reads" + str(cluster_num) + "." + file_type                                    
                    line_1 = f.readline()
                    genome_size = genome_size + len(line_1)
                    with open(fname, 'a') as clusterfile:
                        clusterfile.write(line_0)
                        clusterfile.write(line_1)

    print("Cluster file " + str(cluster_num) + " generated, genome size is:" +str(genome_size))
    with open("genesize.log", 'a+') as gensizefile:
        gensizefile.write(str(genome_size)+":"+str(cluster_num)+"\n")
    #gsizelist.append(str(genomesize)+":"+(cluster_num))



def thread_division(whole_cluster, process_num, list_num):
    divided_cluster = []
    #print(Thread_num)
    for i in range(process_num-1):
        divided_cluster.append(whole_cluster[(list_num*i):(list_num*(i+1))])
        #print(0)
    divided_cluster.append(whole_cluster[list_num*(process_num-1):])
    return divided_cluster


def execute_file_select(file_name, file_type, divided_cluster, block_num, process_num):
    i = divided_cluster[process_num]
    for j in range(len(i)):
        print("Process NO." + str(process_num) + " is generating clustered read " + str(j + block_num) + "\n")
        file_select(file_name, i[j], j + block_num, file_type)

def file_property_judge(input_file):
    if input_file[-5:] == 'fastq':
        print("Input file is: " + input_file +", this is a fastq file.")
        return "fastq"
    elif input_file[-5:] == 'fasta':
        print("Input file is: " + input_file +", this is a fasta file.")
        return "fasta"
    elif input_file[-2:] == 'fa':
        print("Input file is: " + input_file +", this is a fa file.")
        return "fa"
    elif input_file[-2:] == 'fq':
        print("Input file is: " + input_file +", this is a fq file.")
        return "fq"
    else:
        print("File is not neither fastq format nor fasta format, please open this program and change the file name!")
        sys.exit(0) 

if __name__ == "__main__":
    paf_file_name = "overlaps.paf"                                                                              #NOTICE! This should be modified accroding to input file
    overlap_info_file_name = "olapdic.log"                                                                      #This is a temp file.
    original_read_file_name = '/data/tjiang/dbNSF/real_reads/na12878_pacbio_unmapped.fasta'                     #NOTICE! This should be modified accroding to input file
    minimap2_thread_num = 32                                                                                    #NOTICE! This should be modified accroding to different machines.
    process_num = 32                                                                                            #NOTICE! This should be modified accroding to different machines.  Thread num must be larger than 2.

    file_type = file_property_judge(original_read_file_name)                                                    #detect file type

    if not os.path.exists(paf_file_name):                                                                       #if there is no overlap information, execute minimap2 command
        execute_minimap_command(original_read_file_name, paf_file_name, minimap2_thread_num)                    

    if not os.path.exists(overlap_info_file_name):
        start = time.time()
        #generate_clustered_file(generate_overlap_dic(delete_useless_info(read_paf_file(paf_file_name))), overlap_info_file_name)
        file_list = read_paf_file(paf_file_name)
        #output_list = multiprocessing.Array("i", range(len(file_list)))
        file_list_length = int(len(file_list)/process_num)
        process_list = []
        for i in range(process_num-1):
            process_list.append(multiprocessing.Process(target=delete_useless_info, args=(file_list[file_list_length*i:file_list_length*(i+1)],i)))
            #print("Deleting process "+str(i)+" added into pool.")
        for i in range(process_num-1):
            process_list[i].start()
        for i in range(process_num-1):
            process_list[i].join()
        print("All Finished!")
        d = generate_overlap_dic()
        with open("dictionary.dict", 'w+') as f:
            f.write(str(d))
        generate_passed_list(d, overlap_info_file_name)
        end_1 = time.time()
        print("Cluster generation time used: " + str(end_1-start))
    else:
        print("The cluster file have already been generated, producing clusterd reads now...")
    preprocess_clustered_file(overlap_info_file_name, file_type)                                                #preprocess clustered files
    list_length = int(len(NAME_CLUSTER)/process_num)
    d_cluster = thread_division(NAME_CLUSTER, process_num, list_length)

    for i in range(process_num):
        multiprocessing.Process(target=execute_file_select, args=(original_read_file_name, file_type, d_cluster, list_length*i, i)).start()
