# coding: utf-8

''' 
 * All rights Reserved, Designed By HIT-Bioinformatics   
 * @Title:  Reads_Clustering.py
 * @Package: 
 * @Description: This is a clustering script for clustering minimap2's reads' overlaps results.
 * @author: Fu Yilei
 * @date: Augst 17 2018
 * @version V1.0.1     
'''

from collections import defaultdict
import json
import heapq
import time
import os
import sys
import multiprocessing
import logging
import argparse
# import progressbar

PASSED_LIST = [[]]
WHOLE_PASSED_LIST = list()
ALL_PASSED = list()
PASSED = [[]]
COUNT = 0
NAME_CLUSTER = list()
OVERLAP_LIST = multiprocessing.Manager().list()

#FILE_LIST = []


def execute_minimap_command(input_file_name, output_file_name, threads, seq_type):
    logging.info("No paf file found. Executing minimap2 overlap detecting commands...")
    os.system("minimap2 -x %s "%(seq_type)+input_file_name+" " +input_file_name+ " > "+output_file_name+" -t "+ str(threads))
    logging.info("Minimap2 command executed, paf file generated.")

def read_paf_file(file_name):
    #Open overlap file
    f = open(file_name)
    logging.info("Reading paf files...")
    file_list = f.readlines()
    f.close()
    logging.info("Finished")
    return file_list

def delete_useless_info(file_list, process_id):
    #Store overlap information, delete unless information
    length = len(file_list)
    logging.info("Process "+ str(process_id) +" is deleting usless information in the list...")
    for i in range(0, length):
        temp = file_list[i].split('\t')
        if (float(temp[9])/float(temp[10]) >= 0.9) and (2*temp[10] >= temp[6]):
        # if (float(temp[9])/float(temp[10]) >= 0.2):
            OVERLAP_LIST.append((temp[0], temp[5]))
    file_list = []
    logging.info("Process "+  str(process_id) +" is finished.")


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
           ovlpdic [i[0]].append(i[1])
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
'''

def generate_overlap_dic(ovlpl):
    logging.info("Generating overlap dictionary...")
    ovlpdic = {}
    length = len(ovlpl)
    count = 0
    # p = progressbar.ProgressBar(length)
    # p.start()
    for i in ovlpl:
        if i[0] in ovlpdic:
            ovlpdic[i[0]].append(i[1])
        else:
            ovlpdic.update({i[0]:[i[1]]})
        # p.update(count)
        count = count+1
    # p.finish()
    logging.info("Finished")
    return ovlpdic

'''    
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
            logging.info("Clustering reads based on read " + str(i))
            gene_search(overlap_dic, i)
            f.write(str(PASSED_LIST[COUNT])+"\n")
            #print(PASSED_LIST[COUNT],file = f)
            logging.info("Cluster " + str(COUNT) +" generated")
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


def preprocess_clustered_file(dict_name, cluster_count):
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
    for i in cluster:
        if len(i) >= cluster_count:
            temp = []
            for j in i:
                temp.append(">" + j + '\n')                                        
            NAME_CLUSTER.append(temp)
    logging.info("Deleted all clusters having reads larger than %d, the remaining clusters' number is: "%(cluster_count) + str(len(NAME_CLUSTER)))



def file_select(file_path, cluster_list, cluster_num):
    genome_size = 0
    with open(file_path+"unmapped.fasta", 'r') as f: 
        while True:
            line_0 = f.readline()
            if line_0 == '':
                break
            elif line_0 in cluster_list:
                fname = "clustered_reads" + str(cluster_num) + ".fa"                                    
                line_1 = f.readline()
                genome_size = genome_size + len(line_1)
                with open(file_path+fname, 'a') as clusterfile:
                    clusterfile.write(line_0)
                    clusterfile.write(line_1)

    logging.info("Cluster file " + str(cluster_num) + " generated, genome size is:" +str(genome_size))
    with open(file_path+"genesize.log", 'a+') as gensizefile:
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


def execute_file_select(file_path, divided_cluster, block_num, process_num):
    i = divided_cluster[process_num]
    for j in range(len(i)):
        logging.info("Process No." + str(process_num) + " is generating clustered read " + str(j + block_num))
        file_select(file_path, i[j], j + block_num)

def parseArgs(argv):
	parser = argparse.ArgumentParser(prog="Reads_Clustering.py", description=USAGE, formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-f', '--fa', help = "fa file path.",  type = str)
	args = parser.parse_args(argv)
	return args


def cluster_unaligned_reads(input, threads, seq_type, cluster_count):
    if not os.path.exists(input+"overlaps.paf"):
        execute_minimap_command(input+"unmapped.fasta", input+"overlaps.paf", threads, seq_type)
    if not os.path.exists(input+"olapdic.log"):
        file_list = read_paf_file(input+"overlaps.paf")
        file_list_length = int(len(file_list)/threads)
        process_list = list()
        for i in range(threads-1):
            process_list.append(multiprocessing.Process(target=delete_useless_info, args=(file_list[file_list_length*i:file_list_length*(i+1)],i)))
            #print("Deleting process "+str(i)+" added into pool.")
        for i in range(threads-1):
            process_list[i].start()
        for i in range(threads-1):
            process_list[i].join()
        logging.info("All unaligned clusters Finished!")

        d = generate_overlap_dic(OVERLAP_LIST)
        with open(input+"dictionary.dict", 'w+') as f:
            f.write(str(d))
        #generate_passed_list(d, overlap_info_file_name)
        generate_clustered_file(d, input+"olapdic.log")
        # end_1 = time.time()
    else:
        logging.info("The cluster file has already been generated, producing clusterd reads now...")
        
    preprocess_clustered_file(input+"olapdic.log", cluster_count)

    list_length = int(len(NAME_CLUSTER)/threads)
    d_cluster = thread_division(NAME_CLUSTER, threads, list_length)
    for i in range(threads):
        multiprocessing.Process(target=execute_file_select, args=(input, d_cluster, list_length*i, i)).start()


def run(argv):
    args = parseArgs(argv)

    paf_file_name = "overlaps.paf"                                                                              #NOTICE! This should be modified accroding to input file
    overlap_info_file_name = "olapdic.log"                                                                      #This is a temp file.
    original_read_file_name = args.fa
    #original_read_file_name = '/data/tjiang/rCANID/real_reads/na12878_pacbio_unmapped.fa'                      #NOTICE! This should be modified accroding to input file
    minimap2_thread_num = 32                                                                                    #NOTICE! This should be modified accroding to different machines.
    process_num = 32                                                                                            #NOTICE! This should be modified accroding to different machines.  Thread num must be larger than 2.
    if not os.path.exists(paf_file_name):                                                                       #if there is no overlap information, execute minimap2 command
        execute_minimap_command(original_read_file_name, paf_file_name, minimap2_thread_num, "ava-pb")                    
    if not os.path.exists(overlap_info_file_name):
        start = time.time()
        file_list = read_paf_file(paf_file_name)
        file_list_length = int(len(file_list)/process_num)
        process_list = []
        for i in range(process_num-1):
            process_list.append(multiprocessing.Process(target=delete_useless_info, args=(file_list[file_list_length*i:file_list_length*(i+1)],i)))
            #print("Deleting process "+str(i)+" added into pool.")
        for i in range(process_num-1):
            process_list[i].start()
        for i in range(process_num-1):
            process_list[i].join()
        logging.info("All Finished!")
        d = generate_overlap_dic(OVERLAP_LIST)
        with open("dictionary.dict", 'w+') as f:
            f.write(str(d))
        #generate_passed_list(d, overlap_info_file_name)
        generate_clustered_file(d, overlap_info_file_name)
        end_1 = time.time()
        logging.info("Cluster generation time used: " + str(end_1-start))
    else:
        logging.info("The cluster file have already been generated, producing clusterd reads now...")
    os.system("mkdir contigs")
    preprocess_clustered_file(overlap_info_file_name)                                                #preprocess clustered files
    
    list_length = int(len(NAME_CLUSTER)/process_num)
    d_cluster = thread_division(NAME_CLUSTER, process_num, list_length)
    # print "here"
    for i in range(process_num):
        multiprocessing.Process(target=execute_file_select, args=(original_read_file_name, d_cluster, list_length*i, i)).start()

if __name__ == "__main__":
    run(sys.argv[1:])

    '''
    paf_file_name = "overlaps.paf"                                                                              #NOTICE! This should be modified accroding to input file
    overlap_info_file_name = "olapdic.log"                                                                      #This is a temp file.
    original_read_file_name = '/data/tjiang/rCANID/real_reads/na12878_pacbio_unmapped.fa'                       #NOTICE! This should be modified accroding to input file
    minimap2_thread_num = 32                                                                                    #NOTICE! This should be modified accroding to different machines.
    process_num = 32                                                                                            #NOTICE! This should be modified accroding to different machines.  Thread num must be larger than 2.
    if not os.path.exists(paf_file_name):                                                                       #if there is no overlap information, execute minimap2 command
        execute_minimap_command(original_read_file_name, paf_file_name, minimap2_thread_num)                    

    if not os.path.exists(overlap_info_file_name):
        start = time.time()
        file_list = read_paf_file(paf_file_name)
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
        d = generate_overlap_dic(OVERLAP_LIST)
        with open("dictionary.dict", 'w+') as f:
            f.write(str(d))
        #generate_passed_list(d, overlap_info_file_name)
        generate_clustered_file(d, overlap_info_file_name)
        end_1 = time.time()
        print("Cluster generation time used: " + str(end_1-start))
    else:
        print("The cluster file have already been generated, producing clusterd reads now...")
    preprocess_clustered_file(overlap_info_file_name)                                                #preprocess clustered files
    list_length = int(len(NAME_CLUSTER)/process_num)
    d_cluster = thread_division(NAME_CLUSTER, process_num, list_length)

    for i in range(process_num):
        multiprocessing.Process(target=execute_file_select, args=(original_read_file_name, d_cluster, list_length*i, i)).start()
    '''