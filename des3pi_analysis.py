#!/usr/bin/env python

import os
import math
from operator import itemgetter


start_path = os.getcwd()
os.chdir(start_path)


####### IDENTIFY THE IDENTICAL HOT SPOTS GROUPS
peptide_centroid = []
maximum_groups = 0

dico_clusters = {}

def compute_distance_cluster(x1,y1,z1,x2,y2,z2):
    distance = math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
    return distance

for i in range(1,21):
    list_centres = []
    if i <= 9:
        os.chdir("0{}/docking_results".format(i))
    else: 
        os.chdir("{}/docking_results".format(i))
    
    try:
        fichier = open("peptide_centroid.txt","r")
        fichier.close()
        fichier = open("cluster_order.txt","r")
        fichier.close()
    except:
        os.chdir(start_path)
        continue
    with open('peptide_centroid.txt','r') as centroid:
        centres = centroid.readlines()
        for centre in centres:
            list_centres.append(float(centre))

    if maximum_groups >= len(list_centres):
        pass
    else:
        maximum_groups = len(list_centres)
    peptide_centroid.append(list_centres)
    os.chdir(start_path)
index_to_test = []
for i in range(0,maximum_groups):
    for j,peptide in enumerate(peptide_centroid):
        try:
            peptide[i]
            index_to_test.append(j)
            break
        except:
            continue    
    
dico_groups = {}
for i in range(0,maximum_groups):
    dico_groups[i+1] = []

cb = 0
for groups, index in enumerate(index_to_test):
    for k,peptide in enumerate(peptide_centroid):
        for l,centre in enumerate(peptide):
            if math.isclose(centre,peptide_centroid[index][groups],abs_tol=4) == True:
                cb = cb +1
                dico_groups[groups+1].append((k+1,l+1))

##### END IDENTIFYING GROUPS

#### GET CLUSTERS

os.chdir(start_path)
dico_clusters = {}
for i in range(1,21):
    clusters_for_a_run = []
    if i <= 9:
        os.chdir("0{}/docking_results".format(i))
    else: 
        os.chdir("{}/docking_results".format(i))

    try:
        fichier = open("cluster_order.txt","r")
        fichier.close()
    except:
        continue

    with open('cluster_order.txt','r') as order:
        lines = order.readlines()
    for line in lines:
        line2 = line.split("+")
        for group in line2:
            single_group = group.split("cluster")[1:]
            if single_group == []:
                continue
            else:
                clusters_for_a_run.append(single_group)
    dico_clusters[i] = clusters_for_a_run
    os.chdir(start_path)

##### END GET CLUSTERS
##### REGROUP IN A GROUP BY LENGTH
dico_groups_length = {}
lengths = []
for group in dico_groups.keys():
    dico_length = {}
    runs = dico_groups[group]
    for run in runs:
        cluster = dico_clusters[run[0]][run[1]-1]
        lengths.append(len(cluster))

    list_lengths_unique = list(set(lengths))
    for length in list_lengths_unique:
        list_length = []
        for run in runs:
            cluster = dico_clusters[run[0]][run[1]-1]
            if length == len(cluster):
                list_length.append((run[0], cluster))
                
        dico_length[length] = list_length
    dico_groups_length[group] = dico_length

    
##### END REGROUPING BY LENGTH

#### COMPUTE IDENTICAL CLUSTERS...

def store_coordinates(run,cluster,start_path=start_path):
    if run <=9:
        run = "0" + str(run)
    coordinates_1 = []
    with open(start_path + "/{}/".format(run) + "docking_results/" + "cluster{}.xyz".format(cluster)) as cluster1:
        str_coord = cluster1.readlines()[2][2:].split()
        for i in str_coord:
            coordinates_1.append(float(i))
    return coordinates_1



def test_is_identical(clusters1,start_path=start_path,run1 = 6, run2=13,clusters2=0):
    nb_pairs = 0
    correct_clusters1 = []
    correct_clusters2 = []
    for i in clusters1:
        correct_clusters1.append(int(i))
    for j in clusters2:
        correct_clusters2.append(int(j))   

    coord_clst_1 = store_coordinates(run1,correct_clusters1[0])
    coord_run_1 = []
    for i in correct_clusters1:
        coord_run_1.append(store_coordinates(run1,i))
    coord_run_2 = []
    for i in correct_clusters2:
        coord_run_2.append(store_coordinates(run2,i))

    for j, coord1 in enumerate(coord_run_1):
        for i,coord2 in enumerate(coord_run_2):
            if math.sqrt((coord2[0]-coord1[0])**2+(coord2[1]-coord1[1])**2+(coord2[2]-coord1[2])**2) <= 1.75:
                nb_pairs = nb_pairs + 1 
                #print("{} and {} are identical".format(correct_clusters1[j], correct_clusters2[i]))
                break
            else:
                continue   
    if nb_pairs == len(clusters1) and nb_pairs == len(clusters2):
        return True
    else:
        return False
   
    
dico_unic_runs = {}
stop = False
for group_index in dico_groups_length.keys():
    group = dico_groups_length[group_index]
    for length_index in group.keys(): 
        length = group[length_index]
        runs_to_test = []
        for run_hotspots in length:
            runs_to_test.append(run_hotspots[0])
        for clst,run_to_test in enumerate(runs_to_test):
            for i in dico_unic_runs.keys():
                if length[clst] in dico_unic_runs[i]:
                    stop = True
                    break
                else:
                    stop = False
            if stop == False:
                pass
            else:
                continue
            liste_tempo = []
            liste_tempo.append(length[clst])
            if runs_to_test[clst+1:] != []:
                for clst2,tested in enumerate(runs_to_test[clst+1:]):
                    test_results = test_is_identical(clusters1=length[clst][1], clusters2=length[clst+clst2+1][1],run1 = run_to_test, run2 = tested)
                    if test_results == True:
                        liste_tempo.append(length[clst+clst2+1])
            try:
                dico_except = dico_unic_runs[0]
                new = list(dico_unic_runs.keys())[-1]+1
                dico_unic_runs[new] = liste_tempo
            except:
                dico_unic_runs[0] = liste_tempo
for i in dico_unic_runs.keys():
    for j in dico_unic_runs[i]:
        print(j)
    print("END GROUP")
#### END COMPUTING IDENTICAL CLUSTERS
#### COMPUTING SEQUENCES...


def identical_sequences(clusters1,start_path=start_path,run1 = 6, run2=13,clusters2=0):
    dict_sequences = {}
    nb_pairs = 0
    correct_clusters1 = []
    correct_clusters2 = []
    for i in clusters1:
        correct_clusters1.append(int(i))
    for j in clusters2:
        correct_clusters2.append(int(j))   

    coord_clst_1 = store_coordinates(run1,correct_clusters1[0])
    coord_run_1 = []
    for i in correct_clusters1:
        coord_run_1.append(store_coordinates(run1,i))
    coord_run_2 = []
    for i in correct_clusters2:
        coord_run_2.append(store_coordinates(run2,i))
    for j, coord1 in enumerate(coord_run_1):
        dict_sequences[j] = []
        for i,coord2 in enumerate(coord_run_2):
            if math.sqrt((coord2[0]-coord1[0])**2+(coord2[1]-coord1[1])**2+(coord2[2]-coord1[2])**2) <= 1.75:
                nb_pairs = nb_pairs + 1 
                #print("{} and {} are identical".format(correct_clusters1[j], correct_clusters2[i]))
                if run2 <=9:
                    accurate_run2 = "0" + str(run2)
                else:
                    accurate_run2 = run2
                with open(start_path + "/{}/docking_results/".format(accurate_run2) + "best_aa_per_cluster.txt") as best_aa:
                    lines = best_aa.readlines()
                    for line in lines:
                        if "cluster {}:".format(int(correct_clusters2[i])) in line:
                            dict_sequences[j].append(line[-4] + line[-3] + line[-2])
                break
            else:
                continue   


    return dict_sequences

dict_sequences = {}
for key in dico_unic_runs.keys():
    dict_sequences_1_group = {}
    ref = dico_unic_runs[key][0]
    dict_sequences_1_group = identical_sequences(ref[1],run1=ref[0],run2=ref[0],clusters2=ref[1])
    for run in dico_unic_runs[key][1:]:
        dict_sequences_1_run = identical_sequences(ref[1],run1=ref[0],run2=run[0],clusters2=run[1])
        for cluster in dict_sequences_1_group.keys():
            dict_sequences_1_group[cluster] = dict_sequences_1_group[cluster] + dict_sequences_1_run[cluster]
    dict_sequences[key] = dict_sequences_1_group

print(dict_sequences)

##### END COMPUTING SEQUENCES


##### WRITING RESULTS IN CORRECT FILES

## Writing cluster position

for group in dico_unic_runs.keys():
    os.chdir(start_path)
    os.system("mkdir class{}".format(group))
    run = dico_unic_runs[group][0][0]
    if run <=9:
        run_to_write = "0" + str(run)
    else:
        run_to_write = str(run)
    clusters_to_write = dico_unic_runs[group][0][1]
    for i,cluster in enumerate(clusters_to_write):
        os.system("cp {}/{}/docking_results/cluster{}.xyz {}/class{}/cluster{}.xyz".format(start_path,run_to_write,int(cluster),start_path,group,i+1))

## End writing cluster position

## Writing sequences per run

aa_library = {"ala":"A","asn":"N","gln":"Q","ile":"I","lys":"K","phe":"F","thr":"T","arg":"R","cys":"C","gly":"G","leu":"L","met":"M","ser":"S","pro":"P","asp":"D","glu":"E","his":"H","val":"V","tyr":"Y","trp":"W"}


length_peptide = []
with open("{}/des3pi_all_sequences.txt".format(start_path),"w") as final_result:
    for groups in dict_sequences.keys():
        final_result.write("class {}:\n".format(groups))
        with open("{}/class{}/peptides.txt".format(start_path,groups),"w") as peptides:
            for i in range(0,len(dict_sequences[groups][0])):
                for cluster in dict_sequences[groups]:
                    peptides.write(aa_library[dict_sequences[groups][cluster][i]])
                peptides.write("\n")
        with open("{}/class{}/best_5_peptides.txt".format(start_path,groups),"w") as best_peptides:
            with open("{}/class{}/peptides.txt".format(start_path,groups),"r") as peptides:
                lines = peptides.readlines()

                clean_lines = []
                for i in lines:
                    clean = i.split("\n")[0]
                    clean_str = ""
                    for j in clean:
                        if j.isupper() == True:
                            clean_str = clean_str + j
                    clean_lines.append(clean_str)

                counts = {}
                global_counts ={}
                for i in range(1,len(clean_lines[0])+1):
                    counts[i] = []
                    global_counts[i] = ""

                for i in clean_lines:
                    for position,j in enumerate(i):
                        counts[position+1].append(j)

                import itertools
                list_counts = []
                for i in counts.keys():
                    count_aa = []
                    test = list(set(counts[i]))
                    for j in test:
                        aa_count = j,counts[i].count(j)
                        count_aa.append(aa_count)
                        if j not in global_counts[i]:
                            global_counts[i] = global_counts[i] + j
                    list_counts.append(count_aa)

                nb_hotspots = len(lines[0].split("\n")[0])
                print(nb_hotspots)
                if str(nb_hotspots) == "4":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4])) 
                if str(nb_hotspots) == "5":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5]))
                if str(nb_hotspots) == "6":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5],global_counts[6]))
                if str(nb_hotspots) == "7":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5],global_counts[6],global_counts[7]))
                if str(nb_hotspots) == "8":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5],global_counts[6],global_counts[7],global_counts[8]))
                if str(nb_hotspots) == "9":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5],global_counts[6],global_counts[7],global_counts[8],global_counts[9]))
                if str(nb_hotspots) == "10":
                    list_combinations = list(itertools.product(global_counts[1],global_counts[2],global_counts[3],global_counts[4],global_counts[5],global_counts[6],global_counts[7],global_counts[8],global_counts[9],global_counts[10]))

                list_list_combinations = []
                for i in list_combinations:
                    list_list_combinations.append(list(i))


                list_seq_ordered = []
                for i in list_list_combinations:
                    count_seq = 0
                    for j,k in enumerate(i):
                        for element in list_counts[j]:
                            if element[0] == k:
                                count_seq = count_seq + element[1]
                                break
                    list_seq_ordered.append(("".join(i),count_seq))

                list_sorted = sorted(list_seq_ordered, key=itemgetter(1), reverse = True)
                best_5_peptides = list_sorted[0:5]
                print(best_5_peptides)
                for best_pep in best_5_peptides:
                    str_pep = ""
                    for i in range(0,nb_hotspots):
                        str_pep = str_pep + str(best_pep[0][i])
                        if i == nb_hotspots -1:
                            continue
                        else:
                            clst_neighbour = i+2
                            clst1 = open("{}/class{}/cluster{}.xyz".format(start_path,groups,str(i+1),"r"))
                            clst2 = open("{}/class{}/cluster{}.xyz".format(start_path,groups,str(clst_neighbour),"r"))
                            lines1 = clst1.readlines()
                            lines2 = clst2.readlines()
                            x1,y1,z1 = float(lines1[2].split()[1]),float(lines1[2].split()[2]),float(lines1[2].split()[3])
                            x2,y2,z2 = float(lines2[2].split()[1]),float(lines2[2].split()[2]),float(lines2[2].split()[3])
                            distance = compute_distance_cluster(x1,y1,z1,x2,y2,z2)
                            if distance >=4.3:
                                str_pep = str_pep + math.ceil(distance/4.3 -1)*'g'
                    first_clst = 1
                    clst1 = open("{}/class{}/cluster{}.xyz".format(start_path,groups,str(nb_hotspots),"r"))
                    clst2 = open("{}/class{}/cluster{}.xyz".format(start_path,groups,str(first_clst),"r"))
                    lines1 = clst1.readlines()
                    lines2 = clst2.readlines()
                    x1,y1,z1 = float(lines1[2].split()[1]),float(lines1[2].split()[2]),float(lines1[2].split()[3])
                    x2,y2,z2 = float(lines2[2].split()[1]),float(lines2[2].split()[2]),float(lines2[2].split()[3])
                    distance = compute_distance_cluster(x1,y1,z1,x2,y2,z2)
                    if distance >=4.3:
                        str_pep = str_pep + math.ceil(distance/4.3 -1)*'g' 
                    best_peptides.write(str_pep + "\n")
                    final_result.write(str_pep + "\n")
    

                        
    
    
