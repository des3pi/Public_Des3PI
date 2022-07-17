#!/usr/bin/env python

import math
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import sklearn
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def compute_distance_cluster(x1,y1,z1,x2,y2,z2):
    distance = math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2))
    return distance

def compute_angle_bac(dist_bc,dist_ab,dist_ac):
    angle = math.acos((dist_ab*dist_ab+dist_bc*dist_bc-dist_ac*dist_ac)/(2*dist_ab*dist_bc))*180/math.pi
    return angle

def compute_more_represented_aa(poses_df,df_clusters,nb_dock):

    # Compute the more represented amino-acids for each cluster
    #print(poses_df['cluster'].value_counts().max())
    perc_cluster = pd.DataFrame()
    perc_cluster['aa'] = list(poses_df['aa'].unique())
    percent_aa = []
    perc_cluster = perc_cluster.set_index('aa')

    for j in list(df_clusters.index):
        listpercent = []
        clust_interest= int(j.split(" ")[1])
        for i in range(0,nb_dock*20*9,nb_dock*9):
            df_percent = poses_df.iloc[i:i+nb_dock*9,:]
            df_percent_cluster = df_percent.loc[df_percent["cluster"]==clust_interest]
            nb_aa_cluster = len(list(df_percent_cluster["cluster"]))
            listpercent.append((nb_aa_cluster))
        perc_cluster["cluster {}".format(clust_interest)] = listpercent

    aa_per_cluster = {}
    for cluster in list(perc_cluster):
        criteria = perc_cluster[cluster].mean() + perc_cluster[cluster].std()
        transi_df = perc_cluster[cluster].loc[perc_cluster[cluster]>criteria]
        aa_per_cluster[cluster] = list(transi_df.index)

    # Store the best of them in a dictionary

    aa_per_cluster = {}
    for cluster in list(perc_cluster):
        criteria = perc_cluster[cluster].mean() + perc_cluster[cluster].std()
        transi_df = perc_cluster[cluster].loc[perc_cluster[cluster]>criteria]
        aa_per_cluster[cluster] = list(transi_df.index)

    best_aa_per_cluster = {}
    best_aa_nb_poses = []
    for cluster in list(perc_cluster):
        criteria = perc_cluster[cluster].max()
        transi_df = perc_cluster[cluster].loc[perc_cluster[cluster]==criteria]
    
        best_aa_per_cluster[cluster] = list(transi_df.index)
        best_aa_nb_poses.append([cluster,criteria])
    print(best_aa_per_cluster)
    return best_aa_per_cluster, best_aa_nb_poses

def compute_mean_system(df_clusters_group):
    df_clusters_group_processed = df_clusters_group.drop(['cluster'], axis=1)
    df_2D = pd.DataFrame()
    cluster =[]
    x_ca =[]
    y_ca =[]
    acp = PCA(svd_solver='full')
    
    coord = acp.fit_transform(df_clusters_group_processed)

    for i in range(0,len(list(df_clusters_group.index))):
        cluster.append(df_clusters_group.index[i])
        x_ca.append(coord[i,0])
        y_ca.append(coord[i,1])

    df_2D["cluster"] = cluster
    df_2D["x_ca"] = x_ca
    df_2D["y_ca"] = y_ca 
    return df_2D

def get_optimal_order_cluster(results_df,df_clusters):
    cluster_order = {}
    nb_big_group = df_clusters['cluster'].max()
    for main in range(1, nb_big_group+1):
        df_clusters_group = df_clusters.loc[df_clusters["cluster"] == main]
        if len(df_clusters_group["cluster"]) <= 3:
            continue 
        index_list = list(df_clusters_group.index)
        df_clusters_group = compute_mean_system(df_clusters_group)
        df_clusters_group["cluster"] = index_list
        df_clusters_group = df_clusters_group.set_index("cluster")
        print(df_clusters_group.iloc[0])
        points = {}
        points_new_dim = {}
        polar_coor = {}

        # DEFINING THE FIRST POINT/ INITIALIZATION
        x_min = df_clusters_group["x_ca"].min()
        x_max = df_clusters_group["x_ca"].max()
        y_min = df_clusters_group["y_ca"].min()
        y_max = df_clusters_group["y_ca"].max()
        distance_x = math.sqrt((x_min-x_max)**2)
        distance_y = math.sqrt((y_min-y_max)**2)
        first_point = index_list[0]
        print(index_list)
        x_0 = float(df_clusters_group.loc[first_point]["x_ca"])
        y_0 = float(df_clusters_group.loc[first_point]["y_ca"])

        if distance_x >= distance_y:
            angle_x, angle_y = 0,0
        else:
            angle_x, angle_y = math.pi/2, math.pi/2

        # DF without the first point
        df_wo_pt1 = df_clusters_group.drop(index = [first_point])
        cluster_sequence = [first_point]


        # STORE COORDINATES
        for i in list(df_wo_pt1.index):
            x = float(df_wo_pt1.loc[i]["x_ca"])
            y = float(df_wo_pt1.loc[i]["y_ca"])
            points[i] = (x,y)

        # COMPUTE COORDINATES IN THE NEW SYSTEM,
        # using translation + rotation of vectors
        # x' = (x-x0)cos(theta) + (y-y0)sin(theta) ; y' = -(x-x0)sin(theta)+(y-y0)cos(theta)

        for i in list(points.keys()):
            x = (points[i][0] - x_0)*math.cos(angle_x) + (points[i][1] - y_0)*math.sin(angle_y)
            y = -(points[i][0] - x_0)*math.sin(angle_x) + (points[i][1] - y_0)*math.cos(angle_y)
            points_new_dim[i] = (x,y)
    
        # COMPUTE POLAR COORDINATES

        for i in list(points_new_dim.keys()):
            r = math.sqrt(points_new_dim[i][0]**2+points_new_dim[i][1]**2)
            if points_new_dim[i][0] > 0:
                theta = math.atan(points_new_dim[i][1]/points_new_dim[i][0])
            elif points_new_dim[i][0] < 0 and points_new_dim[i][1] >= 0:
                theta = math.atan((points_new_dim[i][1]/points_new_dim[i][0])+math.pi)
            elif points_new_dim[i][0] < 0 and points_new_dim[i][1] <0:
                theta = math.atan((points_new_dim[i][1]/points_new_dim[i][0])-math.pi)

            elif points_new_dim[i][0] == 0 and points_new_dim[i][1] >0:
                theta = math.pi/2
            elif points_new_dim[i][0] == 0 and points_new_dim[i][1] <0:
                theta = -math.pi/2
    
            polar_coor[i] = (r,theta)
        
        list_sorted_polar=[]      
        for i in list(polar_coor.keys()):
            list_sorted_polar.append([i,polar_coor[i]])

    
        list_sorted_polar = sorted(list_sorted_polar, key = lambda list_sorted_polar :(list_sorted_polar[1][1],list_sorted_polar[1][0]))
        first_theta = list_sorted_polar[0][1][1]
        last_theta = list_sorted_polar[-1][1][1]

        new_list_sorted_polar = []
        for i,k in enumerate(list_sorted_polar):
            if k == list_sorted_polar[-1]:
                new_list_sorted_polar.append(k)
                break
            j = i
            current_theta = k[1][1]
            next_theta = list_sorted_polar[j][1][1]
            while math.isclose(current_theta,next_theta,abs_tol=0.5) == True:
                current_theta = next_theta
                j = j +1
                if j > len(list_sorted_polar)-1:
                    break
                else:
                    next_theta = list_sorted_polar[j][1][1]
            new_list_sorted_polar.append([k[0], (k[1][0],current_theta)])

        new_list_sorted_polar = sorted(new_list_sorted_polar, key = lambda new_list_sorted_polar :(new_list_sorted_polar[1][1],new_list_sorted_polar[1][0]))      

        for i in new_list_sorted_polar:
            cluster_sequence.append(i[0])
            x_0 = points[i[0]][0]
            y_0 = points[i[0]][1]
            points.pop(i[0])
            points_new_dim={}
            polar_coor = {}
            break
        while len(points.keys()) > 0:
            for i in list(points.keys()):
                x = (points[i][0] - x_0)*math.cos(angle_x) + (points[i][1] - y_0)*math.sin(angle_y)
                y = -(points[i][0] - x_0)*math.sin(angle_x) + (points[i][1] - y_0)*math.cos(angle_y)
                points_new_dim[i] = (x,y)

            for i in list(points_new_dim.keys()):
                r = math.sqrt(points_new_dim[i][0]**2+points_new_dim[i][1]**2)
                if points_new_dim[i][0] > 0:
                    theta = math.atan(points_new_dim[i][1]/points_new_dim[i][0])
                elif points_new_dim[i][0] < 0 and points_new_dim[i][1] >= 0:
                    theta = math.atan((points_new_dim[i][1]/points_new_dim[i][0])+math.pi)
                elif points_new_dim[i][0] < 0 and points_new_dim[i][1] <0:
                    theta = math.atan((points_new_dim[i][1]/points_new_dim[i][0])-math.pi)

                elif points_new_dim[i][0] == 0 and points_new_dim[i][1] >0:
                    theta = math.pi/2
                elif points_new_dim[i][0] == 0 and points_new_dim[i][1] <0:
                    theta = -math.pi/2
    
                polar_coor[i] = (r,theta)
        
            list_sorted_polar=[]      
            for i in list(polar_coor.keys()):
                list_sorted_polar.append([i,polar_coor[i]])

            list_sorted_polar = sorted(list_sorted_polar, key = lambda list_sorted_polar :(list_sorted_polar[1][1],list_sorted_polar[1][0]))
            first_theta = list_sorted_polar[0][1][1]
            last_theta = list_sorted_polar[-1][1][1]

            new_list_sorted_polar = []
            for i,k in enumerate(list_sorted_polar):
                if k == list_sorted_polar[-1]:
                    new_list_sorted_polar.append(k)
                    break
                j = i
                current_theta = k[1][1]
                next_theta = list_sorted_polar[j][1][1]
                while math.isclose(current_theta,next_theta,abs_tol=0.26) == True:
                    current_theta = next_theta
                    j = j +1
                    if j > len(list_sorted_polar)-1:
                        break
                    else:
                        next_theta = list_sorted_polar[j][1][1]
                new_list_sorted_polar.append([k[0], (k[1][0],current_theta)])

            new_list_sorted_polar = sorted(new_list_sorted_polar, key = lambda new_list_sorted_polar :(new_list_sorted_polar[1][1],new_list_sorted_polar[1][0]))
            for i in new_list_sorted_polar:
                cluster_sequence.append(i[0])
                x_0 = points[i[0]][0]
                y_0 = points[i[0]][1]
                points.pop(i[0])
                points_new_dim={}
                polar_coor = {}
                break
            
        cluster_order["{}".format(main)] = cluster_sequence
    return cluster_order

def write_best_sequences(df_clusters,cluster_sequences,best_aa_per_cluster):
# COMPUTE THE BEST SEQUENCES AND CREATE STR PEPTIDE SEQUENCE

    aa_library = {"ala":"A","asn":"N","gln":"Q","ile":"I","lys":"K","phe":"F","thr":"T","arg":"R","cys":"C","gly":"G","leu":"L","met":"M","ser":"S","pro":"P","asp":"D","glu":"E","his":"H","val":"V","tyr":"Y","trp":"W"}
    nb_big_group = df_clusters['cluster'].max()
    for i in range(1, nb_big_group+1):
        df_clusters_group = df_clusters.loc[df_clusters["cluster"]==i]
        if len(df_clusters_group["cluster"]) <= 3:
            continue
        with open('best_sequences{}.txt'.format(str(i)), 'w') as best_sequences:
            cluster_sequence = cluster_sequences[str(i)]
            index_cluster = 0
            peptide_sequence = ''
            while index_cluster <=len(cluster_sequence)-1:
                aa_3_letters = best_aa_per_cluster[cluster_sequence[index_cluster]][0]
                peptide_sequence = peptide_sequence + aa_library[aa_3_letters]
                x1 = float(df_clusters_group["x_ca"].loc[cluster_sequence[index_cluster]])
                y1 = float(df_clusters_group["y_ca"].loc[cluster_sequence[index_cluster]])
                z1 = float(df_clusters_group["z_ca"].loc[cluster_sequence[index_cluster]])
                if index_cluster == len(cluster_sequence)-1:
                    x2 = float(df_clusters_group["x_ca"].loc[cluster_sequence[0]])
                    y2 = float(df_clusters_group["y_ca"].loc[cluster_sequence[0]])
                    z2 = float(df_clusters_group["z_ca"].loc[cluster_sequence[0]])
                else:
                    x2 = float(df_clusters_group["x_ca"].loc[cluster_sequence[index_cluster+1]])
                    y2 = float(df_clusters_group["y_ca"].loc[cluster_sequence[index_cluster+1]])
                    z2 = float(df_clusters_group["z_ca"].loc[cluster_sequence[index_cluster+1]])

                distance = compute_distance_cluster(x1,y1,z1,x2,y2,z2)
                if distance >=4.3:
                     peptide_sequence = peptide_sequence + math.ceil(distance/4.3 -1)*'g'

                index_cluster = index_cluster+1
            best_sequences.write(peptide_sequence + "\n")
        with open('peptide_centroid.txt','a') as best_positions:
            mean_centroid  = (df_clusters_group["x_ca"].mean() + df_clusters_group["y_ca"].mean() + df_clusters_group["z_ca"].mean())/3
            best_positions.write(str(mean_centroid))
            best_positions.write("\n")
        print(peptide_sequence)

