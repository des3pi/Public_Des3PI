#!/usr/bin/env python

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import seaborn as sns
import math

# FUNCTIONS

def compute_clusters(poses_df,ca_df, csv_output, method='centroid', metric='euclidean', criterion_value=3.5):

    # Compute clusters of the position of the CAs.
    # HCA method : centroid based on euclidean distance
    # Minimal distance of centroids clusters: 3.5 Angstrom

    clusters = []
    Z = linkage(ca_df,method= method,metric= metric)
    plt.title("Classification")
    dendrogram(Z,labels=ca_df.index,orientation='right',color_threshold=0)
    plt.savefig("dendrogramm.pdf")
    clusters_ca = fcluster(Z,t=criterion_value,criterion='distance')
    clusters = clusters + list(clusters_ca)
    index_clusters = np.argsort(clusters_ca)
    #Adding clusters number to every pose
    ca_df['cluster'] = clusters
    poses_df['cluster'] = clusters

    # Create a csv file with the computed clusters
    poses_df.to_csv(csv_output)

    return poses_df

def select_best_clusters(poses_df,ratio=0.001):

    # Compute the number of poses for each selected clusters.
    # Clusters are selected by taking only the best 20 of them AND removing clusters that owns less than 0.1% of the poses.

    clust_percent_list = []
    clust_list = []

    maxi = poses_df['cluster'].max()
    nb_poses = len(list(poses_df['cluster']))
    for j in range(1,maxi+1):
        clust_list.append("{}".format("cluster {}".format(j)))
        df_one_cluster = poses_df.loc[poses_df['cluster']==j,:]
        clust_percent = len(df_one_cluster['cluster'])
        clust_percent_list.append(clust_percent)
    list_combined_percent_clustid = []
    for i in range(0,len(clust_percent_list)):
        list_combined_percent_clustid.append([clust_list[i],clust_percent_list[i]])
    
    list_combined_sorted = sorted(list_combined_percent_clustid, key = lambda pose: pose[1], reverse=True)
    while len(list_combined_sorted)>40000000:
        del list_combined_sorted[-1]
    list_to_del = []
    with open("percent_clust.txt", "w") as perc:
        for i, cluster in enumerate(list_combined_sorted):
            perc.write(str(cluster[1]/nb_poses) + "\n")
    for i, cluster in enumerate(list_combined_sorted):
     
           
        if cluster[1]/nb_poses <= ratio:
            list_to_del.append(i)
    list_to_del = sorted(list_to_del, reverse = True)

    for i in list_to_del:
        del list_combined_sorted[i]

    clust_percent_list = []
    clust_list = []
    for cluster in list_combined_sorted:
        clust_percent_list.append(cluster[1])
        clust_list.append(cluster[0])

    #print(clust_percent_list)
    #print(clust_list)
    results_df = pd.DataFrame()
    results_df['cluster'] = clust_list
    results_df['nb_poses'] = clust_percent_list
    # Plot the number of poses for every selected clusters
    sns.barplot(x="cluster", y="nb_poses", data=results_df)
    plt.savefig("poses_per_clusters.pdf")

    return results_df


def compute_clusters_positions(poses_df,results_df):

    results_list = results_df['cluster'].tolist()

    # Write several .xyz files that represents the coordinates of the mean CA for each selected clusters
    ca_moy_list  = []
    for i,cluster in enumerate(results_list):
        with open('cluster{}.xyz'.format(cluster.split(" ")[1]), 'w') as xyz:
            #print(i,cluster)
            xyz.write('1\n')
            xyz.write('{}\n'.format(cluster))
            aa_df = poses_df.loc[(poses_df['cluster'] == float(cluster.split(" ")[1])),:]
            x_ca,y_ca,z_ca = round(aa_df['x_ca'].mean(),3), round(aa_df['y_ca'].mean(),3), round(aa_df['z_ca'].mean(),3)
            tuple_ca = (x_ca,y_ca,z_ca)
            ca_moy_list.append(tuple_ca)
            xyz.write('C {} {} {}\n'.format(x_ca,y_ca,z_ca))
    # Create a dataframe with the mean coordinates for each clusters
    listxca = []
    listyca = []
    listzca = []
    for ca_moy in ca_moy_list:
        listxca.append(ca_moy[0])
        listyca.append(ca_moy[1])
        listzca.append(ca_moy[2])

    df_clusters = pd.DataFrame()
    df_clusters['cluster'] = results_list
    df_clusters['x_ca'] = listxca
    df_clusters['y_ca'] = listyca
    df_clusters['z_ca'] = listzca
    df_clusters = df_clusters.set_index('cluster')

    return df_clusters

