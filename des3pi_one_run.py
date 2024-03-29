#!/usr/bin/env python

import sys
sys.path.insert(0, "/home/delaunay/des3pi_vtest")
import os 
import argparse
import pandas as pd
import numpy as np
import run_simulation
import compute_ca
import clustering
import find_best_sequences


# Arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input pdbqt file of the targeted protein",
                    type=str, required=True)
parser.add_argument("-nb", help="Number of docking per fragment",
                    type=int, default=100)
parser.add_argument("-xs", help="x size of the box",
                    type=float, default=0)
parser.add_argument("-ys", help="y size of the box",
                    type=float, default=0)
parser.add_argument("-zs", help="z size of the box",
                    type=float, default=0)
parser.add_argument("-xc", help="x center of the box",
                    type=float, default=0)
parser.add_argument("-yc", help="y center of the box",
                    type=float, default=0)
parser.add_argument("-zc", help="z center of the box",
                    type=float, default=0)
parser.add_argument("-analysis_only", help="analysis without docking",
                    type=str, default="False")

args = parser.parse_args()

target = args.i
nb_dock = args.nb
xsbox = args.xs
ysbox = args.ys
zsbox = args.zs
xcbox = args.xc
ycbox = args.yc
zcbox = args.zc
analysis_only_bool = str(args.analysis_only)
print("The selected target is "+target)

start_dir = os.getcwd()
os.chdir(start_dir)

os.system("mkdir docking_results")

### Compute the parameters in case of blind docking
if xsbox == 0 and ysbox == 0 and zsbox == 0 and xcbox ==0 and ycbox ==0 and zcbox ==0:
    xcbox,ycbox,zcbox,xsbox,ysbox,zsbox = run_simulation.compute_blind_docking_parameters(start_dir, target)

### Run the simulation
os.chdir(start_dir)
if analysis_only_bool == "False":
    run_simulation.run_docking(start_dir,target,nb_dock,xcbox,ycbox,zcbox,xsbox,ysbox,zsbox)

### Convert pdbqt to pdb for the next analysis
os.chdir(start_dir)
if analysis_only_bool == "False":
    run_simulation.pdbqt_to_pdb(start_dir,nb_dock)

### Compute CAs 

compute_ca.compute_ca(start_dir,nb_dock)

### Run Clustering

#Create a dataframe with only the position of the fragments CAs
os.chdir(start_dir + "/docking_results")
poses_df = pd.read_csv("fragment_ca.csv")
poses_df = poses_df.set_index("fragment_docking_state")
ca_df = poses_df.drop(['score'], axis=1)
ca_df = ca_df.drop(['aa'], axis = 1)

# Add cluster id to each pose
poses_df = clustering.compute_clusters(poses_df,ca_df, csv_output='fragment_ca_clusters.csv')

# Create a df with only the best clusters
results_df = clustering.select_best_clusters(poses_df)

# Create a df with every centroid position for each cluster of interest
df_clusters = clustering.compute_clusters_positions(poses_df,results_df)

# Re cluster to have the main possible peptides 

df_clusters = clustering.compute_clusters(df_clusters,df_clusters,csv_output='main_clusters.csv',criterion_value=12.6)

# Compute best aa for each cluster

best_aa_per_cluster, best_aa_nb_poses = find_best_sequences.compute_more_represented_aa(poses_df,df_clusters,nb_dock)
with open("best_aa_per_cluster.txt","w") as best:
    for i in best_aa_per_cluster.keys():
        chaine = i + ": " + best_aa_per_cluster[i][0] + "\n"
        best.write(chaine)
# Find optimal cluster order for each big group
print(df_clusters)
cluster_sequences = find_best_sequences.get_optimal_order_cluster(results_df,df_clusters)
print(cluster_sequences)
if cluster_sequences == []:
     pass
else:
    with open("cluster_order.txt","w") as orders:
        for key in list(cluster_sequences.keys()):
            toadd = " ".join(cluster_sequences[key])
            orders.write(toadd)
            orders.write("+")
# Write and display possible peptides

find_best_sequences.write_best_sequences(df_clusters,cluster_sequences,best_aa_per_cluster)

# Build the peptides

os.chdir(start_dir + "/docking_results")

#with open('best_sequences1.txt','r') as best_sequences:
#    list_sequences = best_sequences.readlines()
#    for sequence in list_sequences:
#        build_peptide.build_peptide(sequence.split("\n")[0])

os.chdir(start_dir)

