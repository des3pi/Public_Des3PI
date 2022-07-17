#!/usr/bin/env python

import os
import csv

def compute_ca(start_dir,nb_dock):
    os.chdir(start_dir)
    aa_pdbqt_list = os.listdir(os.path.dirname(__file__) + "/aa_library")
    print("Found {} elements in library".format(len(aa_pdbqt_list)))

    os.chdir(start_dir + "/docking_results")

    list_ca_fragments = [["fragment_docking_state","aa","score","x_ca","y_ca","z_ca"]]

    for fragment in aa_pdbqt_list:
        fragment_name = fragment.split(".")[0]
        os.chdir(fragment_name)
        for docking in range(1,nb_dock+1):
            scores = []
            x_ca = []
            y_ca = []
            z_ca = []
            os.chdir("docking_{}".format(docking))
            fichier = open("{}.pdbqt".format(fragment_name),"r")
            liste_ligne = fichier.readlines()
            for ligne in liste_ligne:
                if "CA " in ligne:
                    x_ca.append(float(ligne[31:39]))
                    y_ca.append(float(ligne[39:47]))
                    z_ca.append(float(ligne[47:55]))
                if "VINA RESULT" in ligne:
                    scores.append(ligne[25]+ligne[26]+ligne[27]+ligne[28])
            fichier.close()
            for state in range(1,len(x_ca)+1):
                list_temporaire_aa_ca = [fragment_name + "_" + str(docking) + "_" + str(state), fragment_name]
                list_temporaire_aa_ca.append(scores[state-1])
                list_temporaire_aa_ca.append(x_ca[state-1])
                list_temporaire_aa_ca.append(y_ca[state-1])
                list_temporaire_aa_ca.append(z_ca[state-1])
                list_ca_fragments.append(list_temporaire_aa_ca)

            os.chdir(start_dir + "/docking_results/" + fragment_name)
        os.chdir(start_dir + "/docking_results")


    with open('fragment_ca.csv', 'w') as comfile:
        writer = csv.writer(comfile)
        [writer.writerow(r) for r in list_ca_fragments]



