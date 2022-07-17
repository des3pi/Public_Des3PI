#!/usr/bin/env python

from Bio import PDB
import os
import sys 

def compute_blind_docking_parameters(start_dir, target):
# Calculate box size and center in case of a blind docking 
    x = []
    y = []
    z = []
    parser = PDB.PDBParser()
    io = PDB.PDBIO()
    struct = parser.get_structure('{}'.format(target.split(".")[0]), target)
    for model in struct:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    x.append(atom.get_coord()[0])
                    y.append(atom.get_coord()[1])
                    z.append(atom.get_coord()[2])

    xsbox = max(x) - min(x) +2
    ysbox = max(y) - min(y) +2
    zsbox = max(z) - min(z) +2
    xcbox = (max(x) + min(x))/2
    ycbox = (max(y) + min(y))/2
    zcbox = (max(z) + min(z))/2
    return xcbox,ycbox,zcbox,xsbox,ysbox,zsbox


def run_docking(start_dir,target,nb_dock,xc,yc,zc,xs,ys,zs):
# Create the config file for each fragment
    os.chdir("docking_results")
    aa_pdbqt_list = os.listdir(os.path.dirname(__file__) + "/aa_library")
    for fragment in aa_pdbqt_list:
        os.system("mkdir {}".format(fragment.split(".")[0]))
        os.chdir("{}".format(fragment.split(".")[0]))
        config_file = open("config.txt", "w")
        config_file.write("receptor = "+ start_dir + "/{}\n".format(target))
        config_file.write("ligand = " + os.path.dirname(__file__) + "/aa_library/{}\n".format(fragment))
        config_file.write("\n")
        config_file.write("center_x = {}\n".format(str(xc)))
        config_file.write("center_y = {}\n".format(str(yc)))
        config_file.write("center_z = {}\n".format(str(zc)))
        config_file.write("\n")
        config_file.write("size_x = {}\n".format(str(xs)))
        config_file.write("size_y = {}\n".format(str(ys)))
        config_file.write("size_z = {}".format(str(zs)))
        config_file.close()
        os.chdir(start_dir)
        os.chdir("docking_results")
# Run docking
    for fragment in aa_pdbqt_list:
        os.chdir(fragment.split(".")[0])
        print("starting with {}".format(fragment.split(".")[0]))
        for dock in range(1,nb_dock+1):
            print("docking {}".format(dock))
            os.system("mkdir docking_{}".format(dock))
            os.chdir("docking_{}".format(dock))
            os.system("vina --config ../config.txt --out ./{}".format(fragment))
            os.chdir(start_dir)
            os.chdir("docking_results")
            os.chdir(fragment.split(".")[0])
        os.chdir(start_dir)
        os.chdir("docking_results")


def pdbqt_to_pdb(start_dir, nb_dock):
# CONVERT PDBQT RESULTS FILE INTO PDB FOR NEXT STEPS

    aa_library=["ala","trp","thr","tyr","phe","ser","val","gly","glu","asp","asn","arg","lys","his","met","gln","cys","leu","ile","pro"]


    for aa in aa_library:
        os.chdir(start_dir + "/docking_results/{}".format(aa))
        for docking in range(1,nb_dock+1):
            os.chdir("docking_{}".format(docking))
            os.system("babel -ipdbqt {}.pdbqt -opdb {}.pdb".format(aa, aa))
            os.chdir(start_dir + "/docking_results/{}".format(aa))
        os.chdir(start_dir)








