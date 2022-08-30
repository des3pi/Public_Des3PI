#!/usr/bin/env python

import os
import argparse

file_path = os.path.realpath(__file__)
des3pi_path = os.path.dirname(file_path)
start_path = os.getcwd()
os.chdir(start_path)



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


args = parser.parse_args()

target = args.i
nb_dock = args.nb
xsbox = args.xs
ysbox = args.ys
zsbox = args.zs
xcbox = args.xc
ycbox = args.yc
zcbox = args.zc

os.system("mkdir des3pi_calculations")
os.chdir(start_path + "/des3pi_calculations")

for nb_dir in range(1,21):
    if nb_dir <10:
        name_dir = "0{}".format(str(nb_dir))
        os.system("mkdir {}".format(str(name_dir)))
    else:
        name_dir = "{}".format(str(nb_dir))
        os.system("mkdir {}".format(str(name_dir)))
    os.system("cp ../{} {}".format(target,str(name_dir)))
    os.chdir(start_path + "/des3pi_calculations/{}".format(str(name_dir)))
    os.system("python3 {}/des3pi_one_run.py -i {} -nb {} -xc {} -yc {} -zc {} -xs {} -ys {} -zs {}".format(des3pi_path,target,str(nb_dock),str(xcbox),str(ycbox),str(zcbox),str(xsbox),str(ysbox),str(zsbox)))
    os.chdir(start_path + "/des3pi_calculations")

os.chdir(start_path + "/des3pi_calculations")
os.system("python3 {}/des3pi_analysis.py".format(des3pi_path))

