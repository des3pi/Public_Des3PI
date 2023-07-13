# Des3PI

Des3PI proposes peptides sequences targeting a protein-protein interface.

To run Des3PI you need a PDBQT file of your targeted protein, the coordinates of the center of VINA docking box surrounding the interface (xc,yc,zc) and the size of each side of the box (xs,ys,zs).

des3pi.py is the main file to run the program. Open a terminal on the directory of your PDBQT file and execute des3pi.py:

python3 /YOUR/PATH/TO/DES3PI/des3pi.py -i your_pdbqt.pdbqt -nb 50 -xc xc -yc yc -zc zc -xs xs -ys ys -zs zs

The output are generated by class of peptides.

If you proceed to ALL the docking steps previously, you can run this command line to only analysing the data without redoing all the docking steps:

python3 /YOUR/PATH/TO/DES3PI/des3pi.py -i your_pdbqt.pdbqt -nb 50 -xc xc -yc yc -zc zc -xs xs -ys ys -zs zs -analysis_only True
