#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Marta Garcia Mondejar
'''


import os
import shutil
import sys

from Bio import SeqIO

import parseandofiles as ps
import maketree
import domains


def HELP():
    print("""\n\033[1;37m\033[4;37mWELCOME TO THE HELP FUNCTION\033[0m\n
          With this program we will make:
          
          - A Blastp from parsed genbanks and a file in fasta format introduced
           separately. This files must be in the same directory as the script.             
          Then the result of that process will be filtered with specific
          Identity and Coverage values that you can provide. If not, the 
          predetermined values will be 30 for identity and 50 for coverage.
          As a result you will have one filtered blastp file for each query
          and will be named with the name of the query_filtrado.tsv
          
          - A filogenetic tree (Neighbor-Joining) by doing an alignment
           with MUSCLE. We will have one tree for each query, saved as 
           queryname_tree.nw files. 
          
          - Finding the protein domains in Prosite database. This file must be
          in the same folder as the scripts. We will have one file for each 
          protein in a proteinname_domain.txt file.
          This file contains domains, accession, description of pattern found.
          
          \033[1;37m\033[4;37mHow to execute the program:\033[0m\n """)

    print("\033[1;37m" + "          ./mainmarta.py" + "\033[1;32m" +
          " [genbank] " + "\033[1;33m" + " [query] " + "\033[1;34m" +
          " [Identity] " + "\033[1;35m" + " [Coverage] " + "\033[0m")

    print(""" 
          \033[1;32m[genbank]\033[0m: The name of the folder with the genbanks.
          
          \033[1;33m[query]\033[0m: Query file in fasta format.
          
          \033[1;34m[Identity]\033[0m: Identity value between 1-100. Optional.
          
          \033[1;35m[Coverage]\033[0m: Coverage value between 1-100. Optional.
          
          All the files will be saved in _Results folder.
          
          \033[1;37m\033[4;37mTHANK YOU FOR USING THIS PROGRAM\033[0m""")



    sys.exit()


# CONTROL OF ARGUMENTS AND PRINT THE HELP

if len(sys.argv) < 3:
    print("\033[1;31mERROR: you are not introducing all "
          "the required files.\033[0m\n"
          "Do you want to see the help of this program? (Yes/No):")
    ayuda = input()
    if ayuda.lower() == "yes":
        HELP()

    elif ayuda.lower() == "no":
        sys.exit()
else:
    pass


# CREATE RESULTS FOLDER WITH PERSONALIZED NAME


print("\nName the folder _Results:\n")
usr_ip=input()
try:
    newfolder = str(usr_ip)+"_Results"
    os.mkdir(newfolder)
except OSError:
    print("\033[1;31mERROR: that folder already exists.\033[0m ")
    sys.exit()


# ARGUMENTS


folder = sys.argv[1]
archivo_query = sys.argv[2]


# FUNCTION PARSER, BLASTP AND FILTER THE BLASTP


ps.Parser(folder)
ps.BlastP(archivo_query)
ps.FILTRAR(archivo_query,'ID_usr', 'COV_usr')


# FUNCTION TO DO ALIGMENT AND THE FILOGENETIC TREE


maketree.NewFasta(archivo_query)
maketree.ALIGN(archivo_query)
maketree.Maketree(archivo_query)


# FUNCTION FIND DOMAINS


domains.dominios()
domains.buscardominios(archivo_query)


#COPY RESULTS FILES TO RESULTS FOLDER



for record in SeqIO.parse(archivo_query, "fasta"):
    files = [record.id + "_results.tsv", record.id + "_filtrado.tsv",
             record.id + "_primermuscle.fasta", record.id + "_align.fasta",
             record.id + "_tree.nw", record.id + "_dominios.txt"]
    for f in files:
        shutil.move(f, str(usr_ip) + "_Results" + "/")

shutil.move("prosite_parsed.tsv", str(usr_ip) + "_Results" + "/")
shutil.move('parserfile.fasta', str(usr_ip)+"_Results"+"/")

print("The files have been copied successfully to _Results folder.\n"
      "Thank you for using this program.")

