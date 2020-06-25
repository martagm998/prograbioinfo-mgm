#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Marta Garcia Mondejar
'''

''' 
SCRIPT INFO:
This script will create a filogenetic tree Neighbor-Joining.
We will use MUSCLE and will have one document per query.
The results document for the tree will be saved in *_tree.nw
'''


from Bio import SeqIO
from subprocess import Popen, PIPE


def NewFasta(archivo_query):

    '''
    Function to convert .tsv file from the blastp
    to a fasta file
    '''

    for record in SeqIO.parse(archivo_query,"fasta"):
        test = open(record.id+"_filtrado.tsv","r")
        muscle_out = open(record.id+"_primermuscle.fasta","w")
        muscle_out.write(f">{record.id}\n{record.seq}\n")

        with open(record.id+"_primermuscle.fasta","a+") as muscle_out:
            for row in test.readlines():
                fields = row.rstrip().split("\t")
                if fields[6] != 'SUBJECTseq':
                    subject_ID = fields[5]
                    SUBJECTseq = fields[6]
                    muscle_out.write(f">{subject_ID}\n{SUBJECTseq}\n")
            muscle_out.close()


def ALIGN(archivo_query, file='record.id+"_primermuscle.fasta"'):

    '''
    Function to do the align between the querys and subjects
    '''

    for record in SeqIO.parse(archivo_query, "fasta"):
        with open(record.id+"_primermuscle.fasta","r") as file, \
                open(record.id+"_align.fasta", "a") as align:
            process = Popen(['muscle','-in', record.id+"_primermuscle.fasta",
                             "-out", record.id+"_align.fasta"], stdout=PIPE, stderr=PIPE)
            salida = process.stdout.read().decode("utf-8")
            process.stdout.close()
            align.write(salida)
        align.close()
        file.close()
    print("Alignment done successfully")


def Maketree(archivo_query, input_m='record.id+"_align.fasta"'):

    '''
    Function to create the filogenetic tree, one for each query.
    We use the alignment files from the above function
    '''

    for record in SeqIO.parse(archivo_query, "fasta"):
        with open(record.id+"_align.fasta", "r") as input_m,\
                open(record.id+"_tree.nw","w") as output_tree:
            tree = Popen(['muscle', '-maketree', '-in', record.id+"_align.fasta",
                          '-out',record.id+"_tree.nw",'-cluster', 'neighborjoining'],
                         stderr=PIPE)
        input_m.close()
        output_tree.close()
    print("Tree done successfully")
