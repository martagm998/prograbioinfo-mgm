#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Marta Garcia Mondejar
'''

'''
 SCRIPT INFO:
Find domains in Prosite data base between the hits of the blastp
You will have one result for each query.
The results document *_domains.txt.
Files contain: domain name, accession, description and pattern.
'''


import re
from Bio import SeqIO
from Bio.ExPASy import Prosite


def dominios():

    '''
    FUNCTION PARSER THE prosite.dat
    Obtain the patterns
    '''
    h_dom = open("prosite.dat","r")
    result_dom = open("prosite_parsed.tsv","w")
    result_dom.write("NAME"+"\t"+"ACCESSION"+
                     "\t"+"DESCRIPTION"+"\t"+"PATTERN"+"\n")
    records = Prosite.parse(h_dom)

    for record in records:
        n_dom = record.name
        acc_dom = record.accession
        descrip_dom = record.description
        p_dom = record.pattern
        result_dom.write(f"{n_dom}\t{acc_dom}\t{descrip_dom}\t{p_dom}\n")
    records.close()
    result_dom.close()


def prosite_re(result_dom = "prosite_parsed.tsv"):

    '''Function to replace some character of the pattern to re lenguage    '''


    patronfinal = result_dom.replace(".","")
    patronfinal = patronfinal.replace("x",".")
    patronfinal = patronfinal.replace("-","")
    patronfinal = patronfinal.replace("{","[^")
    patronfinal = patronfinal.replace("}","]")
    patronfinal = patronfinal.replace("(","{")
    patronfinal = patronfinal.replace(")","}")
    patronfinal = patronfinal.replace("<","^")
    patronfinal = patronfinal.replace(">","$")

    # Return pattern understandable for re

    return(patronfinal)


def buscardominios(archivo_quey, result_dom="prosite_parsed.tsv"):

    '''Function to search for the domains in the prosite.dat that
    are similar to the ones of our hits of the blastp'''

    for record in SeqIO.parse(archivo_quey, "fasta"):
        entradaf = open(record.id+"_primermuscle.fasta", "r")
        entradaf = entradaf.read()
        salidaf = open("prosite_parsed.tsv","r")
        salidaf = salidaf.read()

        with open(record.id+"_dominios.txt", "w") as lastfile:
            parte1 = entradaf.split("\n")
            parte2 = salidaf.split("\n")

            # Parse the file two lines by two lines

            for i in range(len(parte1) //2):
                secuenciasf = parte1[2*i+1]
                IDsubects = parte1[2*i]
                s_id = str("PATTERNS IN SEQ_"+IDsubects[1:]+":\n\n")
                lastfile.write(s_id)

                # Parse the file with the patterns ignoring the header

                for pattern in parte2[1:]:
                    if pattern != '':

                        # Separate each line of the patterns file

                        pattern = pattern.split("\t")

                        # Transform the pattern into re with FUNCTION: prosite_re

                        patron_re = prosite_re(pattern[3])

                        # If you find the pattern in the sequence write it on the final file

                        if pattern[3] != "":
                            if re.search(patron_re,secuenciasf):
                                lastfile.write(f"\tNAME:{pattern[0]}\n\tACCESSION:{pattern[1]}"
                                               f"\n\tDESCRIPTION:{pattern[2]}\n\tPATTERN:{pattern[3]}\n\n")
                    else:
                        pass
            lastfile.close()
    print("The domains were found successfully")





