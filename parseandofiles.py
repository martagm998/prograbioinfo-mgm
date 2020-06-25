#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Marta Garcia Mondejar
'''

'''
SCRIPT INFO:
This script will make BlastP for the sequences of a genbank.
You will also need a file with the proteins query in fasta format.
The final results will be filtered and save in *_filtrado.tsv
'''


import os
import pandas as pd
import sys
from Bio import SeqIO
from subprocess import Popen, PIPE


def Parser(folder):

    '''Function to check the genbank format and
    create a multifasta'''

    directorio = os.getcwd()

    try:
        # If the results file already exists, delete it

        if os.path.isfile("parserfile.fasta") == True:
            os.remove("parserfile.fasta")
        else:
            pass

        output_handle = open("parserfile.fasta", "a")
        os.chdir(folder)

        # Change to work in the folder where the Genbanks are

        os.listdir()

        #Check if the files are genbank

        for file in os.listdir():
            with open(file, "r") as input_handle:
                for seq_record in SeqIO.parse(input_handle, "genbank"):
                    for seq_feature in seq_record.features:
                        try:
                            if seq_feature.type == 'CDS':
                                output_handle.write(f">{seq_feature.qualifiers['locus_tag'][0]}"
                                                    f"@{seq_record.name}\n{seq_feature.qualifiers['translation'][0]}\n")
                        except:
                            pass

            input_handle.close()
        output_handle.close()

        # Get back to main folder

        os.chdir(directorio)

    except:
        print("\033[1;31m"+"Wrong format. Please introduce the correct "
                           "Genbank folder"+"\033[0m")
        sys.exit()


def BlastP(archivo_query, subject = 'parserfile.fasta'):

    '''
    FUNCTION BLASTP: two needed files
    archivo_query (fasta format) incorporated by the user
    subject is the file obtained from the Parser function
    '''


    # Check if the input file is in fasta format

    with open(archivo_query, "r") as archivo:
        for record in SeqIO.parse(archivo, "fasta"):
            seq = str(">%s\n%s" %(record.id, record.seq))
            temp = open(record.id+"temporal.fasta", "w")
            temp.write(seq)
            temp.close()
            with open(record.id+"_results.tsv", "w") as blastp_results:

                # Process of blastp

                blastp = Popen(['blastp', '-query', record.id+"temporal.fasta",
                                '-subject', subject, '-evalue', '0.00001', '-outfmt',
                                 "6 qseqid qcovs pident evalue sseqid sseq "],stdout=PIPE, stderr=PIPE)
                cabecera_blastp = str("queryID\tCover\tIdent\tevalue\tsubjectID\tSUBJECTseq\n")

                # Save the results of the blastp with header

                result_blastp = blastp.stdout.read().decode("utf-8")
                blastp_results.write(cabecera_blastp)
                blastp_results.write(result_blastp)
                blastp_results.close()
                os.remove(record.id + "temporal.fasta")
        print("\nBlastp done successfully")



def FILTRAR(archivo_query, ID_usr, COV_usr, blastp="blastp_result.tsv"):

    '''
    Function to filter blast results with specific identity and coverage
    values. If the inputs values are out of range (1-100) or arent added,
    the predetermined values will be 30 for ID and 50 for coverage
    '''

    # IDENTITY

    try:
        if int(sys.argv[3]) >= 0 and int(sys.argv[3]) <= 100:
            ID_usr = sys.argv[3]
            pass
        else:
            print("\033[1;31m"+"Error: ID value out of range (0-100). Please try again"+"\033[0m")
            sys.exit()
    except:
        ID_usr = 30.0


    # COVERAGE

    try:
        if int(sys.argv[4]) >= 0 and int(sys.argv[4]) <= 100:
            COV_usr = sys.argv[4]
            pass
        else:
            print("\033[1;31m"+"Error: coverage value out of range (1-100). Please try again"+"\033[0m")
            sys.exit()
    except:
        COV_usr = 50.0


    # Filter process

    for record in SeqIO.parse(archivo_query, "fasta"):
        with open(record.id+"_results.tsv") as tsvfile, \
                open(record.id+"_filtrado.tsv","w") as tsv_filtred:
            tsvreader = pd.read_csv(tsvfile, delimiter='\t')
            trying = tsvreader.loc[(tsvreader['Ident']>=int(ID_usr)) &
                                   (tsvreader['Cover']>=int(COV_usr)), :]
            trying.to_csv(tsv_filtred, sep='\t')
            tsvfile.close()
            tsv_filtred.close()

