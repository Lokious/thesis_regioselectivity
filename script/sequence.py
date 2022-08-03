#!/usr/bin/env python3
"""
Author:         Yingjie Shao
Description:
Dependencies:   Python3.9
                numpy
                pandas
This code is to parse sequence data
"""
import sys
from Bio import SeqIO
class sequences():

    def __init__(self, msa_model=None):
        self.model = msa_model

    def get_acative_site(self,seq = "",num_AA:int=0):
        """
        This function is for return the sepreate aa sequences around acative site
        :param seq:
        :return:
        """
        print("unfinihed function")


    def remove_duplicate(self,files):
        for file in files:
            with open('data/{}_rm.fasta'.format(file), 'a') as outFile:
                record_ids = list()
                for record in SeqIO.parse("data/{}.fasta".format(file), 'fasta'):
                    if record.id not in record_ids:
                        record_ids.append(record.id)
                        SeqIO.write(record, outFile, 'fasta')
        else:
            print("finished removing duplicate")
