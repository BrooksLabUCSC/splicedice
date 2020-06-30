#!/usr/bin/env python3
"""
PATH=/private/groups/brookslab/bin/bedtools-2.25.0/bin/:$PATH

python3 createIRTables_sigmaPSI.py -cov coverage_manifest.txt -m DMSO_v_SSA100_allPSI.tsv -c me_all_clusters2.tsv -o mesa_IR_sigmaPSI.tsv -p percentiles.txt
"""

import time
start_time = time.time()
import pybedtools
import sys
import argparse
import os
import csv
from numpy import median


class CommandLine() :
    """
    Takes command line options
    """
    def __init__(self, inOpts=None) :
            """
            CommandLine constructor.
            Implements a parser to interpret the command line argv string using argparse.
            """
            self.parser = argparse.ArgumentParser(description = 'runRTest.py arguments',
                                                 epilog = 'For more help contact: ',
                                                 add_help = True, #default is True
                                                 prefix_chars = '-',
                                                 usage = '%(prog)s [options] -option1[default]'
                                                 )
             # Add args
            self.parser.add_argument('-cov', '--coverageManifest', action = 'store', required=True, help='manifest file containing paths for *sigmaPSI_percentiles.txt files from sigma_psi.py')
            self.parser.add_argument('-m', '--mesaTable', action = 'store', required=True, help='mesa allPSI.tsv output from quantMESA.py')
            self.parser.add_argument('-c', '--clusters', action = 'store', required=True, help='me_clusters.tsv output from constructMESA.py')
            self.parser.add_argument('-o', '--output', action = 'store', required=True, help='output PSI table with % IR')
            self.parser.add_argument('-p', '--percentiles', action = 'store', required=True, help='percentiles.txt file from sigma_psi.py')

            if inOpts is None :
                self.args = self.parser.parse_args()
            else :
                self.args = self.parser.parse_args(inOpts)

class createIRTable():

    def __init__(self, myCommandLine):
        """
        Instantiate object attributes
        """
        self.coverage = myCommandLine.args.coverageManifest
        self.mesa = myCommandLine.args.mesaTable
        self.output = myCommandLine.args.output
        self.cluster = myCommandLine.args.clusters
        self.percentile = myCommandLine.args.percentiles

    def getCoverageFiles(self):
        """
        Get paths for coverage files
        """

        print('getting paths for coverage files')

        self.abundance = {}

        with open(self.coverage, 'r') as abundance_files:
            for line in abundance_files:
                line = line.strip().split('\t')
                sample = line[0]
                path = line[1]
                self.abundance[sample] = path

        self.header_list = [key+'%IR' for key in self.abundance.keys()]


    def createCoverageDicts(self):
        """
        abundance ={sample1:{chr1:10-100: 99, ...}, sample2:{...}}
        """

        print('creating coverage dictionaries')

        self.abundance_dict = {key:{} for key in self.abundance.keys()}

        for key in self.abundance:
            new_dict = {}
            file_path = self.abundance[key]
            with open(file_path, 'r') as abundance:
                abundance_csv = csv.reader(abundance, dialect='excel-tab')
                header = next(abundance_csv, None)
                for row in abundance_csv:
                    junc_key = row[0] + ':' + row[2]
                    new_dict[junc_key] = float(row[-1])

            self.abundance_dict[key] = new_dict

    def getClusters(self):
        """
        get clusters from mesa cluster file
        """
        print('getting clusters')
        self.clusters = {}

        with open(self.cluster, 'r') as cluster_file:
            for line in cluster_file:
                line = line.strip().split(' ')
                #print (line)
                key = line[0]
                mxes = line[1].split(',')
                self.clusters[key] = mxes

        #print(self.clusters['chr1:24294213-24298339'])

    def getPercentiles(self):
        """
        get percentile file created from sigma_psi.py (?), might be able to remove
        """
        self.percentiles = {}

        with open(self.percentile, 'r') as percentile_file:
            percentile_csv = csv.reader(percentile_file, dialect='excel-tab')
            for row in percentile_csv:
                key = row[0]
                percentiles = row[1:-1]
                self.percentiles[key] = percentiles

        #print(self.percentiles['chr1:24294213-24298339'])

    def getIRTable(self):
        print('creating IR table')
        with open(self.mesa, 'r') as mesa_file:
            mesa_csv = csv.reader(mesa_file, dialect='excel-tab')
            header = next(mesa_csv, None)
            new_header = header[:-1] + self.header_list + [header[-1]]
            with open(self.output, 'w') as out_file:
                print('\t'.join(new_header), file=out_file)
                for row in mesa_csv:
                    #line = line.strip().split('\t')
                    columns = {key:0.0 for key in self.header_list}
                    new_row = self.calculateIR(row, columns)
                    print('\t'.join(new_row), file=out_file)


    def calculateIR(self, row, columns):
        cluster = row[-1]
        mxes = self.clusters[cluster]

        for key in columns:
            key_name = key[:-3]
            try:
                abundance = self.getSigmaPSI(key_name, cluster)
                mxe_abundance = 0
                for mxe in mxes:
                    mxe_a = self.getSigmaPSI(key_name, mxe)
                    mxe_abundance += mxe_a

                IR_val = abundance / (abundance + mxe_abundance)

            except ZeroDivisionError:
                IR_val = 'nan'

            columns[key] = IR_val

        new_row = row[:-1] + [str(columns[key]) for key in self.header_list] + [row[-1]]

        return new_row

    def getSigmaPSI(self, sample, junc):
        sigma_val = []
        chr = junc.split(':')[0]
        for percentile in self.percentiles[junc]:
            try:
                percentile_key = chr + ':' + percentile
                sigma_val.append(self.abundance_dict[sample][percentile_key])
            except KeyError:
                pass

        #print(sigma_val)
        #print(median(sigma_val))
        return median(sigma_val)

    def caseSensitivity(self):
        """
        handling case case sensitivity issues
        """
        print('handling case sensitivity')

        with open('tmp.fasta.fa', 'r') as read_file:
            with open(self.output, 'w') as out_file:
                for line in read_file:
                    if '>' in line:
                        #skip headers
                        pass
                    else:
                        line = line.upper()

                    print(line.strip(), file=out_file)

        os.system('rm tmp.fasta.fa')

def main(myCommandLine=None):
    myCommandLine = CommandLine(myCommandLine)
    abundance = createIRTable(myCommandLine)
    abundance.getCoverageFiles()
    abundance.createCoverageDicts()
    abundance.getClusters()
    abundance.getPercentiles()
    abundance.getIRTable()

    print("Your runtime was %s seconds." % (time.time() - start_time))

main()
