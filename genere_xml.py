"""
Ecrit un  arbre au format XML .
Usage:
  genere_xml.py <treeFile> <alignmentFile> <resultsFile>

Positional arguments:
  treeFile              family name / phylogenic trees in newick format
  alignmentFile         multiple alignment in fasta format
  resultsFile           statistics

"""

#
# Modified by Genna BEN HASSEN, Camille SIHARATH and Grégoire ALIZADEH NOIRET, 16 December 2021
# Modified by Grégoire ALIZADEH NOIRET, 19 July 2022
# Modified by Laurent GUÉGUEN, 6 March 2023
#

import argparse
import json
import sys
import re
import random
import os
from Bio import Phylo
from io import StringIO
from math import *

from lxml import etree
from xml.etree import ElementTree
from xml.dom import minidom

from lxml.etree import XMLParser

parser = argparse.ArgumentParser(description='Generate an XML file from \
    a tree, an alignment and statistics.')
## Files
parser.add_argument('-t', '--tree', dest='treeFile', action='store',\
    required=True,\
    help='File containing a tree in Newick format')
parser.add_argument('-a', '--align', dest='alignmentFile', action='store',\
    required=True,\
    help='File containing an alignment in FASTA format')
parser.add_argument('-r', '--results', dest='resultsFile', action='store',\
    required=True,\
    help='File containing the results')
parser.add_argument('-o', '--output', dest='output', action='store',\
    help='Name of the output XML file (if not specified, the XML will have \
    the same name as the tree file)')
## Other parameters
# parser.add_argument('-c', '--col', dest='statcol', action='store', type=int,\
#     default=1,\
#     help='Index of the results column to use')
parser.add_argument('-e', '--exp', dest='exprcol', action='store', type=str,\
    default=1,\
    help='Formula to be applied on the column values')
parser.add_argument('-n', '--nostat', dest='nostat', action='store', type=float,\
    default=-1.0,\
    help='Value to use in case there is no statistic associated\
    with a site in the sequence')
## Mode
parser.add_argument('-b', '--branchsite', dest='isBranchsite', action='store_true',\
    required=False,\
    help='View site-branch data')
## Toggles
parser.add_argument('--skipmissing', dest='skipMissingSites', action='store_true',\
    required=False,\
    help='Prevent the addition of special values (-n, --nostat) for sites that are absent from the results file. \
        This results in sites being next to each other on the graph even though their positions are distant.')
args = parser.parse_args()


def loadDico(fileDico):
    """
    Fonction qui prend un fichier d'alignement en entree et retourne un dictionnaire des noms d'espece des sequences
    """

    file = open(fileDico,"r")
    speciesDico = {}
    for line in file:
        if line[0]==">":
            speciesDico[line.split('>')[1].split('\n')[0]] = line.split('>')[1].split('\n')[0]
            # tline = re.split('_',line)
            # #key = f"{tline[0]} {tline[1]}"
            # key = line.split('\n')[0]
            # fin = tline[4].split("\n")[0]
            # speciesDico[key.split(">")[1]] = f"{tline[2]} {tline[3]} {fin}"
    file.close()
    return speciesDico


def geneticCode():
  matches = {
    # DNA
    'AAA':'K', 'AAC':'N', 'AAG':'K', 'AAT':'N',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AGA':'R', 'AGC':'S', 'AGG':'R', 'AGT':'S',
    'ATA':'I', 'ATC':'I', 'ATG':'M', 'ATT':'I',
    'CAA':'Q', 'CAC':'H', 'CAG':'Q', 'CAT':'H',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'GAA':'E', 'GAC':'D', 'GAG':'E', 'GAT':'D',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'TAA':'*', 'TAC':'Y', 'TAG':'*', 'TAT':'Y',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TGA':'*', 'TGC':'C', 'TGG':'W', 'TGT':'C',
    'TTA':'L', 'TTC':'F', 'TTG':'L', 'TTT':'F',

    # RNA
    'AAU':'N',
    'ACU':'T',
    'AGU':'S',
    'AUA':'I', 'AUC':'I', 'AUG':'M', 'AUU':'I',
    'CAU':'H',
    'CCU':'P',
    'CGU':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'GAU':'D',
    'GCU':'A',
    'GGU':'G',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'UAA':'*', 'UAC':'Y', 'UAG':'*', 'UAU':'Y',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UGA':'*', 'UGC':'C', 'UGG':'W', 'UGU':'C',
    'UUA':'L', 'UUC':'F', 'UUG':'L', 'UUU':'F',

  }
  return(matches)


def nucToAmino(nuc_seq:str):
    '''Converts a DNA/RNA sequence into a proteic sequence.'''

    matches=geneticCode()
    matches["---"] = "-";
    aa = ''
    codons = [nuc_seq[i:i+3].upper() for i in range(0, len(nuc_seq), 3)]
    for codon in codons:
        if len(codon) == 3:
            try:
                aa += matches[codon]
            except KeyError:
                print(f"Unknown codon: {codon}, adding 'X' to sequence")
                aa += 'X'
        else:
            print(f"Ignored: {codon} (not a codon)")
    return aa


def loadAlignment(alignmentFile):
    '''Loads a FASTA alignment.'''

    alignmentDict = {}
    maxSeqIdLength = 0
    with open(alignmentFile, 'r') as af:
        for line in af:
            if line[0] == '>': # sequence id
                seq_id = line.strip('\n')[1:]
                if len(seq_id) > maxSeqIdLength:
                    maxSeqIdLength = len(seq_id)
                alignmentDict[seq_id] = ''
            else: # sequence
                aligned_seq = line.strip('\n')
                if len(alignmentDict[seq_id]) == 0:
                    alignmentDict[seq_id] = aligned_seq
                else:
                    alignmentDict[seq_id] += aligned_seq
    return alignmentDict, maxSeqIdLength


def loadResultsSites(resultsFile, exprcol="'[1]'", nostat_value=-1.0, skipMissingSites=False):
    '''Loads site results from formula exprcol on columns.'''

    if len(exprcol)>30:
      return
    expr=exprcol[1:-1].replace("[","float(line[")
    expr=expr.replace("]","])")
    resultsDict = {}
    with open(resultsFile, 'r') as f:
        i = -1
        for line in f:
            line = line.strip().split()
            if re.search('^[0-9]+$', line[0]): # results line
                try:
                    site = int(line[0])
                    res = float(eval(expr))
                except ValueError:
                    print(f"Conversion failed in line {line}")
                else:
                    if i == -1:
                        i = site
                    while i < site and i < MAX_SITE:
                        # in case there is no statistic for
                        # a site, give it a default value
                        if not skipMissingSites:
                            # print(f'Site {i} missing ({site})')
                            resultsDict[i] = nostat_value
                        i += 1
                    resultsDict[i] = res
            i += 1

    resultsList = []
    for key, item in resultsDict.items():
        resultsList.append(item)
    return resultsList


def loadResultsBranchSite(resultsFile, nostat_value=-1.0, skipMissingSites=False):
    '''Loads branch-site results.'''

    col_lists = []
    with open(resultsFile, 'r') as f:
      col_headers=f.readline().strip().split()[1:]
      nbcol=len(col_headers)
      col_lists=[[] for i in range(nbcol)]

      numl = -1
      for line in f:
        line = line.strip().split()
        site = int(line[0].strip("[]"))

        ## line: ['0', '0.25866598', '0.66856214', '0.19517522', ...]
        if numl==-1:
          numl=site
        if not skipMissingSites:
          while numl < site and numl < MAX_SITE:
          # in case there is no statistic for
          # a site, give it a default value
            for i in range(nbcol):
              col_lists[i].append(nostat_value)
            numl += 1
        else:
          numl=site

        for i in range(nbcol):
          col_lists[i].append(line[i+1])

        numl+=1

    d_cols = {}
    for i in range(len(col_lists)):
      d_cols[col_headers[i]] = col_lists[i]

      ## col_lists: [['0.00885554', '0.25866598', '0.03920189'], ...]
      ## d_cols {'38': ['0.00885554', '0.25866598', '0.03920189'], ...}

    d_cols_2 = dict(d_cols)
    nbs=len(d_cols['0'])
    for col_key in d_cols:
      ## For every branch
      if re.search(r'^[0-9]+$', col_key):
        col_text = ''
        ## Add each result to the text
        for i in range(nbs):
          col_text += f'{d_cols[col_key][i]}, '
          ## col_text: '0.1248, 0.12381, ...'
          ## Add brackets for JSON format
        col_text = '['+col_text[:-2]+']' # [:-2] to remove last space and comma
        d_cols_2[col_key] = col_text
        ## d_cols_2: {1: '[0.1248, 0.12381, ...]', ...}

    nb_branches = len(col_lists)
    print(nb_branches, 'branches found in results')

    return d_cols_2


def getColnames(file):
    '''Returns a list of headers from a file.'''

    with open(file, 'r') as f:
        headers =  f.readline().rstrip().split()
    return headers


def cleanTree(tree:str):
    '''Adds branch numbers to a tree.'''

    tree = re.sub(r'([\),])([0-9]+):', r'\1:', tree)
    tree = re.sub(r'([\),])(:[0-9]+):', r'\1:', tree)

    new_tree = ''
    branch_id = 0
    if ':' in tree:
        for i in range(len(tree)):
            if tree[i] in ':':
                new_tree += ':' + str(branch_id) + tree[i]
                branch_id += 1
            else:
                new_tree += tree[i]
    else:
        for i in range(len(tree)):
            if tree[i] in ',)':
                new_tree += ':' + str(branch_id) + ':' + str(1) + tree[i]
                branch_id += 1
            else:
                new_tree += tree[i]
    return new_tree


def createPhyloXML(fam,newick,results):
    newick = cleanTree(newick)

    # Parse and return exactly one tree from the given file or handle
    if not ':' in newick:
        nv_arbre = ""
        for i in range(len(newick)):
            if (newick[i] in ',)' and newick[i-1] != ')'):
                nv_arbre += ":0.4"
                nv_arbre += newick[i]
            elif (newick[i] in ',)' and newick[i-1] == ')'):
                nv_arbre += ":0.2"
                nv_arbre += newick[i]
            else:
                nv_arbre += newick[i]
        newick = nv_arbre

    handle = StringIO(newick)
    trees = Phylo.read(handle, 'newick')

    ## Branch IDs match the order in which <clade></clade> tags are closed in the PhyloXML tree,
    ## the first ID being 0.

    # Write a sequence of Tree objects to the given file or handle
    rd = str(random.randint(0,10000))
    Phylo.write([trees], 'tmpfile-'+rd+'.xml', 'phyloxml')
    file = open('tmpfile-'+rd+'.xml', 'r')
    # Copies all objects in a variable and removes the created file
    text = file.read()
    file.close()
    os.remove('tmpfile-'+rd+'.xml')

    p = XMLParser(huge_tree=True)
    text = text.replace("phy:", "")

    text = re.sub("b'([^']*)'", "\\1", text)
    text = re.sub('branch_length_attr="[^"]+"', "", text)
    header = "<phyloxml>"

    text = re.sub('<phyloxml[^>]+>', header, text)
    text = text.replace('Phyloxml', 'phyloxml')

    def count_repl(match):
        global current_branch
        current_branch += 1
        return f'<closing_order>{current_branch}</closing_order></clade>'
    text = re.sub(r"</clade>", repl=count_repl, string=text)

    tree = etree.fromstring(text, parser=p)
    treename = etree.Element("name")
    treename.text = fam

    ins = tree.find('phylogeny')
    ins.append(treename)

    clade = tree.xpath("/phyloxml/phylogeny/clade")
    subtree = tree.xpath("/phyloxml")
    nbfeuille = 0
    famspecies = {}
    lenseq = 0

    # Check results

    if type(results)==type([]): # Site results
      lenres=len(results)
    else:
      lenres=[len(eval(v)) for v in results.values()][0]

    for element in clade[0].iter('clade'):
        # look for a <name> element in the current <clade> element
        enom = element.find('name')
        if (enom is not None):
            # if there is a <name> element, it means we're in a leaf
            nbfeuille = nbfeuille + 1
            cds = enom.text
            sp = dico.get(cds)
            if (not  sp):
                print ("undefined species for "+ cds)
                sp = "undefined"

            famspecies[sp] = 1

            ## Find sequence for current leaf name
            seq_alg = alignmentDict.get(cds)
            if not seq_alg:
                print ("undefined alignment for "+ cds)
                seq_alg = ""
            else:
              if lenseq==0:
                lenseq=len(seq_alg)
              elif lenseq!=len(seq_alg):
                print ("sequences of different length for "+ cds)
                seq_alg = ""

            evrec = etree.Element("eventsRec")
            leaf = etree.Element("leaf")
            leaf.set('speciesLocation', sp)
            if 'seqdefdico' in globals():
                if cds in seqdefdico:
                    leaf.set('definition', seqdefdico[cds])

            ## Add sequence to leaf
            if lenres==lenseq/3:
                isCodon="true"
                leaf.set('dnaAlign', seq_alg) # coding sequence
                leaf.set('aaAlign', nucToAmino(seq_alg)) # translated sequence
            elif lenres==lenseq:
                isCodon="false"                
                leaf.set('aaAlign', seq_alg) # raw sequence (any type)
            else:
                sys.exit("Mismatch lengths between alignment " + str(lenseq) + " and sites " + str(lenres))

            evrec.append(leaf)
            element.append(evrec)

    ## Match branch IDs and branch results, then add branch info
    if args.isBranchsite:
        for element in tree.iter('clade'):
            branch_id = element.find('closing_order').text
            try:
                ## Look for the column with header <branch_id> in the results
                branch_info = etree.Element('branch_info')
                branch_info.set('id', branch_id)
                branch_info.set('results', str(dict_results[branch_id]))
                element.append(branch_info)
            except KeyError:
                ## If there is no matching column is results, add one with negative results
                print(f'Branch ID {branch_id} not found, adding placeholder')
                branch_info = etree.Element('branch_info')
                branch_info.set('id', branch_id)
                dummy_col = json.loads(dict_results['0'])
                dummy_col = [-1.0 for _ in range(len(dummy_col))]
                col_str = json.dumps(dummy_col)
                branch_info.set('results', col_str)
                element.append(branch_info)

    print ("Number of leaves : ")
    print (nbfeuille)
    print ("Length of sequences : ")
    print (lenseq)
    print ("Length of results : ")
    print (lenres)

    LengthMaxSeqID = etree.Element('maxSeqIdLength')
    LengthMaxSeqID.text = str(maxSeqIdLength)

    treesize =  etree.Element("size")
    treesize.set('leaves',str(nbfeuille))
    treesize.set('isCodon',str(isCodon))
    e=subtree[0].find('phylogeny')
    e.append(treesize)
    e.append(LengthMaxSeqID)

    ## Add global results to tree
    if not args.isBranchsite:
      globalResultsElement = etree.Element('global_results')
      globalResultsElement.set('results', str(results))
      e.append(globalResultsElement) # add the tag containing results

    geneticcode = etree.Element("geneticCode")
    gC = geneticCode()
    for cod,aa in gC.items():
      geneticcode.set(cod,aa)

    e.append(geneticcode) # add the tag containing geneticcode


    text =  minidom.parseString(ElementTree.tostring(subtree[0])).toprettyxml()
    # remove blank lines
    cleantext = "\n".join([ll.rstrip() for ll in text.splitlines() if ll.strip()])
    # print(cleantext)
    return cleantext

MAX_SITE = 100000
sys.setrecursionlimit(15000)

print ("Loading alignment... ")
loadedAlignment = loadAlignment(args.alignmentFile)
alignmentDict, maxSeqIdLength = loadedAlignment[0], loadedAlignment[1]
print ("OK")

print ("Loading results... ")
if not args.isBranchsite:
  results =  loadResultsSites(args.resultsFile, args.exprcol, args.nostat, skipMissingSites=args.skipMissingSites)
else:
  results = loadResultsBranchSite(args.resultsFile)

if args.isBranchsite:
  dict_results = loadResultsBranchSite(args.resultsFile)
print("OK")

# Creates empty phyloxml document
# project = Phyloxml()  # uncomment to have a unique xml file

# Loads Species name dico
dico = loadDico(args.alignmentFile)

# Loads newick tree
treefile = open(args.treeFile, "r")
if args.output:
    output_name = args.output
else:
    if '.' in args.treeFile:
        output_name = args.treeFile[::-1].split('.', 1)[1][::-1]
    else:
        output_name = args.treeFile
xmloutputfile = open(output_name,"w")

for line in treefile:
    # tline = re.split(' ',line)
    # if len(tline) > 1:
    #     newick=tline[1]
    #     fam=tline[0]
    # else:
    #     newick = tline[0]
    #     fam = ''
    current_branch = -1
    try:
      phyloxmltree = createPhyloXML("",line,results)
      xmloutputfile.write(phyloxmltree)
      print("Tree OK")
    except ValueError:
      sys.exit("Mismatch lengths")

treefile.close()
xmloutputfile.close()
