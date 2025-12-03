# Methods to perform ancestral sequence reconstruction on a tree, given sequences at leaves
#
# Up to now, only maximum parsimony. Indels are minimized (ie most
# recent as possible).
#
# To be done, indel mgmt

from tree import *
import argparse
import sys

def SubsCost(a,b):
  """Costs for substitutions. Indels cost likewise.

  """

  if a=="*" or b=="*":
    return 1
  if a!=b:
    return 1
  else:
    return 0
  
def add_cost_to(cost1, cost2):
  """Add disctionary cost2 in cost1 (that is changed).
  
  """
  for k,v in cost2.items():
    if k in cost1:
      cost1[k] += v
    else:
      cost1[k] = v

      
def loadAlignment(alignmentFile):
  """Loads a FASTA alignment.

  """
  alignmentDict = {}
  with open(alignmentFile, 'r') as af:
        for line in af:
            if line[0] == '>': # sequence id
                seq_id = line.strip('\n')[1:]
                alignmentDict[seq_id] = ''
            else: # sequence
                aligned_seq = line.strip('\n').upper()
                if len(alignmentDict[seq_id]) == 0:
                    alignmentDict[seq_id] = aligned_seq
                else:
                    alignmentDict[seq_id] += aligned_seq

  return alignmentDict


class ASR_Node(Node):
  """Cost of states at leaves, given list of conditional states at a
  node.
  
  Costs:  list of dictionnaries {"condition state": cost of transition}
  
  """

  def __init__(self, keep_comments = False, **kw):
    Node.__init__(self, keep_comments, **kw)
    self.__costs = []
    self.__bestseq = ""
    self.__lenseq = 0
    
  def newnode(self):
    return ASR_Node()

  def set_sequence(self, seq):
    self.__costs=[]
    for s in seq:
      self.__costs.append({s:0})
    self.__lenseq = len(seq)
    self.__bestseq = seq
    
  def get_sequence(self):
    if len(self.__bestseq)==0:
      lseq=[]
      for pos in range(self.__lenseq):
        lseq.append(min(self.__costs[pos], key=self.__costs[pos].get))
      self.__bestseq = "".join(lseq)
    return self.__bestseq
  
  def len_seq(self):
    return self.__lenseq

  ### Parsimony part
  def states(self, pos):
    return self.__costs[pos].keys()

  def cost_from(self, state, pos):
    return self.__costs[pos].get(state, 0)
  
  def parsimony(self):
    """
    Computes conditional costs from children of the node.
    Each substitution costs SubsCost.
    """

    self.__compute_forward()
    self.__compute_backward()


  def __get_costs(self, pos):
    return self.__costs[pos]

  def __compute_forward(self):
    """
    Forward recursion of downward parsimony costs.
    Olnly known states are built as they are found.
    """

    if self.is_leaf(): ##Leaf
      return
    else:
      for child in self.get_children():
        child.__compute_forward()

    ### Set here costs
    self.__costs= []
    # if self.is_root():
    #   leaves = self.get_leaves()

    self.__lenseq=self.get_children()[0].len_seq()

    for pos in range(self.__lenseq):
      self.__costs.append({})

      ## check if unknown 

      herecost = self.__costs[pos]

      ## check if unknown 
      herecost = self.__costs[pos] 
      for child in self.get_children():
        for st in child.states(pos):
          if not st in herecost:
            herecost[st] = 0
      if not "*" in herecost:
        herecost["*"] = 0
      
      # if self.is_root():
      #   for leaf in leaves:
      #     for st in leaf.states(pos):
      #       if not st in herecost:
      #         herecost[st] = 0
        
      for child in self.get_children():
        for k in herecost.keys():
          vok = 1000000
          for st in child.states(pos):
            vst = child.cost_from(st, pos) + SubsCost(k,st)
            if vst < vok:
              vok = vst
          herecost[k] += vok
          
  def __compute_backward(self, upcost = []):
    """Backward recursion of upward parsimony costs.
    A dictionnary of up costs are transmitted downward.

    Substitution costs are recomputed intead of stored previously. To
    be ameliorated later.

    """

    if self.is_leaf():
      return
    
    ### First compute all upcosts of children, then recursion, to avoid to double use of down costs
    upcostchildren={}
    for child in self.get_children():
        upcostchildren[child]=[]

    for pos in range(self.__lenseq):
      for child in self.get_children():
        upcostchildpos={}

        ## check if unknown 
        for other in self.get_children():
          if other!=child:
            for st in other.states(pos):
              if not st in upcostchildpos:
                upcostchildpos[st] = 0
        if upcost != []:
          for st in upcost[pos]:
            if not st in upcostchildpos:
              upcostchildpos[st] = 0
        if not "*" in upcostchildpos:
          upcostchildpos["*"]=0
        
        ## then add/min
        for other in self.get_children(): # neighbours
          if other==child:
            continue
          for k in upcostchildpos.keys():
            vok = 1000000
            for st in other.states(pos):
              vst = other.cost_from(st, pos) + SubsCost(k,st)
              if vst < vok:
                vok = vst
            upcostchildpos[k] += vok

        # father
        if upcost!=[]:
          for k in upcostchildpos.keys():
            vok = 1000000
            for st in upcost[pos].keys():
              vst = upcost[pos][st] + SubsCost(k,st)
              if vst < vok:
                vok = vst
          upcostchildpos[k] += vok

        upcostchildren[child].append(upcostchildpos)
        
   ## Then recursion
    for child in self.get_children():
      child.__compute_backward(upcostchildren[child])

    ## Then update self costs
    if upcost!=[]:
      for pos in range(self.__lenseq):
        for k,v in upcost[pos].items():
          if k in self.__costs[pos]:
            self.__costs[pos][k]+=v
          else:
            self.__costs[pos][k]=v+self.__costs[pos]["*"]
        for k in self.__costs[pos]:
          if not k in upcost[pos]:
            self.__costs[pos][k]+=upcost[pos]["*"]
            
if __name__ == "__main__":
  parser = argparse.ArgumentParser()

  parser.add_argument('-t', '--tree', dest='tree', action='store',\
                      required=True,\
                      help='File containing a tree in Newick format')
  parser.add_argument('-a', '--align', dest='align', action='store',\
                      required=True,\
                      help='File containing an alignment in FASTA format')
  parser.add_argument('-o', '--output', dest='output', action='store',\
                      required=True,
                      help='Name of the output Fasta file')
  parser.add_argument('-c', '--costs', dest='costs', action='store',\
                      required=False,
                      help='Name of the costs tsv file. Node numbers as columns & sites as lines.')
  
  args = parser.parse_args()

  tree = ASR_Node()

  tree.read_nf(args.tree, True)
  tree.intersect_ancestral_labels()
  
  align = loadAlignment(args.align)

  ## set sequences
  leaves = tree.get_leaves()
  for leaf in leaves:
    if not leaf.label() in align:
      print("Missing seq " + leaf.label() + " in " + args.align)
    else:
      leaf.set_sequence(align[leaf.label()])

  ## Compute asr
  tree.parsimony()

  ### output all sequences

  fout=open(args.output, "w")
  lnodes = tree.get_all_children()
  for node in lnodes:
    fout.write(">"+node.label()+"\n")

    seq = node.get_sequence()

    fout.write(seq)
    fout.write("\n")

  fout.close()

  ### output all transition costs

  if args.costs:
    fout=open(args.costs, "w")
    lnodes = tree.get_all_children()

    fout.write("[Site]\t")
    fout.write("\t".join([node.label() for node in lnodes]))
               
    lnodes = tree.get_all_children()
    lseq = len(lnodes[0].get_sequence())

    for i in range(lseq):
      fout.write("[%d]"%(i+1))
    
      for node in lnodes:
        if not node.up:
          cost = 0
        else:
          seqi = node.get_sequence()[i]
          upi = node.up.get_sequence()[i]

          cost = SubsCost(upi, seqi)

      fout.write("\t%f"%(cost))

    fout.write("\n")
    fout.close()

