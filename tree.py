#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Newick format parser module.

Manipulates trees in the Newick format. This format is mainly used to
describe phylogenetic binary trees but can have a much wider use. Some
details on the format can be found here:

http://evolution.genetics.washington.edu/phylip/newicktree.html

Examples of Newick format:

* (Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;

* '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'

Note that bracked delimited nested and unnested comments are ignored.

"""

import string
import copy
import re
from math import log
import numpy

#######################################################################
#######################################################################
########  Node class: attributes and instance creation
def find_block(s, delimin, delimout):
    """ Finds in s the first block delimited with delimiters delimin-delimout.
Returns the indexes of these delimitors in s.
    """

class Node(object):
  """Node is the smallest unit of definition for a tree.

    A Node is linked to its father-Node and its children-Nodes: it
    defines a sub-tree for which it is the root. All connected Nodes
    (including leaves and root) make up for the whole tree.
  """

  def __init__(self, keep_comments=False, **kw):
    """Create a Node.
    
    Keyword argument to build directly the Node from a Newick string: [newick=string].
    Examples of Newick format:
    * '(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268,Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;'
    * '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);' 
    """
    
    self.__l=None         # length of the branch to the father-Node
    self.__lab=""         # Node label
    self.__psize=None       # population size up the Node
    self.__comment=None     # nested comment attached to Node
    self.__father=None      # father Node instance
    self.__children=[]      # list of child Node instances
    if 'newick' in kw:
      self.parser(s=kw['newick'], keep_comments=keep_comments)
    elif 'fic' in kw:
      self.read_nf(a=kw['fic'], keep_comments=keep_comments)

      
#####################################################
############## methods for access to object atributes:

  def newnode(self):
    """class-specific instance generator for construction of trees of Nodes"""
    return Node()

  def __getitem__(self,n):
    """Return the Node with this name."""
    if self.__lab==n: return self
    for s in self.__children:
      if s.label()==n:
        return s
      else:
        t=s[n]    # <==> t = s.__getitem__(n)   (recursive function)
        if t:
          return t
    return None
        
  def __str__(self):
    """Return printable string in Newick format of the sub-tree defined by the Node."""
    return self.newick()
              
  def __iter__(self):
    return self.generator()
              
  def generator(self):
    children = self.get_all_children()
    for child in children:
      yield child

  def lg(self):
    """Return the length of the edge to the father."""
    return self.__l

  def label(self):
    """Return the label."""
    return self.__lab

  def psize(self):
    return self.__psize
  
  def comment(self):
    """Return the comment."""
    return self.__comment
              
  def get_all_children(self):
    """Return the list of all nodes below the Node, including itself"""
    a=[]
    if self.__children!=[]:
      for i in self.__children:
        a+=i.get_all_children()
    a+=[self]
    return a

  def get_sorted_children(self):
    """Return the list of all nodes below the Node, including itself, ordered by decreasing depth (i.e. increasing node distance from root)."""
    a = self.get_all_children()
    a.sort(key=lambda x: x.depth())
    return a
              
  def get_all_parents(self):
    """Return the list of all nodes above the Node, excluding itself, ordered by increasing depth (i.e. decreasing node distance from root)."""
    lparents = []
    f = self.go_father()
    while f:
      lparents.append(f)
      f = f.go_father()
    return lparents
              
        
                
  #####################################################
  ############## methods for editing object atributes:
    
  def __setattr__(self, name, value):
     self.__dict__[name] = value
          
  def set_lg(self,l):
    """Set the length of the edge to the father to $2 if it is >=0."""
    lg = float(l)
    if lg>=0:
      self.__l = lg
      
  def set_psize(self, ps):
    boot = float(ps)
    if boot>=0:
      self.__psize = boot

  def add_label(self, lab):
    self.__lab = lab
          
  def edit_label(self, label, mode="write"):
    """edits the label attached to the node"""
    if mode=="write" or self.__lab=="":
      self.__lab = str(label)
    elif mode=="append":
      self.__lab += ', '+str(label)
    else:
      raise ValueError    
          
  def edit_comment(self, comment, mode="write"):
    """edits the comment attached to the node"""
    if mode=="write" or not self.__comment:
      self.__comment = str(comment)
    elif mode=="append":
      self.__comment += ', '+str(comment)
    else:
      raise ValueError

        
  def restrictToLeaves(self, lleaves):
     """returns a copy of the tree restricted to the input leaf set"""
     c = copy.deepcopy(self)
     # maps to the common ances tor of all leaves
     mrca = c.map_to_node(lleaves)
     # extract this subtree
     if mrca == c:
       st = c
     else:
       st = c.pop(mrca.label())
     # remove the other leaves
     for leaf in st.get_leaf_labels():
       if leaf not in lleaves:
         st.pop(leaf)
     return st

  def complete_label(self, prefix='N', labels=[], force=True):
    """gives label to an internal node given a set of pre-existing labels in the tree"""
    children = self.get_sorted_children()
    for c in children:    # first builds a list of existing labels to avoid redundancy of new names
      cl = c.label()
      if (cl!="") and (not cl in labels):
        labels.append(cl)
        
    n = 1
    if self.label()=="" or force:
      lab = "%s%d"%(prefix, n)
      while lab in labels:
        n += 1
        lab = "%s%d"%(prefix, n)
        self.add_label(lab)

  def complete_internal_labels(self, prefix='N', labels=[], force=False):
    """recursively gives label to internal nodes that lack one given a set of pre-existing labels in the tree"""
    children = self.get_sorted_children()
    for c in children:    # first builds a list of existing labels to avoid redundancy of new names
      cl = c.label()
      if (cl!="") and (not cl in labels):
        labels.append(cl)

    n = 1
    for c in children:
      try:
        label = float(c.label())
      except ValueError:
        label = c.label()
      if label=="" or force or type(label) is float:
        lab = "%s%d"%(prefix, n)
        while lab in labels:
          n += 1
          lab = "%s%d"%(prefix, n)
        c.add_label(lab)
        labels.append(lab)
        n += 1    
                      
#####################################################
############## methods for access to object's child atributes:

  def get_children(self):
    """Return the list of direct children of the Node"""
    return self.__children

  def nb_children(self):
    """Return the number of children."""
    return len(self.__children)


  def children_labels(self):
    """Return the list of direct child labels."""
    llc = []
    for c in self.__children:
      llc.append(c.label())
    return llc

  def children_comments(self):
    """Return the list of child comments."""
    lcc = []
    for c in self.__children:
      lcc.append(c.comment())
    return lcc
              
  def get_children_labels(self, sorted=True):
    """Return the list of labels of all nodes below the Node, including itself"""
    if sorted:
      children = self.get_sorted_children()
    else:
      children = self.get_all_children()
    l =[]
    for c in children:
      l.append(c.label())
    return l
              
  def get_children_node_ids(self):
    lid = []
    children = self.get_all_children()
    if children[0].nodeid() == None:
      self.complete_node_ids()
    for child in children:
      lid.append(child.nodeid())
    return lid
              
#####################################################
############## methods for access to object's leaf atributes:
        
  def nb_leaves(self):
    """Return the number of leaves under this node."""
    return len(self.get_leaves())  
              
  def get_leaf_labels(self, comments=False):
    """Return the list of labels of the leaves defined by the Node."""
    
    a=[]
    if self.__children!=[]:
      for i in self.__children:
        a+=i.get_leaf_labels(comments)
    else:
      text = self.__lab
      if comments==True and self.__comment:
        text += '[%s]'%self.__comment
        a+=[text]
      else:
        a+=[text]
    return a
        
  def get_leaves(self):
    """Return the list of leaves defined by the Node."""
    
    a=[]
    if self.__children!=[]:
      for i in self.__children:
        a+=i.get_leaves()
    else:
      a+=[self]
    return a
          
  
#####################################
############## file input methods:    
  
  def read_nf(self, a, keep_comments=False):
    """Read the $1 file containing a unique tree in Newick format and builds the Node from it."""

    try:
      f=open(a,'r')
      l=f.readline()
      f.close()
      self.parser(l, keep_comments)
    except OSError:
      print("Invalid file name: " + str(a))
    return

#####################################
############## file output methods:

        
  def newick(self, NHX=False):
    """Returns (extented) Newick string of the sub-tree defined by the Node.

    Follows the (extended) Newick/New Hampshire (NHX) standard for coding trees. Supports bracketed comments, filled with tree.Node attribute as specified by 'comment'.
    Newick format example:
                '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
    """
    s=''
    end=''
    if self.__father==None:
      end=';'
    i=0
    while i<len(self.__children):
      if i==0:
        s+='('
      s+=self.__children[i].newick(NHX)
      if i!=(len(self.__children)-1):
        s+=','
      else:
        s+=')'
      i+=1
      flag=False

    if self.is_leaf():
      s+=self.__lab
      flag=True

    if self.__l!=None:
      s+=":"+str(self.__l)

    if NHX:
      if self.__comment:
        lc=[self.__comment]
      else:
        lc=[]
      if self.__psize:
        lc.append("W=%f"%self.__psize)
      if self.__lab and not self.is_leaf():
        lc.append("ND=%s"%self.__lab)

      if len(lc)!=0:
        s+='[&&NHX:%s]'%(":".join(lc))

    s+=end

    return s

  def matrix(self, leavesOnly=False, header=True):
    """Return additive distances matrix of the leaves."""
    
    d=self._matrix(leavesOnly=leavesOnly)
    n=d.keys()
    n.sort() #order lines and columns
    s=''
    if header:
      s += '\t'.join(n) + '\n'
    for i in n:
      if header:
        s += i+'\t'
      for j in n:
        s+=str(d[i][j])+'\t'
        s+='\n'
    s=string.strip(s)
    return s
              
  def dichotomic_list(self, sort=None):
    """returns an ordrer list of node labels in the tree, with father node loacted between its two children.
    
    by default order the nodes by decrerasing size; specifying 0 or 1 makes them appear as read in the tree.
    """
    l = []
    children = self.__children
    if children :
      if len(children)==2:
        if not sort:
          if children[0].nb_leaves() <= children[1].nb_leaves():
            sort = 1
          else:
            sort = 0        
          l += children[0-sort].dichotomic_list(sort=sort)
          l.append(self.__lab)
          l += children[1-sort].dichotomic_list(sort=sort)
    else:
      l.append(self.__lab)
    return l
              
        
  def arborescence_ASCII(self, depth=0, uncleBranches=[0], lastSon=False, out="print"):
    """prints a multi-line ASCII drawing of the tree, in a hirearchical arborescence form.
    
    If 'out' argument is "print" (default), the ouput is written to STDOUT ; else, the output is APPENDED to file named 'out'.
    """
    if out=="print":
      from sys import stdout
      Out = stdout
    else :
      Out = open(out, 'a')
      
    branchlength = self.go_root().max_leaf_depth()
    lenpat = (max(uncleBranches))
    pattern = ["    "]*(lenpat+2)
    for branch in uncleBranches:
      pattern[branch+1] = "|   "
      Out.write( "".join(pattern)+'\n' )
    if self.is_leaf():
      Out.write( "".join(pattern[:-1])+"+----"+"----"*(branchlength-lenpat)+" "+self.label()+'\n' )
      if not out=="print":
        Out.close()
    else:
      # Out.write( "".join(pattern[:-1])+"+---+"+"----"*(branchlength-lenpat)+" "+self.label()+'\n' )
      # if not out=="print":
      #   Out.close()
      children = self.get_children()
      if lastSon:
        del uncleBranches[-1]
      for child in children[:-1]:
        child.arborescence_ASCII(depth=depth+1, uncleBranches=uncleBranches+[depth], out=out)
        children[-1].arborescence_ASCII(depth=depth+1, uncleBranches=uncleBranches+[depth], lastSon=True, out=out)
              
    
  
#####################################
########### distance measure methods:

  def depth(self):
    """Return the depth of this Node as its 'distance' (in number of branches) to the root.
                
    WARNING: Branch lengths are not used.
    """
    
    d=0
    if self.__father:
      d=self.__father.depth()+1
    return d
              
  def max_leaf_depth(self):
    """returns the maximum of all leaf depths (from the tree root) under the Node."""
    mld = 0
    for leaf in self.get_leaves():
      if leaf.depth() > mld :
        mld = leaf.depth()
    return mld

  def distance_root(self, nullBranchesAsZeroLength=False):
    """Return the distance (sum of branch lengths) from this Node to the root."""
    if not self.__l:
      if nullBranchesAsZeroLength:
        l = float(0)
      else:
        raise ValueError("No branch length at node %s"%self.__lab)
    else:
      l = self.__l
    if self.__father:
      return l+self.__father.distance_root(nullBranchesAsZeroLength=nullBranchesAsZeroLength)
    else:
      return l
                      

#####################################
############## miscellaneous methods:
  def get_parents(self):
    lparents = []
    f = self.go_father()
    while f:
      lparents.append(f)
      f = f.go_father()
    return lparents        
              
  def get_parent_labels(self):
    """Return the list of labels of all nodes above the Node, excluding itself, ordered by increasing depth (i.e. decreasing node distance from root)."""
    lparents = []
    f = self.go_father()
    while f:
      lparents.append(f.label())
      f = f.go_father()
    return lparents  

  def get_dicCladeToLeaves(self, prefix):
    """Returns a dictionary of all leaves below each node under Node"""
    self.complete_internal_labels(prefix)
    dic = {}
    children = self.get_sorted_children()
    for c in children:
      dic[c.label()] = c.get_leaf_labels()
    return dic
        
  def __imul__(self, x):
    """ Recursively multiplies the length of the branches that are under this node by factor x (>0)."""
    if x>0:
      for i in self.__children:
        i.__l*=x
        i*=x
    return self
              
  def __idiv__(self, x):
    """ Recursively divides the length of the branches that are under this node by factor x (>0)."""
    if x>0:
      for i in self.__children:
        i.__l/=x
        i/=x
    return self
              
########################################
### Topology / taxonomic comparison methods:
  
  def is_bifurcated(self):
    if len(self.__children) == 2:
      return True
    else:
      return False

  def is_leaf(self):
    """returns bollean stating if Node is a terminal node/leaf"""
    if self.get_leaves() == [self]:
      return True
    else:
      return False
                      
  def is_child(self, node1):
    """returns boolean stating if the Node is below node1 in the tree"""
    if not self.__father:
      return False
    if self.__father==node1:
      return True
    else:
      return go_father().is_child(node1)      

  
  def is_brother(self, node1):
    """returns bollean stating if the Node is neighbor of node1 in the tree"""
    try:
      bro = self.go_brother()
    except ValueError:
      raise IndexError("Node %s has no brother"%self.__lab)
    if node1==bro:
      return True
    else:
      return False
                
#####################################     
### 'walking along the tree' methods:

  def go_root(self):
    """Return Root Node."""  
    if self.__father:
      self=self.__father.go_root()
    return self

  def is_root(self):
    """Test Root Node."""  
    return not self.__father

  def go_father(self):
    """Return father Node of this Node."""
    return self.__father
              #      raise IndexError, 'Already at root. No more father to go to.'
    

  def go_child(self,n):
    """Return $1 (number) child of this Node."""
    if self.__children and n<len(self.__children):
      self=self.__children[n]
    else:
      raise IndexError('Already at leaf. No more children to go to.')
    return self
        
  def go_brother(self):
    """return brother node (other child from biffurcated father node)"""
    fat = self.go_father()
    if fat:
      c = fat.get_children()

      if len(c)==2:
        if self == c[0]:
          return c[1]
        elif self == c[1]:
          return c[0]
        else:
          raise ValueError("Node %s not in his father's child list"%self.__lab)
      else:
        raise ValueError("Node %s's father has %s child(ren)"%(self.label(),len(c)))
    else:
      raise ValueError("Node %s has no father"%self.__lab)
          
  def go_uncle(self):
    """return uncle node (other child from biffurcated grandfather node)"""
    fat = self.go_father()    
    if fat:
      uncle = fat.go_brother()
      return uncle  
          
  def path_to(self, n, returnLabels=False):
    """returns list of nodes on the path from the Node to node $1 (including both)"""
    if self.go_root() != n.go_root():
      raise IndexError("Nodes are not parts of the same tree.")
    elif self==n:
      if returnLabels:
        return [self.label()]
      else:  
        return [self]
    else:
      if n in self.get_all_children():
        path = n.path_to(self, returnLabels=returnLabels)
        path.reverse()
        return path
      elif self in n.get_all_children():
        f = self.go_father()
        if returnLabels:
          return [self.label()] + f.path_to(n, returnLabels=returnLabels)
        else:  
          return [self] + f.path_to(n, returnLabels=returnLabels)
      else:
        f = self
        path = []
        while (n not in f.get_all_children()):
          if returnLabels:
            path.append(f.label())
          else:
            path.append(f)  
            f = f.go_father()
        return path + f.path_to(n, returnLabels=returnLabels)
            
  def lengths_on_path_to(self, n, excludeLeaves=False):
    """returns list of branch bootstraps on the path from the Node to node $1"""
    treePath = self.path_to(n)
          #    for node in treePath:
          #      print node.label(),
    llg = []
    for i in range(1, len(treePath)):
      if treePath[i-1].depth() > treePath[i].depth():
        if not (excludeLeaves & treePath[i-1].is_leaf()):
          llg.append(treePath[i-1].lg())
      elif treePath[i].depth() > treePath[i-1].depth():
        if not (excludeLeaves & treePath[i].is_leaf()):
          llg.append(treePath[i].lg())
      else:
        raise IndexError('Neighbour nodes %s and %s should be at different depth'%(treePath[i-1].label(), treePath[i].label()))
    return llg    
        
  def hierachical_path(self, labels=[], labprefix='N'):
    root = self.go_root()
    if not labels: labels = root.get_children_labels()
    sufix = ""
    if self!=root:
      path = root.path_to(self)
      for node in path:
        if node.label() in labels:
          sufix += '.'+node.label().lstrip(labprefix)
    sufix += '.'+self.label().lstrip(labprefix)
    return sufix

#####################################
#################### parsing methods:
  
  def clean(self,s):
    """Clean the string in Newick format.
    
    Bracked delimited nested and unnested comments are eliminated.
    """
    s=s.strip()
    if s[len(s)-1]==';':
      if s.count('(')!=s.count(')'):
        raise ValueError("Opening parenthesis do not match closing parenthesis.")
      else:
        brackopen=s.count('[')
        brackclose=s.count(']')
        if brackopen!=0 or brackclose!=0:
          if brackopen==brackclose:
            open=s.find('[')
            x=1
            for i in range(open+1,len(s)):
              if s[i]=='[':
                x+=1
              elif s[i]==']':
                x-=1
              if x==0:
                break
            s=s[:open]+s[i+1:]
            s=self.clean(s)
          else:
            raise ValueError("Opening brackets do not match closing brackets.")
    else:
      raise ValueError("Missing ';' at end.")
    return s

  def read_commented_lab(self, s):
    """deals with nested comments located next to labels"""
    lbrack = s.find('[')
    rbrack = s.find(']')
    if lbrack==-1 and rbrack ==-1:
      self.__lab = s        # clean case
    else:
      lp=s[lbrack:rbrack].split(":")
      nn=[]
      for n in lp:
        if n[:2]=="W=":
          self.__psize=float(n[2:])
        elif n[:3]=="ND=":
          self.__lab=n[3:]
        else:
          nn.append(n)
      self.__comment=":".join(nn[1:])
      
  def _parser(self,s,nodeId=[]):
    """Should not be directly used. Use parser() instead.

    Length of nodeId is node Id

    """
    
    x=0
    pcomas=[] #cutting points for nodes of the same depth
    for i in range(len(s)):
      if s[i]=='(':
        x+=1
      elif s[i]==')':
        x-=1
      elif x==0 and s[i]==',':
        pcomas+=[i]

    if pcomas==[]: # Comment level
      pare=s.rfind(")")
      semicol=s.find(':',pare+1)
      bracke=s.rfind(']')
      brackd=s.rfind('[',0,bracke)

      if bracke>pare:  # there are comments
        self.read_commented_lab(s[brackd:bracke+1])
      else:
        brackd=len(s)

      if semicol!=-1 and semicol<brackd : # there is a length
        try:
          self.__l=float(s[semicol+1:brackd])
        except ValueError:
          raise ValueError("Incorrect branch value -> must be numerical.")
      else:
        semicol=len(s)

      if not self.__lab:
        self.__lab=s[pare+1:min(semicol,brackd)]

      try:
        x=float(self.__lab)
        self.__lab="N"+str(len(nodeId))
      except ValueError:
        pass
      
      if s[0]=="(":  #fetch children
        nodeId.append(0)
        self._parser(s[1:pare],nodeId)

      return

    ### children level
    if len(pcomas)==0:
      return

    pcomas.append(len(s))
    
    d=0
    for pcom in pcomas:
      child=s[d:pcom]
      d=pcom+1
      nnode = self.newnode()
      nnode.__father=self
      nodeId.append(0)
      nnode._parser(child,nodeId)
      self.__children.append(nnode)

        
  def parser(self, s, keep_comments=False):
    """Fill the Node's attributes from parsing $1 string.
    
    Follows the Newick standard for coding trees.
    Bracked delimited nested and unnested comments are ignored.
      
    Newick format example:
         '((A:0.1,B:0.2,C:0.1)ABCnode:0.2,(D:0.4,E:0.1)99:0.1);'
    
    WARNING: all interior nodes (if given) should be strings ->
    numerical values will be interpreted as the node's bootstrap
    value.
    """
    
    if not keep_comments:
      s=self.clean(s)

    s="".join(s.split()).strip(';') # remplace s=s.strip(';')
    # if no branche length defined at root, insert zero branch legnth value
    parent = s.rfind(')')
    if not ':' in s[parent+1:]:
      s += ":0"
                          
    if s:
      self._parser(s)
      #      if self==self.go_root():
      #        self.__l = 0 # sets distance above root at 0.
    else:
      raise ValueError("This should not have happened! Check behind your back.")


                      
  def add_child(self, newchild, silent=True):
    self.__children.append(newchild)
    if silent==False:
      return("%s -< %s"%(self.label(), newchild.label()))
            
  def rm_child(self, child, silent=True):
    try:
      self.__children.remove(child)
    except ValueError:
      return("rm_child : Except no child %s in %s "%(self.label(), child.label()))
    if silent==False:
      return("%s -/- %s"%(self.label(), child.label()))
