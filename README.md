### Web interface for the visualization of branch-site positive selection

### Université Claude Bernard Lyon 1

# Introduction

This interface was developed to visualize data pertaining to branch-site information in genes.
As long as the data contains a Newick formatted tree file, a FASTA alignment file and in option a results file
organized with branch ids in columns and site in lines, the Python script used to generate
the XML file used by the interface will work.

## Prerequisites for a local installation

- Python
    - biopython, io
    - lxml, xml
      
- Node.js
    - npm

# Usage

## Uploading your files

On the homepage, upload your three files on the required fields:

- Newick-formatted tree
- FASTA alignment
- In option statistical results

 | sites | 8 | 0 | 7 | 1 | 6 | 
 | ----- | --- | --- | --- | --- | --- |
 | [1] | 0.323695788103582 | 0 | 0 | 0 | 0 | 
 | [2] | 0.445776741626087 | 0.44878335841538 | 0.408887020201172 | 0 | 0 | 
 | [3] | 0.366470636875379 | 0.358224955556533 | 0.34148165722248 | 0.32657043259573 | 0.31074948178259 | 
 | [4] | 0 | 0.387294788874982 | 0.318187116484806 | 0.673384398569693 | 0 | 
 | [5] | 0 | 0.345311018159442 | 0.326011772922341 | 0 | 0 | 
 | [6] | 0.329950491998222 | 0.328879470064924 | 0.319574257399895 | 0.304528333261124 | 0.302141100257283 | 
 | [7] | 0.326998331811734 | 0.323310067236089 | 0.313183991818365 | 0 | 0 | 
 | [8] | 0.361985904669695 | 0 | 1.08536617952866 | 0.335431204836716 | 0.308539664507734 | 
 | [9] | 0.337258050881582 | 0.328157824123904 | 0.318720486268041 | 0.305833875355349 | 0.302609308226057 | 
 | [10] | 0.44260781332429 | 0 | 0.384769524486936 | 0 | 0 | 
 | [11] | 0.424406720630182 | 0 | 0.363677085259181 | 0.364045980870431 | 0.321583080347298 | 
 | [12] | 0.322704805534562 | 0 | 0 | 0 | 0 | 

Once you are done, just click the *Display* button and the visualization will be available.

## Visualization page

![Image formulaire](BranchSight.png)

The layout is as follows:

- Tree on the left
- Results histogram and alignment on the right

Buttons include:
- Toggle **Histogram** to show/hide the histograms
- **Precision** to give more importance to sites near the middle of the window: if
the middle of the window is at position $m$, and the precision is $p$,
the importance of a site $s$ inside the visible window is proportional
to $\exp(p . (s-m)^2)$. When $p=0$ (default value), all sites in the
window have equal importance.
- Several tree manipulation options through clicking on nodes:
    - **View Branch Results**
    - **Switch children** to switch orders in tree branches
    - **Collapse/Expand** to hide/show clades
    - **Sub Tree/Upper Tree** to show only a subtree
- **Expand Tree** shows all hidden branches
- **Reset Tree** resets tree with all starting nodes
- **Width, Logarithmic, Height** buttons change appearance sizes of the tree
- **Alignment Type** In case of coding sequences, chooses which sequence type to display (amino acids / codons)
- **Redundancy** with codons toggles fading of synonymous codons
- **Export** to export Tree/Alignment in svg format
- **Thresholds** to choose between thresholds for coloring:
    - "<" will look at values below the lower threshold
    - ">" will consider values above the upper threshold

