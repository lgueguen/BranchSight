// Variables globales
// ------------------
var treeRoot;
var treeRootIni;
var cladeRoot;
var layout = d3.layout.cladogram();   // Fonction  D3 qui transforme les donnees en cladograme
var nodeWidth = 80;                   // Largeur de l'arbre
var nodeHeight = 20;                  // Hauteur de l'arbre
var hLetter = 15;                     // Hauteur des blocs de couleurs
var hUnderline = 15;                  // Position du soulignage
var fitNormalPoliceSize = 15;         // Taille de la police des aa et dna en mode normal
var fitSmallPoliceSize = 8;           // Taille de la police des aa et dna en mode reduit
var fitPoliceSize = fitNormalPoliceSize;  // Taille de la police des aa et dna 
var stringWidth;                      // Largeur affichage arbrevar svg;                              // La div SVG
var g;                                // Conteneur principal
var diagonal = svgLinks.shoulder;     // cf links.js
var symbol = d3.symbol().size(128);   // Symbol des noeuds
var color = d3.scaleOrdinal(d3.schemeCategory20); // Echelle des couleurs
var options;                          // Options d'affichage
var action;                           // Action associee au clic
var selectedCrossref;                 // Reference croisée choisie (nom de la sequence)
var selectedCrossrefURL;              // URL associe
var crossref;                         // Objet decrivant les differentes ref croisees associees a la famille
var div;
var myStorage = window.localStorage;  // utilisation de localStorage pour consrever les preferences
var displaySeqType = 'AA';            // type de sequence
var hStepNormal = 12;                 // Epaisseur du carré de couleur en mode normal
var hStepSmall = 6;                   // Epaisseur du carré de couleur en mode reduit
var hStepBase = hStepNormal;          // Epaisseur du carré de base
var hStep;                            // Epaisseur du carré ( different si aa ou dna)
var alignmentLetterNormalSpacing = hStepBase - 10;
var alignmentLetterSmallSpacing = 0; 
var alignmentLetterSpacing =  alignmentLetterNormalSpacing;
var psThresholdMaxUp = parseFloat($('#ps-threshold-max-up').val()) || parseFloat(1.0);
var psThresholdMinUp = parseFloat($('#ps-threshold-min-up').val()) || parseFloat(0.9);
var psThresholdMinDown = parseFloat($('#ps-threshold-min-down').val()) || parseFloat(-0.9);
var psThresholdMaxDown = parseFloat($('#ps-threshold-max-down').val()) || parseFloat(-1.0);
var psThresholdMaxUpBgColor = "#FF0000";
var psThresholdMinUpBgColor = "#FFFF00";
var psThresholdMaxDownBgColor = "#0000FF";
var psThresholdMinDownBgColor = "#00FFFF";
var psNormalThresholdBgColor = "#DDDDDD";
var psMissingBgColor = "#faf0e6";
var seq_lg = xmlparser.flatTree(recTree.phyloxml.phylogeny.maxSeqIdLength)*10;
var isDNA = (xmlparser.flatTree(recTree.phyloxml.phylogeny.size.isDNA)=="true");
var isCodon = (xmlparser.flatTree(recTree.phyloxml.phylogeny.size.isCodon)=="true");
var collapseMode = 'collapseOff';
var navWindowColor = '#000000';
var alignmentWindow = $('#alignment-container');
var selposWindow = $('#selpos');
var compactWindow = $('#navigator-container');
var resultsJSON = null;
var resultsShown;
var maxSeqLen = 0;
var updateNodesTimeout = null;
var selectedNode = null;
var globalI = 0;
var leftBoundary = 0;
var rightBoundary = 1;
var firstLoad = 1;
var automaticBoundaries = true;
var precisionWindowStat = 0;
var siteOffset = 0;
var nbBranches = 0;
var nbLeaves = 0;
var nbSites = 0;
var upperThresholdMode = true;
var modStepX = 150;               // Pas de la variation en largeur
var modStepY = 5;                 // Pas de la variation en hauteur
var modeSite = false;             //false : mode fenetre, true : mode site
var selectedSite = 0;

var hslaa = new Object(); // Dictionnaire des aa

  hslaa["K"]= [8,80,50];
  hslaa["R"]= [8,80,50];

  hslaa["A"]= [210,80,50];
  hslaa["F"]= [210,80,50];
  hslaa["I"]= [210,80,50];
  hslaa["L"]= [210,80,50];
  hslaa["M"]= [210,80,50];
  hslaa["V"]= [210,80,50];
  hslaa["W"]= [210,80,50];

  hslaa["N"]= [120,78,45];
  hslaa["Q"]= [120,78,45];
  hslaa["S"]= [120,78,45];
  hslaa["T"]= [120,78,45];

  hslaa["H"]= [180,75,40];
  hslaa["Y"]= [180,75,40];

  hslaa["C"]= [0,67,70];

  hslaa["D"]= [300,55,55];
  hslaa["E"]= [300,55,55];

  hslaa["P"]= [60,100,40];

  hslaa["G"]= [30,75,60];

  hslaa["-"] = [0,0,100];
  hslaa["*"] = [0,0,0];

var flag_redond = 1;
var flag_histogram = 1;
var colordna = new Object();
var colordna_red = new Object();
define_colordna() // Remplit colordna et colordna_red
var  dico_aa = new Object();
var  dico_dna = new Object();
if (!isStorageSupported(myStorage)) {
  console.log("window.localStorage problem");
  $(testnav).append('It seems that <b>window.localStorage</b> is full or is not supported by your browser<br>Please check this <a href="/clear/">here</a> ');
}
else {
  main_display();
}

// Fonction d'affichage principale
// ===============================
function  main_display() {
  $('#redondance').hide();
    cladeRoot = recTree.phyloxml.phylogeny.clade;
  // Creation du SVG
  div = d3.select("body").append("div")  // Definit  la div de classe tooltip
    .attr("class", "tooltip")
    .style("opacity", 1)
    .style('position', 'absolute');
  svg1 = d3.select("#svg1");    // selectionne la div svg1
  g = svg1.append("g")          // ajoute l'element treeAlign au svg1
    .attr("id", "treeAlign");
  var transition = d3.transition();
  margin = {
    top: 100,
    down: 20,
    left: 200,
    right: 300
  }
  // Option d'affichage
  stringOptions = myStorage.getItem("options");
  if (stringOptions === null || stringOptions === undefined || stringOptions === "undefined"){
    stringOptions = "name,species";
  }
  options = stringOptions.split(",");
  // Action associee au click
  action = svgEvents.showBranchResults;
  var stringAction = "View Branch Results";
  myStorage.setItem("action",stringAction);
  stringWidth = myStorage.getItem("width");
  if ((stringWidth === null )||(stringWidth === undefined)||(stringWidth === "undefined")){
    stringWidth = nodeWidth;
    myStorage.setItem("width",stringWidth);
  }
  // Concerne les arbres recphyloxml
  var flatTreeConfig = {
    transferBack :false,
    speciationLoss : false,
    speciationOutLoss : false,
  }
  var _cladeRoot = null;
  var _Parent = null;
  treeRoot = d3.hierarchy(cladeRoot, function(d) {
    return d.clade;
  });  
  addNumberSeqSpec(treeRoot);
  treeRootIni = treeRoot;
  updateLayout(cladeRoot,1,1);
}