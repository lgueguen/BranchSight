// Fonction recursive qui renvoie les noeuds ancestraux communs a
// une liste de feuilles
// ---------------------------------------------------------------
function recuAncestors (todo, done) {
  if (todo.length > 1) {
    var n1 = todo.shift();
    var n2 = todo.shift();
    var anc = n1.path(n2);
    var maxHeight=0;
    var n3=n1;
    anc.forEach(function (d ){
      if (d.height > maxHeight) {
        maxHeight = d.height;
        n3=d;
      }
    });
    todo.splice(0, 0,n3);
    done.push(n1);
    done.push(n2);
    anc = n1.path(n3);
    anc.forEach(function (d ){
      done.push(d)
    });
    anc = n2.path(n3);
    anc.forEach(function (d ){
      done.push(d)
    });
    return recuAncestors (todo, done)
  }
  else {
    var n1  = todo.shift();
    done.push(n1);
    return (done, [])
  }
}

// Fonction de calcul de la position  d'un noeud a partir des longueurs
// de branch_length
// --------------------------------------------------------------------
function phylo(n) {
  var dist = 0.0;
  if (n && n.data.branch_length != null) {
    var p = n.parent;
    if (logBranchLength) {
      dist = Math.log(parseFloat(n.data.branch_length)+1.1) + phylo(p);
    } else {
      dist = parseFloat(n.data.branch_length) + phylo(p);
    }
    return dist;
  }
  return dist;
}
// Fonction de calcul de la longueur max d'un arbre
// ------------------------------------------------
function getmaxlength(treeRoot,max) {
  treeRoot.each(function (d) {
    var phylodist = phylo(d);
    if (phylodist > max) {
      max = phylodist;
    }
  });
  return max;
}
// Fonction pour utiliser les longueurs de branches lors de l'affichage
// d'un arbre
// --------------------------------------------------------------------
function phylogeny(treeRoot,offset) {
  treeRoot.each(function (d) {
    var phylodist = phylo(d);
    d.x = phylodist*offset;
  });
}

// Fonction ajoute nb sequences et especes
// ---------------------------------------
function  addNumberSeqSpec(treeRoot) {
  treeRoot.eachAfter(function (d) {
    var specs = [];
    if (!d.children) {
      d.data.nbseq =  {nbSeq : 1};
    }
    else {
      var nbseq = 0;
      var fils = d.children;
      fils.forEach(function (d){
        nbseq = nbseq + d.data.nbseq.nbSeq;
      })
      d.data.nbseq =  {nbSeq :nbseq};
    }
  });
}
// Fonction qui ouvre les noeuds d'une certaine profondeur
// ----------------------------------------------------------
function expandTree(treeRoot) {
  treeRoot.each(function (d) {
    if (d.data.nodeinfo) {
      if (d.data.nodeinfo.status === "collapsed") {
        if (d.data._clade) {
          d.data.clade = d.data._clade;
          d.data._clade = null;
          d.data.nodeinfo = {status : "extended"};
        }
      }
    }
  });
}
// Fonction qui reinitialise l'arbre
// -------- -------------------------
function resetTree(treeRoot) {
  var fils = treeRoot.descendants();
  fils.forEach(function (d) {
        if (!d.data.clade) {
      d.data.clade = d.data._clade;
      d.data._clade = null;
    }
    if (d.data.nodeinfo) {
      if (d.data.nodeinfo.status === "collapsed") {
        d.data.nodeinfo = {status : "extended"};
      }
    }
    // if   (!d.data.clade) {
    //   d.data.clade = d.data._clade;;
    // }
  })
}
