
// Autres fonctions utilitaires
// ----------------------------
// Fonction tableau hsl -> chaine de caractere hsl
// -----------------------------------------------
function  hsl2col(hsl) {
  var col = "hsl("+hsl[0]+","+hsl[1]+"%,"+hsl[2]+"%)";
  return col
}
// Definit les couleurs des codons a partir de la couleur des aa
// -------------------------------------------------------------
function define_colordna (){
  var geneticcode = xmlparser.flatTree(recTree.phyloxml.phylogeny.geneticCode);
  geneticcode["---"] = "-";
  var redond_dna = new Object();
  for (d in geneticcode){
    var col = "grey";
    var _aa = geneticcode[d];
    if (_aa != undefined) {
      var _hsl = hslaa[_aa];
      if (_hsl != undefined) {
        col = hsl2col(_hsl)
      }
      if (_aa in redond_dna) {
        redond_dna[_aa].push(d);
      }
      else {
      redond_dna[_aa] = [d];
      }
    }
    colordna[d] = col;
  }
  /*Redefinit les couleurs des codons*/
  var step_satur = 50;
  for (d in redond_dna){
    var _codons = redond_dna[d];
    if (d != "*" && d!= "-") {
      for (i = 0; i < _codons.length; i++) {
        var _hsl =  hslaa[d];
        var new_hsl = [_hsl[0],_hsl[1],_hsl[2] + Math.floor(i* step_satur / _codons.length)];
        //console.log(_codons[i]+ ": avant "+ colordna[_codons[i]]+ ", apres "+ hsl2col(new_hsl));
        colordna_red[_codons[i]] = hsl2col(new_hsl);
      }
    }
  }
}
// Renvoie la couleur associee au taux
// -----------------------------------
function getColorFromResult(res) {
  if (upperThresholdMode) {
    if (res >= psThresholdMaxUp) {
      return psThresholdMaxUpBgColor;
    } else if (res >= psThresholdMinUp) {
      return psThresholdMinUpBgColor;
    } else {
      return psNormalThresholdBgColor;
    }
  } else {
    if (res <= psThresholdMaxDown) {
      return psThresholdMaxDownBgColor;
    } else if (res <= psThresholdMinDown) {
        return psThresholdMinDownBgColor;
    } else {
        return psNormalThresholdBgColor;
    } 
  }
}
  
// Renvoie l'opacite associee au taux
// ----------------------------------
function getOpacFromResult(res) {
  if (upperThresholdMode) {
    if (res >= psThresholdMaxUp) {
      return 1.0;
    } else if (res >= psThresholdMinUp) {
      return 0.5;
    } else {
      return 0.0;
    }
  } else {
    if (res <= psThresholdMaxDown) {
      return 1.0;
    } else if (res <= psThresholdMinDown) {
      return 0.5;
    } else {
      return 0.0;
    }
  }
}
  
// Fonction de calcul de la statistique maximale dans les résultats
// ----------------------------------------------------------------
function maxStatJSON(jsonResultsList) {
    var currentmax = jsonResultsList[0];
    jsonResultsList.forEach(function(res) {
      if (res > currentmax) {
        currentmax = res;
      }
    });
    return currentmax;
    }
    
// Fonction de calcul de la statistique minimale dans les résultats
// ----------------------------------------------------------------
function minStatJSON(jsonResultsList) {
    var currentmin = jsonResultsList[0];
    jsonResultsList.forEach(function(res) {
      if (res < currentmin) {
        currentmin = res;
      }
    });
    return currentmin;
    }
    
// Couleur des noeuds et branches de l'arbre selon le taux
// -------------------------------------------------------
function colorByRate(rate) {
    if (upperThresholdMode) {
        return hsl2col([0, 100, 85 * (1  - rate**0.3)]);
    } else {
        return hsl2col([250, 100, 85 * (1  - rate**0.3)]);
    }
}
 
// Fonction de test du support local.storage
// ------------------------------------------
// https://michalzalecki.com/why-using-localStorage-directly-is-a-bad-idea/
function isStorageSupported(storage) {
  try {
    const key = "__some_random_key_you_are_not_going_to_use__";
    storage.setItem(key, key);
    storage.removeItem(key);
    return true;
  } catch (e) {
    return false;
  }
  }
  
