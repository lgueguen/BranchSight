// Get the ponderated rate of positive sites within the indicated boundaries
//   Given a preci
// --------------------------------------------------------------
function getPositiveRateJSON(textResultsList, leftB, rightB) {
  var jsonResultsList = JSON.parse(textResultsList);

  var mid = (leftB + rightB)/2;
  var nbPos = 0.0;
  var sPos = 0.0;
  var fact = 1.0;

  jsonResultsList.forEach((value, site) => {
      if (site >= leftB && site <= rightB) {
          if (precisionWindowStat != 0) {
              fact = Math.exp(-((site - mid) * (site - mid) * precisionWindowStat));
          }
          if (upperThresholdMode) {
               if (value >= psThresholdMaxUp) {
                   nbPos += fact;
               }
          } else {
               if (value <= psThresholdMaxDown) {
                   nbPos += fact;
               }
          }
          sPos += fact;
      }
  })
  return nbPos/sPos;
}

  
// Get the value  at a  site 
// --------------------------------------------------------------
function getValueJSON(textResultsList, position_x) {
    var jsonResultsList = JSON.parse(textResultsList);
    valuePos = 0;
    jsonResultsList.forEach((value, site) => {
        if ( site == position_x) {
            valuePos = value;
        }
        
    });
    //console.log("value = "+valuePos);
    return valuePos;
}

// Get the transition on a ste & branch"
  
  function getTransitionJSON(node, target, position_x) {
      if (displaySeqType == "Nuc") {
        if (isCodon)
          {
              var targetseq = target.dnaAlign[3*position_x]+target.dnaAlign[3*position_x+1]+target.dnaAlign[3*position_x+2];
              var nodeseq = node.dnaAlign[3*position_x]+node.dnaAlign[3*position_x+1]+node.dnaAlign[3*position_x+2];
          } else {
              var targetseq = target.dnaAlign[position_x];
              var nodeseq = node.dnaAlign[position_x];
          }
      } else {
          var targetseq = target.aaAlign[position_x];
          var nodeseq = node.aaAlign[position_x];
      }


      return nodeseq + ":" + targetseq;
  }

// Reload 
// ------
function reloadGlobalResults() {
  resultsJSON = JSON.parse(xmlparser.flatTree(recTree.phyloxml.phylogeny.global_results.results));
  updateLayout(cladeRoot);
}