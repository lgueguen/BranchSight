function hideTooltip() { $('.tooltip').hide(); }
function showTooltip() { $('.tooltip').show(); }
function writeTooltip(text) { $('.tooltip').text(text); }


// Navigation on the small graph
// -----------------------------
function scrollAlignment(event) {
  // Get the site's number (mouse position in pixels) and navigate to this site in the large window
  var siteNumber = event.pageX + compactWindow.scrollLeft() -  window.innerWidth *0.3; // factor by 0.5 because the div are located at  30% 
  var leftScroll = 50 + seq_lg + 10 + siteNumber*hStep - alignmentWindow.innerWidth()/2;
  var rightScroll = 50 + seq_lg + 10 + siteNumber*hStep + alignmentWindow.innerWidth()/2;
  var siteStart = Math.trunc((leftScroll - (60 + seq_lg))/hStep)+2;
  var siteEnd = Math.trunc((rightScroll - (60 + seq_lg))/hStep)+2;
  alignmentWindow.scrollLeft(leftScroll, 0);
  selposWindow.scrollLeft(leftScroll, 0);
  if (automaticBoundaries) {
    leftBoundary = siteStart;
    rightBoundary = siteEnd;
  }
}
// Keep the navigation window visible after a branch change
// --------------------------------------------------------
function updateScroll() {
  alignmentWindow
    .scrollLeft(alignmentWindow.scrollLeft()+1)
    .scrollLeft(alignmentWindow.scrollLeft()-1);
  selposWindow
      .scrollLeft(selposWindow.scrollLeft()+1)
      .scrollLeft(selposWindow.scrollLeft()-1);
}



// Actions sur les boutons
// -----------------------
$('#selectCrossref').on('click', function(event) {
  selectedCrossref=$( "#selectCrossref option:selected" ).text();
  selectedCrossrefURL=crossref[selectedCrossref].url;
  myStorage.setItem("crossref",selectedCrossref);
});
$('#selectCrossref').ready(function() {
  });
$('#action').ready(function() {
  switch (myStorage.getItem("action")) {
    case "View Branch Results":
      $('#branch-res').css('color','#FF0000');
      break;
    case "Switch children":
      $('#switch').css('color','#FF0000');
      break;
    case "Collapse/Expand":
      $('#collapse').css('color','#FF0000');
      break;
    case "Sub Tree/Upper Tree":
      $('#subtree').css('color','#FF0000');
      break;
    default:
  }
});


$('#stats').on('click', function(event, dropdownData) {
  event.preventDefault();
  var target = event.toElement || event.relatedTarget || event.target || function () { throw "Failed to attach an event target!"; }
  myStorage.setItem("stats",target.text);
  //console.log("DEBUG STATS "+target.text);
  $('#button_stats').text(target.text);
  switch (target.text) {
    case "Proportion":
      stats_type = "proportion";
      $('#branch-stats').css('color','#000000');
      $('#rate').css('color','#FF0000');
      $('#extreme').css('color','#000000');      
      break;
    case "Extreme":
      stats_type = "extreme";
      $('#branch-stats').css('color','#000000');
      $('#rate').css('color','#000000');
      $('#extreme').css('color','#FF0000');      
      break;
    default:
  }
if (target.text != undefined) {
  // updateLayout(cladeRoot,0,0);
  if (modeSite) {
    updateNodeColorsOnClick(selectedSite);
  }
  else {
    updateNodeColors();
  }
}
});


// Gestion de la taille de la fenetre
$('#stat_precision').on('change', function(event) {
    precisionWindowStat = parseFloat($('#stat_precision').val());
//    console.log(" DEBUG PrecisionWindowStat "+ precisionWindowStat);
//    updateLayout(cladeRoot);
    if (modeSite) {
        updateNodeColorsOnClick(selectedSite);
    }
    else {
        updateNodeColors();
    }
});


$('#action').on('click', function(event, dropdownData) {
  event.preventDefault();
  var target = event.toElement || event.relatedTarget || event.target || function () { throw "Failed to attach an event target!"; }
  myStorage.setItem("action",target.text);
  $('#button_action').text(target.text);
  switch (target.text) {
    case "View Branch Results":
      action = svgEvents.showBranchResults;
      $('#branch-res').css('color','#FF0000');
      $('#switch').css('color','#000000');
      $('#collapse').css('color','#000000');
      $('#subtree').css('color','#000000');
      break;
    case "Switch children":
      action = svgEvents.switchChildren;
      $('#branch-res').css('color','#000000');
      $('#switch').css('color','#FF0000');
      $('#collapse').css('color','#000000');
      $('#subtree').css('color','#000000');
      break;
    case "Collapse/Expand":
      action = svgEvents.collapse;
      $('#branch-res').css('color','#000000');
      $('#switch').css('color','#000000');
      $('#collapse').css('color','#FF0000');
      $('#subtree').css('color','#000000');
      break;
    case "Sub Tree/Upper Tree":
      action = svgEvents.focus;
      $('#branch-res').css('color','#000000');
      $('#switch').css('color','#000000');
      $('#collapse').css('color','#000000');
      $('#subtree').css('color','#FF0000');
      break;
    default:
  }
if (target.text != undefined) {
  var node_to_change = g.selectAll(".node");
  node_to_change.on("click", action);
}
});
$('#globalexpand').on('click', function(event, dropdownData) {
  treeRoot = d3.hierarchy(cladeRoot, function(d) {
    return d.clade;
  });
  expandTree(treeRoot);
  updateLayout(cladeRoot);
});
$('#resettree').on('click', function(event, dropdownData) {
  cladeRoot = recTree.phyloxml.phylogeny.clade;
  treeRoot = treeRootIni;
  resetTree(treeRoot);  
  updateLayout(cladeRoot);
});
$('#moreWidth').on('click', function(event, dropdownData) {
  nodeWidth += modStepX;
  myStorage.setItem("width",nodeWidth);
  updateLayout(cladeRoot);
});
$('#lessWidth').on('click', function(event, dropdownData) {
  nodeWidth -= modStepX;
  if  (nodeWidth < modStepX ) {
    nodeWidth = modStepX;
  }
  myStorage.setItem("width",nodeWidth);
  updateLayout(cladeRoot);
});
$('#loga').on('click', function(event, dropdownData) {
  if (logBranchLength) {
    logBranchLength = false;
    var max_tree_length =  getmaxlength(treeRoot,0);
    if (max_tree_length > 0) {
      nodeWidth = Math.floor(80 / max_tree_length) * 3.0;
    }
  }
  else {
    logBranchLength = true;
    var max_tree_length =  getmaxlength(treeRoot,0);
    if (max_tree_length > 0) {
      nodeWidth = Math.floor(80 / max_tree_length) * 3.0;
    }
  }
  updateLayout(cladeRoot);
});
$('#moreHeigth').on('click', function(event, dropdownData) {
  $('#lessHeigth').css('cursor','pointer')
  nodeHeight += modStepY;
  hLetter = nodeHeight - 2 ;
  if (nodeHeight > 10) {
    hUnderline = 12 ;
  }
  if (nodeHeight > 15) {
      hLetter = 15 ;
      hUnderline = 15 ;
      hStepBase = hStepNormal;
      fitPoliceSize = fitNormalPoliceSize;
      alignmentLetterSpacing = alignmentLetterNormalSpacing;
    } 
  //console.log("node = "+nodeHeight+"; letter = "+hLetter);
  updateLayout(cladeRoot);
});
$('#lessHeigth').on('click', function(event, dropdownData) {
  if (nodeHeight - modStepY >= 10 ) {
  // if (true) {
    nodeHeight -= modStepY;
    if (nodeHeight <= 15) {
      hLetter = nodeHeight - 2 ;
      hUnderline = 12 ;
      if (nodeHeight <= 10) {
        hUnderline = 9 ;
        $('#lessHeigth').css('cursor','not-allowed')
      }
      hStepBase = hStepSmall;
      alignmentLetterSpacing = alignmentLetterSmallSpacing;
      fitPoliceSize = fitSmallPoliceSize;
    }    
      //console.log("node = "+nodeHeight+"; letter = "+hLetter);
    updateLayout(cladeRoot);
  }
});
$('#redondance').on('click', function(event, dropdownData) {
    flag_redond = - flag_redond;
    //console.log("REDOND "+flag_redond);
    updateLayout(cladeRoot,2);
});
$('#histogram').on('click', function(event, dropdownData) {
    flag_histogram = - flag_histogram;
    //console.log("HISTO "+flag_histogram);
    updateLayout(cladeRoot);
});
$('#export').on('click', function(event, dropdownData) {
  var target = event.toElement || event.relatedTarget || event.target || function () { throw "Failed to attach an event target!"; }
  myStorage.setItem("export",target.text);
  switch (target.text) {
    case "Tree":
      $(".internal").css('font-size', '0.8rem');
      saveSVGTree();
      break;
    case "Alignment":
      saveSVGAlignment();
      // Correction des styles de la page après l'export
      $('.seqblock text').css({
        'font-family':'Courier',
      });
      $('.seqblock g text').css({
        'letter-spacing': (alignmentLetterSpacing)+'px',
        'font-size': fitPoliceSize+'px',
      });
      break;
    default:
  }
});
// The page has loaded
$(document).ready(function() {
  $('#upperSelected').change(function() {
    upperThresholdMode = true;
    updateLayout(cladeRoot);
  });
  $('#lowerSelected').change(function() {
    upperThresholdMode = false;
    updateLayout(cladeRoot);
  });
  if (!isCodon) {
      $('#alignment-type').hide();
  }
  $('#left-boundary, #right-boundary').prop('disabled', true);
  
  $('input#automatic-boundaries').change(function (e) {
    if ($(this).is(':checked')) {
      automaticBoundaries = true;
    } else {
      automaticBoundaries = false;
    }
  });
  var tooltips = document.querySelectorAll('.tooltip');
  window.onmousemove = function (e) {
      var x = (e.clientX + 16) + 'px';
      var y = (e.clientY + 16 + window.scrollY) + 'px';
      for (var i = 0; i < tooltips.length; i++) {
          tooltips[i].style.top = y;
          tooltips[i].style.left = x;
      }
  };
  $('button#reload-global-results').click(reloadGlobalResults);
  switch (displaySeqType) {
    case 'AA':
      $('a#aminoacid-display').css('color', '#FF0000');
      $('a#codons-display').css('color', '#000000');
      break;
    case 'Nuc':
      $('a#aminoacid-display').css('color', '#000000');
      $('a#codons-display').css('color', '#FF0000');
      break;
    default:
      $('a#aminoacid-display').css('color', '#000000');
      $('a#codons-display').css('color', '#000000');
      break;
  }
  $('main').css('margin-top', $('header').outerHeight());
  // Gestion du redimensionnement de la fenêtre
  $(window).on('ready resize', function() {
    $('main').css('margin-top', $('header').outerHeight());
  });
  // Gestion du bouton pour le type d'alignement
  $('#alignment-type').on('click', function(event, dropdownData) {
    var target = event.toElement || event.relatedTarget || event.target || function () { throw "Failed to attach an event target!"; }
    myStorage.setItem("alignment-type",target.text);
    switch (target.text) {
      case "Codons": // Codons to Amino Acids
        if (displaySeqType != 'Nuc') {
          displaySeqType = "Nuc";
          $(".aa-sequence").css("display", "none");
          $(".dna-sequence").css("display", "block");
          $("#toggle-seqtype").html("Current display: "+displaySeqType);
          $("#codons-display").css('color', '#FF0000');
          $("#aminoacid-display").css('color', '#000000');
          $('#redondance').show();

          updateLayout(cladeRoot);
          alignmentWindow.scrollLeft(alignmentWindow.scrollLeft()*3);
          selposWindow.scrollLeft(selposWindow.scrollLeft()*3);
        }
        break;
      case "Amino Acids": // Amino Acids to Codons
        if (displaySeqType != 'AA') {
          displaySeqType = "AA";
          $(".dna-sequence").css("display", "none");
          $(".aa-sequence").css("display", "block");
          $("#codons-display").css('color', '#000000');
          $("#aminoacid-display").css('color', '#FF0000');
          $('#redondance').hide();
          alignmentWindow.scrollLeft(alignmentWindow.scrollLeft()/3);
          selposWindow.scrollLeft(selposWindow.scrollLeft()/3);
          updateLayout(cladeRoot);
        }
        break;
      default:
    }
  });

  // Gestion de la modification des seuils
  $('#ps-threshold-max-down, #ps-threshold-min-down, #ps-threshold-min-up, #ps-threshold-max-up').on('change', function(event) {
    // On a modifié un des seuils
    // psThresholdsChanged = true;
    // Correction des valeurs en présence d'incohérences
    var psMaxUp = parseFloat($('#ps-threshold-max-up').val());
    var psMinUp = parseFloat($('#ps-threshold-min-up').val());
    var psMaxDown = parseFloat($('#ps-threshold-max-down').val());
    var psMinDown = parseFloat($('#ps-threshold-min-down').val());
    if (psMaxUp < psMinUp) {
      if (event.currentTarget.id == 'ps-threshold-high') {
        // On a descendu le seuil supérieur plus bas que le seuil inférieur
        // ex : psHighVal = 0.95 | psLowVal = 1.0
        // On veut donc redecsendre le seuil inférieur en conséquence
        $('#ps-threshold-low').val(psHighVal);
        psMaxUp = psMinDown;
      } else {
        // On a monté le seuil inférieur plus haut que le seuil supérieur
        // ex : psHighVal = 1.0 | psLowVal = 1.05
        // On veut donc remonter le seuil supérieur en conséquence
        $('#ps-threshold-high').val(psLowVal);
        psMinUp = psMaxDown;
      }
    }
    // Réaffectation des variables globales
    psThresholdMaxUp = psMaxUp;
    psThresholdMaxDown = psMaxDown;
    psThresholdMinUp = psMinUp;
    psThresholdMinDown = psMinDown;
//    console.log(psThresholdMaxUp);
//    console.log(psThresholdMaxDown);
//    console.log(psThresholdMinUp);
//    console.log(psThresholdMinDown);
    updateLayout(cladeRoot);
  });
})