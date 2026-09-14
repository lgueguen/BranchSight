// Fonction de mise a jour de l'affichage
// ======================================
function updateLayout(cRoot,firstLoad=0,scrolling=1) {
  var treeRoot = d3.hierarchy(cRoot, function(d) {
    return d.clade;
  });
  treeRoot.each(function (d) {
    var eventsRec = d.data.eventsRec;
    // Test si on a un noeud enrichi avec l'info eventsRec
    if (eventsRec){
      d.data.lastEvent = eventsRec[eventsRec.length - 1];
    }
    else {
      // Lors de la premiere lecture de l'arbre phyloxml, il n'y a pas de data.lastEvent defini aux noeuds
      // (Sauf les noeuds ecrases par la monophylie)
      // Mais ensuite ils sont definis : on choisit "speciation" meme si il ne s'agit pas d'une speciation en
      // realité, mais cela permet d'utiliser les codes graphiques de speciation. eventType peut etre "speciation".
      // Par contre il ne doit plus etre   "collapsed" car c'est lie a la representation dans le DOM. Utiliser nodeinfo.status
      if (d.children) {
        if (!d.data.lastEvent) {
          d.data.lastEvent =  {eventType : "speciation"};
        }
      }
    }
  });
  layout.nodeSize([nodeWidth, nodeHeight]); // *i*
  layout(treeRoot);
  // Traitement des distances dans le cas d'un arbre de gene de type phyloxml

      if (firstLoad == 1) {
        var max_tree_length =  getmaxlength(treeRoot,0);
        //console.log("MAX TREE LENGTH SIMPLE "+max_tree_length);
        //console.log("MAX TREE nodeWidth 1 "+nodeWidth);
        if (max_tree_length > 0) {
          nodeWidth = Math.floor(nodeWidth / max_tree_length) * 3.0;
        }
        //console.log("MAX TREE nodeWidth 2"+nodeWidth);
    };
    phylogeny(treeRoot,nodeWidth);
 
  // Remplissage du dico nom de sequence => tableau de couleurs
  if ((firstLoad == 1) ||  (firstLoad == 2)) {
    var nodes_ini = treeRoot.descendants()
    // var leaves_ini = nodes_ini.filter(function (e) {
    //   return !e.children;
    // });
    // On regarde tous les noeuds, pas que les feuilles
    var leaves_ini = nodes_ini.filter(function (e) {
      return true;
    });
    var leaves_sequences = leaves_ini.filter(function(d) {return (d.data.lastEvent.aaAlign);});
    leaves_sequences.forEach (function (d,i) {
      var leaf_name = d.data.name;
      var col_aas = [];
      var col_dnas = [];
      //console.log("SEQUENCE NAME = " + leaf_name);
      var dna_sequence = d.data.lastEvent.dnaAlign;
      var aa_sequence = d.data.lastEvent.aaAlign;
      var aas = aa_sequence.split("");
      //console.log("AA "+aa_sequence)
      aas.forEach (function (aa,j) {
        var col = hslaa[aa];
        if (col == undefined) {
          col  = "grey"
        }
        else {
          col = hsl2col(col);
        };
        col_aas.push([col,aa]);
      })
      dico_aa[leaf_name] = col_aas;
      if (dna_sequence != undefined) {
        var dnas = dna_sequence.split("");
        for (i = 0; i <= dnas.length-3; i += 3) {
          var codon = dnas[i]+dnas[i+1]+dnas[i+2];
          var col = "grey";
          if (flag_redond == -1) {
            col = colordna_red[codon];
              //console.log("redundant color "+ col);
          } else {
             col = colordna[codon];
              //console.log("classical color "+ col);
          }
          if (col == undefined) {col = "grey"};
          col_dnas.push([col,codon]);
        }
      }
      dico_dna[leaf_name] = col_dnas;
    })
  }
  updateSvg(treeRoot,firstLoad);
  //console.log("SCROLLING "+scrolling)
  if (scrolling == 1) {
  updateScroll();
  }
}

// Fonction principale : mise a jour du svg
// ----------------------------------------
function updateSvg(treeRoot,firstLoad ,config = {}) {
  //console.log('Start updatesvg');
  // Suppression de l'ancien affichage de l'alignement et du graphe
  d3.selectAll('#ps-graph, #sequences-container').remove();
  // ====================================================
  // | Début modif : variables pour graphe et séquences |
  // ====================================================
  // Ajustement du pas selon le type de séquences affichées
  hStep = displaySeqType == 'AA' ? hStepBase : hStepBase*3 ;
  var seqblocks = d3.selectAll('.seqblock');
  if (!resultsJSON) {
    showTooltip();
    writeTooltip('Please click any tree node to display branch data');
    resultsJSON = JSON.parse(treeRoot.data.branch_info.results);
  } 
  var psGraphMargins = {top: 20, right: 50, bottom: 50, left: 50+seq_lg};
  var psGraphWidth = resultsJSON.length*hStep;
  var psGraphHeight = 100;
  // ==============================================================
  // Fin modif : variables pour graphe et séquences
  // ==============================================================
  var scrollref=0;
  var configLayout = {
    layout : config.layout || "cladogramSpecial",
    links : config.links || "shoulder",
    symbolSize : config.symbolSize || 128,
    lengthLinkLoss : config.lengthLinkLoss || 10,
    linkStrokeSize : config.linkStrokeSize || 3,
    nodeWidth : config.nodeWidth || 30,
    nodeHeight : config.nodeHeight || 30,
    margin : config.margin || { top: 10, down: 20, left: 500  , right: 50},
    color : config.color || {
      speciation : "#2F4F4F",
      leaf : "#FF7F0E"
    },
    symbols : config.symbols || {
      speciation : "symbolCircle",
      collapsed : "symbolCircle"
    }
  }
  var nodes = treeRoot.descendants();
  var links = treeRoot.links();
  // Dimensions du svg
  var minX = d3.min(nodes, function(d) {
    return d.x;
  });
  var maxX = d3.max(nodes, function(d) {
    return d.x;
  });
  var minY = d3.min(nodes, function(d) {
    return d.y;
  });
  var maxY = d3.max(nodes, function(d) {
    return d.y;
  });
  var widthSVG = maxX - minX;
  var heightSVG = maxY - minY;
  var ajustementX = resultsJSON.length*hStep - 120 ;
  var ajustementY = -200 ;
  var decalageY = -150 ;
  svg1.attr("width", widthSVG + margin.right + margin.left);
  svg1.attr("height", heightSVG + margin.top + margin.down+300      + ajustementY);
  g.attr("transform", "translate(" + (margin.right - minX - 250) + "," + (margin.top - minY + 50 + psGraphHeight + decalageY) + ")");
  //  D3: Objets LINK
  var link = g.selectAll(".link").data(links);  // associe l'element link a la donnee links
  //EXIT gere les elements qui ont disparu
  link.exit().remove();
  //ENTER gere les elements nouveaux
  var linkEnter =
    link
    .enter()
    .append("path")
    .attr("class", "link");    
  linkEnter
    .merge(link)  // fusionne les nouveaux elements avec les anciens
    .attr("fill","none")
    .attr("stroke-width",configLayout.linkStrokeSize)
    .attr("stroke","#0e2e2e")
    .attr("d", diagonal)
            .on('mouseover', function(d) {
          if (selectedNode) {
            showTooltip();
            // Display a tooltip when hovering the cursor over a node
            if (d.target.data.branch_info) {
                // if (d.target.data.positiveRate) {
                //   writeTooltip(d.target.data.positiveRate.toFixed(3));
                // } else {
                if (modeSite) {
                    writeTooltip(getTransitionJSON(d.source.data.lastEvent, d.target.data.lastEvent, selectedSite));
                } else {
                    writeTooltip(getPositiveRateJSON(d.target.data.branch_info.results, leftBoundary, rightBoundary).toFixed(3));
                }
//                }
            } else {
                writeTooltip('');
            }
          }
        })
        .on('mouseout', function() {
          if (selectedNode) {
            hideTooltip();
          }
        });
  // D3:  Objets NODE
  var node = g.selectAll(".node").data(nodes);
  node.exit().remove();
  var nodeEnter =
    node
    .enter()
    .append("g")
    .attr("class", "node");
  nodeEnter
    .append("g")
      .attr("class", "gsymbol")
      .append("path")
        .attr("class", "symbol")
        .style('stroke-width', '0')
        .style('stroke', '#000000')
        .on('mouseover', function(d) {
          if (selectedNode) {
            showTooltip();
            // Display a tooltip when hovering the cursor over a node
            if (d.data.branch_info) {
              if (modeSite) {
                  writeTooltip(getValueJSON(d.data.branch_info.results, selectedSite).toFixed(3));
              } else {
                  writeTooltip(getPositiveRateJSON(d.data.branch_info.results, leftBoundary, rightBoundary).toFixed(3));
              }    
            } else {
              writeTooltip('');
            }
          }
        })
        .on('mouseout', function() {
          if (selectedNode) {
            hideTooltip();
          }
        })
        .on('click', function(d) {
          if (action == svgEvents.showBranchResults && d.data.branch_info) {
            resultsJSON = JSON.parse(d.data.branch_info.results);
            //console.log("DEBUG")
            //console.log(this)
            selectedNode = d3.select(this);
          }
        });
  if (firstLoad) {
    d3.selectAll('path.symbol')
      .style('stroke-width', '0')
      .style('stroke', '#888888')
      .style('fill', '#000000')
      ;
  } else {
    d3.selectAll('path.symbol')
      .style('stroke-width', '0')
      .style('stroke', '#000000')
      .style('fill', '#000000')
      ;
    if (selectedNode) {
      // selected node color and outline
      selectedNode
        .style('stroke-width', '2')
        .style('stroke', '#000000')
        .style('fill', '#FFFFFF')
        ;
    }
  }
  nodeEnter
    .append("text")
    .attr("class", "label");
  nodeEnter
    .append("text")
    .attr("class", "internal");
  var allNodes =
    nodeEnter
    .merge(node)
    .attr("transform", function(d) {
        return "translate(" + [d.x, d.y] + ")";
    })
    .on("click",action);  
  allNodes
    .select(".gsymbol")
    .attr("transform", function(d) {
      var str = "";
      if (d.data.lastEvent && d.data.lastEvent.eventType === "leaf")  {
        str += "rotate(90)";
      } else {
        str += "";
      }
      return str;
    });    
  allNodes
    .select(".symbol")
    .attr("d", function(d) { //Type de symbol
       // Verification si le noeud est collapse
      if (d.data.nodeinfo && d.data.lastEvent) {
          if (d.data.nodeinfo.status === "collapsed") {
            return symbol.type(d3.symbolCircle)()
          }
      }
      if (d.data.lastEvent) {
        switch (d.data.lastEvent.eventType) {
          case "monophyletic":          // cas ou le noeud monophyletic n'est pas pas collapse
            return symbol.type(d3[configLayout.symbols.speciation])();
            break;
          case "speciation":
            return symbol.type(d3[configLayout.symbols.speciation])();
            break;
          case "leaf":
            return symbol.type(d3.symbolTriangle)();
            break;
          default:
            return symbol.type(d3.symbolCross)();
            break;
        }
      } else {
        return symbol.type(d3.symbolCross)();
      }
    });    
  allNodes
    .select(".label")
    .transition()
    .attr("y", function(d) {
      // return d.children ? -8 : 3;
      return d.children ? -8 : 4;
    })
    .attr("x", function(d) {
      // return d.children ? -8 : 8;
      return d.children ? -8 : 12;
    })
    .style("text-anchor", function(d) {
      return d.children ? "end" : "start";
    })
    .text(function(d) {
      var name = "";
      if (d.data.nodeinfo) {
        if (d.data.nodeinfo.status === "collapsed") {
            name += "("+d.data.nbseq.nbSeq + " sequences) " + d.data.name ;
        }
      }
      if (d.data.name && options.includes("name") && d.data.lastEvent.eventType === "leaf") {
          name += d.data.name ;
      }
      return name;
    })
    .style("fill",function (d) {
    })
    .style("font-family","sans-serif")
    .style("font-size",12+"px");
  allNodes
    .select(".internal")
    .transition()
    .attr("y", 11)
    .attr("x", -10)
    .style("text-anchor","end")
    .text(function(d) {
      var name = "";
      if (options.includes("intname")) {
        if (d.data.confidence) {
          // Display confidence below branches
          var name = d.data.confidence._text;
        }
      }
      return name
    });
  // Get tree characteristics (do this only once)
  if (firstLoad) {
    nbBranches = g.selectAll('path.symbol').size();
    g.selectAll('path.symbol').each(function(elem) {
      if (elem.data.lastEvent.aaAlign) {
        nbLeaves ++;
        if (elem.data.lastEvent.aaAlign.length > nbSites) {
          nbSites = elem.data.lastEvent.aaAlign.length;
        }
      }
    });
  }
  // =======================
  // Selection des feuilles
  // =======================
  var leaves = nodes.filter(function (e) {
    return !e.children;
  });
  var internal_nodes = nodes.filter(function (e) {
    flag = !e.children;
    return !flag;
  }); 
  //console.log("internal "+internal_nodes) 
  var widthSVG0 = resultsJSON.length;
  var heightSVG0 = 60;
  var psGraphGlobalWidth = resultsJSON.length;
  var psGraphGlobalHeight = 60;
  var psGraphGlobalMargins = {top:0, right:0, bottom:0, left:0};
  var svg0 = d3.select('#svg0');
  svg0.attr('width', widthSVG0 + 15);
      // .attr("height",psGraphGlobalHeight);
      // .attr('transform', 'translate('+ window.innerWidth*0.5+ ')');
  var svg3 = d3.select('#svg3');
  var svgbidon = d3.select('#svgbidon');
  svg3.attr("width",  50 + seq_lg + 10 + resultsJSON.length*hStep + 10 + seq_lg + 50);
  if (flag_histogram == -1 )  {
    svg3.attr("height",  0);
    svgbidon.attr("height",  0);
    svg0.attr("height",  0);
    d3.select("#tree-container")
    .style("top", "0px");
    d3.select("#alignment-container")
    .style("top", "0px");    
  }
  else {
    svg3.attr("height",  150);
    svgbidon.attr("height",  150);
    svg0.attr("height",  psGraphGlobalHeight);
    d3.select("#tree-container")
    .style("top", "300px");
    d3.select("#alignment-container")
    .style("top", "300px");      
  }
  d3.select('#ps-graph-global').remove();
  if (!document.getElementById('ps-graph-global')) {
    // Ajout du graphe des résultats au SVG
    svg0.append('g')
      .attr('id', 'ps-graph-global')
      .attr('transform', 'translate('+psGraphGlobalMargins.left+', '+psGraphGlobalMargins.top+')');
  }
  var psGraphGlobal = d3.select('#ps-graph-global');
  // Calcul du décalage des sites pour le positionnement des rectangles du graphe
  // Association des données aux futurs rectangles
  var globalRects = psGraphGlobal
    .selectAll('.ps-rects-global')
    .data(function() {
      return resultsJSON;
    });
  globalRects.exit().remove();
  // Affichage des rectangles
  var i = 0;
  globalRects.enter()
    .append('rect')
      .attr('class', 'ps-rects-global')
      .attr('width', function(d) { return 1; })
      .attr('height', function(d) { return psGraphGlobalHeight; })
      .attr('x', function(d) { i ++; return i - siteOffset; })
      .attr('y', function(d) { return 0; })
      .style('fill', function(d) {
        return getColorFromResult(d);
      });          
  // Handle click-drag events on the small graph
  document.getElementById(compactWindow.attr('id')).onmousedown = function(e) {
    scrollAlignment(e);
    compactWindow.on('mousemove', function (e) {
      modeSite = false;
      //console.log("MODE_SITE " + modeSite);
      scrollAlignment(e);
    });
    compactWindow.on('mouseup', function (e) {
      // scrollAlignment(e);
      compactWindow.off('mousemove mouseup');
    });
  };
  document.getElementById(alignmentWindow.attr('id')).onscroll = function() {
    // Generic scrolling (even without the small graph)
    // Draw rectangle at currently visualized site
    $('.navigation-window').remove();
    var siteNumber = - Math.trunc((60 + seq_lg - alignmentWindow.scrollLeft() - alignmentWindow.innerWidth()/2)/hStep) - 2;
    var navWindowWidth = alignmentWindow.innerWidth()/hStep;
    var siteStart = Math.trunc(siteNumber - navWindowWidth/2) + 5;
    var siteEnd = Math.trunc(siteNumber + navWindowWidth/2) + 5;
    var navWindowRect = psGraphGlobal
      .append('rect')
        .attr('class', 'navigation-window')
        .attr('x', siteStart)
        .attr('y', 0)
        .attr('width', navWindowWidth)
        .attr('height', psGraphGlobalHeight)
        .style('fill', 'rgba(0,0,0,0.17)')
        .style('stroke', navWindowColor)
        .style('stroke-width', '0');
    // autre scroll
    var ss =  alignmentWindow.scrollLeft();
    selposWindow.scrollLeft(ss, 0);
    if (automaticBoundaries) {
      leftBoundary = siteStart;
      rightBoundary = siteEnd;
    }
    updateNodeColors();
  };
  document.getElementById(selposWindow.attr('id')).onscroll = function() {
    // Generic scrolling (even without the small graph)
    // Draw rectangle at currently visualized site
    $('.navigation-window').remove();
    var siteNumber = - Math.trunc((60 + seq_lg - selposWindow.scrollLeft() - selposWindow.innerWidth()/2)/hStep) - 2;
    var navWindowWidth = selposWindow.innerWidth()/hStep;
    var siteStart = Math.trunc(siteNumber - navWindowWidth/2) + 5;
    var siteEnd = Math.trunc(siteNumber + navWindowWidth/2) + 5;
    var navWindowRect = psGraphGlobal
      .append('rect')
        .attr('class', 'navigation-window')
        .attr('x', siteStart)
        .attr('y', 0)
        .attr('width', navWindowWidth)
        .attr('height', psGraphGlobalHeight)
        .style('fill', 'rgba(0,0,0,0.17)')
        .style('stroke', navWindowColor)
        .style('stroke-width', '0');
    var ss =  selposWindow.scrollLeft();
    alignmentWindow.scrollLeft(ss, 0);
    if (automaticBoundaries) {
      leftBoundary = siteStart;
      rightBoundary = siteEnd;
    }
    updateNodeColors();
  };
  console.log("mise en place des graphes");
  // ======================================================
  // | Début : graphe des résultats de sélection positive |
  // ======================================================
  var  svg2 = d3.select("#svg2");
  svg2
    .attr("width",  50 + seq_lg + 10 + resultsJSON.length*hStep + 10 + seq_lg + 50)
    .attr("height",  heightSVG + margin.top + margin.down + 100);
  // === Définition de l'échelle pour l'axe x ===
  var xScale = d3.scaleLinear()
    .domain([1, resultsJSON.length+1])
    .range([1, psGraphWidth]);
  // === Définition de l'échelle pour l'axe y ===
  var yScale = d3.scaleLinear()
    .domain([minStatJSON(resultsJSON), maxStatJSON(resultsJSON)])
    .range([psGraphHeight, 0]);
  // === Ajout du graphe des résultats au SVG ===
  var psGraph = svg3.append('g')
    .attr('id', 'ps-graph')
    .attr('transform', 'translate('+(50 + seq_lg)+', '+psGraphMargins.top+')');
  // === Ajout de l'axe x au graphe ===
  psGraph.append('g')
    .attr('transform', 'translate(10,' + psGraphHeight + ')')
    .call(d3.axisBottom(xScale).ticks(resultsJSON.length/5)); // graduations tous les 5 codons / AA
  // === Ajout de l'axe y au graphe ===
  psGraph.append('g')
    .call(d3.axisLeft(yScale));
  // === Calcul du décalage des sites pour le positionnement des rectangles du graphe ===
  // === Association des données aux futurs rectangles ===
  var rects = psGraph
    .selectAll('.ps-rects')
    .data(function() {
      // console.log('big graph: attaching resultsJSON', resultsJSON);
      return resultsJSON;
    });
  var rectsfond = psGraph
  .selectAll('.ps-rects')
  .data(function() {
    // console.log('big graph: attaching resultsJSON', resultsJSON);
    return resultsJSON;
  });   
  
  // === Calcul du décalage des sites pour le positionnement des rectangles du graphe ===
  // === Association des données aux futurs rectangles ===
  // === Affichage des rectangles de fond ===
  var i = 0;
  rectsfond.enter()
    .append('rect')
      // .attr('id', function(d) {
      //   return 'ps-rect-'+d.site;
      // })
      .attr('class', 'ps-rects')
      // .attr('title', function(d) {
      //   return d.site;
      // })
      .attr('width', function(d) { return hStep; })
      .attr('height', 100)
      .attr('x', function(d) { i ++; return xScale(i - siteOffset); })
      .attr('y', 0)
      .style('fill', "white")
      .style('stroke', '#000000')
      .style('stroke-width', '0')
      .on('mouseover', function(d) {
        // console.log(d);
        writeTooltip(d.toFixed(3));
        showTooltip();
      })
      .on('mouseout', hideTooltip)
      .on('click', function(index,d) {
        //console.log("click " + index + "=>"+d);
        modeSite = true;
        selectedSite = d;
        writeTooltip(d+1);
        showTooltip();
        updateNodeColorsOnClick(d);
      })
      ;  

  // === Affichage des rectangles indiquant la valeur===
  var i = 0;
  rects.enter()
    .append('rect')
      // .attr('id', function(d) {
      //   return 'ps-rect-'+d.site;
      // })
      .attr('class', 'ps-rects')
      // .attr('title', function(d) {
      //   return d.site;
      // })
      .attr('width', function(d) { return hStep; })
      .attr('height', function(d) {
      if ( yScale(0)>=yScale(d) ) {
      return yScale(0) - yScale(d);
      }
      else {
      return yScale(d) - yScale(0);
      }
      })
      .attr('x', function(d) { i ++; return xScale(i - siteOffset); })
      .attr('y', function(d) { 
        if ( yScale(0)<=yScale(d) ) {
        return yScale(0);
        }
        else {
          return yScale(d);
        }
      })
      .style('fill', function(d) {
        return getColorFromResult(d);
      })
      .style('stroke', '#000000')
      .style('stroke-width', '0')
      .on('mouseover', function(d) {
        // console.log(d);
        writeTooltip(d.toFixed(3));
        showTooltip();
      })
      .on('mouseout', hideTooltip)
      .on('click', function(index,d) {
        //console.log("click " + index + "=>"+d);
        modeSite = true;
        selectedSite = d;
        writeTooltip(d+1);
        showTooltip();
        updateNodeColorsOnClick(d);
      })
      ;



  // Ajout du seuil supérieur au graphe
  if (upperThresholdMode){
    psGraph.append('line')
      .attr('x1', xScale(1))
      .attr('y1', yScale(psThresholdMaxUp))
      .attr('x2', xScale(resultsJSON.length+1))
      .attr('y2', yScale(psThresholdMaxUp))
      .style('stroke', '#000000')
      // .style('stroke', psThresholdMaxUpBgColor)
      .style('stroke-width', '2')
      .style('stroke-dasharray', ('4, 4'));
    // Ajout du seuil inférieur au graphe
    psGraph.append('line')
      .attr('x1', xScale(1))
      .attr('y1', yScale(psThresholdMinUp))
      .attr('x2', xScale(resultsJSON.length+1))
      .attr('y2', yScale(psThresholdMinUp))
      .style('stroke', '#000000')
      // .style('stroke', psThresholdMinUpBgColor)
      .style('stroke-width', '1')
      .style('stroke-dasharray', ('4, 4'));
  } else {
    psGraph.append('line')
      .attr('x1', xScale(1))
      .attr('y1', yScale(psThresholdMaxDown))
      .attr('x2', xScale(resultsJSON.length+1))
      .attr('y2', yScale(psThresholdMaxDown))
      .style('stroke', '#000000')
      // .style('stroke', psThresholdMaxDownBgColor)
      .style('stroke-width', '2')
      .style('stroke-dasharray', ('4, 4'));
    // Ajout du seuil inférieur au graphe
    psGraph.append('line')
      .attr('x1', xScale(1))
      .attr('y1', yScale(psThresholdMinDown))
      .attr('x2', xScale(resultsJSON.length+1))
      .attr('y2', yScale(psThresholdMinDown))
      .style('stroke', '#000000')
      // .style('stroke', psThresholdMinDownBgColor)
      .style('stroke-width', '1')
      .style('stroke-dasharray', ('4, 4'));
  }
  psGraph.append('br');


  console.log("graphes OK");
  // ====================================================
  // | Fin : graphe des résultats de sélection positive |
  // ====================================================
  // ==============================
  // | Début : blocs de séquences |
  // ==============================
  seq_dataset = [];
  // var leaves_sequences = leaves.filter(function(d) {return (isNuc ? d.data.lastEvent.dnaAlign : d.data.lastEvent.aaAlign);});
  var leaves_sequences = leaves.filter(function(d) {return (d.data.lastEvent.aaAlign);});
  var leaves_sequences = leaves.filter(function(d) {return (true);});
  // console.log("debug leave "+ d.data.name)
    leaves_sequences.forEach (function (d,i) {
      var leaf_name = d.data.lastEvent.name;
      var dna_sequence = d.data.lastEvent.dnaAlign;
      var aa_sequence = d.data.lastEvent.aaAlign;
      var lengthmax=d.data.lastEvent.lengthMax;
      var locnode = {
        name:d.data.name,
        y:d.y,
        dna:d.data.lastEvent.dnaAlign,
        aa:d.data.lastEvent.aaAlign,
      };
      seq_dataset.push(locnode);
    });
    
  var internal_sequences = internal_nodes.filter(function(d) {return (true);});
  internal_nodes.forEach (function (d,i) {
    //console.log("debug inter "+ d.data.name)
    //console.log("debug inter "+ d.data.nbseqspe.nbSeq)
    // if d.data.nodeinfo === "collapsed"
    if (d.data.nodeinfo) {
      //console.log("debug inter "+ d.data.name)
    }
    var leaf_name = d.data.lastEvent.name;
    var dna_sequence = d.data.lastEvent.dnaAlign;
    var aa_sequence = d.data.lastEvent.aaAlign;
    var lengthmax=d.data.lastEvent.lengthMax;
    var locnode = {
      name:d.data.name,
      y:d.y,
      dna:d.data.lastEvent.dnaAlign,
      aa:d.data.lastEvent.aaAlign,
    };
    // seq_dataset.push(locnode);
    
  });  
  var seqDiv = d3.select('#svg2')
    .append('g')
    .attr('id', 'sequences-container')
    //.attr("transform", "translate(" + (margin.right - minX - 300) + "," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")")
    .attr("transform", function(d) {
      //return "translate(" + (margin.right - minX - 300) + "," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")"
      return "translate( 0 ," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")"
    })

  // Association des données aux futurs blocs de séquences
  var seqblocks = seqDiv
    .selectAll('.seqblock')
    .data(seq_dataset);
  // Création et positionnement des blocs de séquences
  var seqblocksEnter = seqblocks.enter().append('g')
    .attr('class', 'node seqblock')
    .attr('transform', function(d) {
      var transx = psGraphMargins.left;
      var transy = d.y; // fonctionne seulement avec nodeHeight = 15 *amod*
      return 'translate('+transx+', '+transy+')';
    });
    seqblocksEnter.append("text")
      .attr('class', 'sequence-name')
      .text(function(d) {
        return d.name;
      })
      .attr('transform', 'translate('+'-'+(Math.floor(seq_lg))+','+5+')');
  var maxYSequences = function(seqDataset) {
    var currentmax = seqDataset[0].y;
    seqDataset.forEach(function(seq) {
      currentmax = seq.y > currentmax ? seq.y : currentmax;
    });
    return currentmax;
  }
  var maxYSeq = maxYSequences(seq_dataset);
  var minYSequences = function(seqDataset) {
    var currentmin = seqDataset[0].y;
    seqDataset.forEach(function(seq) {
      currentmin = seq.y < currentmin ? seq.y : currentmin;
    });
    return currentmin;
  }
  var minYSeq = minYSequences(seq_dataset);
  // Alignement : ajout des conteneurs au SVG
  var seqTextAA = seqblocksEnter.append('g')
      .attr('class', 'aa-sequence')
      .style('display', function() {
        return displaySeqType == 'AA' ? 'block' : 'none';
      });
  var seqTextNuc = seqblocksEnter.append('g')
      .attr('class', 'dna-sequence')
      .style('display', function() {
        return displaySeqType == 'Nuc' ? 'block' : 'none';
      });
  // Alignement : ajout des rectangles
  var aaSequence = d3.select('.aa-sequence');
  var dnaSequence = d3.select('.dna-sequence');
  var targetSequence = displaySeqType == 'AA' ? aaSequence : dnaSequence;
  // Alignement : ajout des couleurs AA
  console.log("Genere aa_dataset and dna_dataset");
  aa_dataset = [];
  dna_dataset = [];
  var opacityAA = 0;
  var opacityNuc = 0;
  if (displaySeqType == 'AA') {
      opacityAA = 0.5;
      opacityNuc = 0.0;
  }
  else {
    opacityAA = 0.0;
    opacityNuc = 0.5;
  }
  d3.select("#aa-container").remove(); // on vire  tout
  d3.select("#dna-container").remove();
  var leaves_sequences = leaves.filter(function(d) {return (d.data.lastEvent.aaAlign);});
  leaves_sequences.forEach (function (d,i) {
    var index = i;
    var seqname =  d.data.name;
    var aas = dico_aa[seqname];
    
    //console.log("AA2 "+aas)
    if (aas != undefined) {
      aas.forEach (function (aa,j) {
        var opac =  getOpacFromResult(resultsJSON[j]);
        var locnode = {
          index:i,
          name:d.data.name,
          y:d.y,
          pos:j,
          col:aa[0],
          opac:opac,
          val:aa[1],
        };
        if (displaySeqType == 'AA') {
          aa_dataset.push(locnode);
        }
        })
    }
      var dnas = dico_dna[seqname];
      if (dnas != undefined) {
        dnas.forEach (function (dna,j) {
          var opac =  getOpacFromResult(resultsJSON[j]);
          var locnode = {
          index:i,
          name:d.data.name,
          y:d.y,
          pos:j,
          col:dna[0],
          opac:opac,
          val:dna[1],
        };
        if (displaySeqType == 'Nuc') {
          dna_dataset.push(locnode);
        }
        })
    }
    });
  //console.log("aa_dataset and dna_dataset OK");
  // seqTextAA.append('text')
  //   .style('letter-spacing', (alignmentLetterSpacing)+'px') // Espacement des lettres pour correpondre au graphe
  //   .style("font-size",fitPoliceSize+"px")
  //   .attr('x', function(d) {return 0;})
  //   .attr('y', function(d) {return Math.floor(hStepBase/3);})
  //   .text(function(d) {return d.aa;});
  // seqTextNuc.append('text')
  //   .style('letter-spacing', (alignmentLetterSpacing)+'px') // Espacement des lettres pour correpondre au graphe
  //   .style("font-size",fitPoliceSize+"px")
  //   .attr('x', function(d) {return 0;})
  //   .attr('y', function(d) {return Math.floor(hStepBase/3);})
  //   .text(function(d) {return d.dna;});
  console.log("sequences  OK");
  if (displaySeqType == 'AA') {
  var aaDiv = d3.select('#svg2')
      .append('g')
      .attr('id', 'aa-container')
      //.attr("transform", "translate(" + (margin.right - minX - 300) + "," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")")
      .attr("transform", "translate( 0," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")")
    var aablock =
      aaDiv.selectAll(".aa_block")
      .data(aa_dataset);
      //EXIT
      aablock.exit().remove();
      //ENTER
      var aablockEnter =
      aablock
      .enter()
      .append("g")
      .attr("class", "aa_block");
      //ENTER + updatesATE
      aablockEnter
      .merge(aablock)
      .append("rect")
      .attr('width', hStep)
      .attr('height', hLetter )
      .style("fill", function(d){ return (d.col)})
      .style("opacity",opacityAA);
      aablockEnter
      .merge(aablock)
      .append("text")
      .text(function(d) { return d.val; })
      .attr("transform","translate(0, "+(fitPoliceSize - 2) +")")
      .style("font-family","Courier")
      .style("font-size",fitPoliceSize+"px");

      //ENTER + updatesATE
      aablockEnter
      .merge(aablock)
      .attr("transform",function(d,i) {
        var transx = (d.pos)*hStep + seq_lg + 50 ;
        if (isDNA == 1 ) {
          transx = transx - 80 ;// ???
        }
        var transy = d.y - 7 ;
        return "translate("+transx+", "+transy+")";
      })
      .append("rect")
      .attr("class","aa_rect")
      .attr("transform","translate(0,"+hUnderline+")" )
      .attr('width', hStep)
      .attr('height', 2)
      .attr("stroke","black")
      .attr("stroke-opacity", 1)
      .style("fill", "black")
      .style("opacity",function(d){
        if (displaySeqType != "AA") {
          return 0
        } else {
        return d.opac
        }
      });
      aablockEnter
      .merge(aablock)
      .attr("transform",function(d,i) {
        var transx = (d.pos)*hStep + seq_lg + 50 ;
        if (isDNA == 1 ) {
          transx = transx - 80 ;// ???
        }
        var transy = d.y - 7 ;
        return "translate("+transx+", "+transy+")";
      });
    }
    if (displaySeqType == 'Nuc') {
      var dnaDiv = d3.select('#svg2')
        .append('g')
        .attr('id', 'dna-container')
       // .attr("transform", "translate(" + (margin.right - minX - 300) + "," + (margin.top - minY + 50 + psGraphHeight  + decalageY) + ")")
       .attr("transform", "translate(0 ," + (margin.top - minY + 50 + psGraphHeight  + decalageY ) + ")")
 
      var dnablock =
        dnaDiv.selectAll(".dna_block")
        .data(dna_dataset);
        //EXIT
        dnablock.exit().remove();
        //ENTER
        var dnablockEnter =
        dnablock
        .enter()
        .append("g")
        .attr("class", "dna_block");
        //ENTER + updatesATE
        dnablockEnter
        .merge(dnablock)
        .append("rect")
        .attr('width', hStep)
        .attr('height', hLetter )
        .style("fill", function(d){ return (d.col)})
        .style("opacity", opacityNuc);
        dnablockEnter
        .merge(dnablock)
        .append("text")
        .text(function(d) { return d.val; })
        .attr("transform","translate(3, "+(fitPoliceSize - 2) +")")
        .style("font-family","Courier")
        .style("font-size",fitPoliceSize+"px");
        // Rectangles qui soulignent les position sel pos
        //ENTER + updates
        dnablockEnter
        .merge(dnablock)
        .attr("transform",function(d,i) {
          var transx = (d.pos)*hStep + seq_lg + 50 ;
          var transy = d.y - 7 ;
          return "translate("+transx+", "+transy+")";
        })
        .append("rect")
        .attr("class","dna_rect")
        .attr("transform","translate(0,"+hUnderline+")" )
        .attr('width', hStep)
        .attr('height', 2 )
        .attr("stroke","black")
        .attr("stroke-opacity", 1)
        .style("fill", "black")
        .style("opacity", function(d){
          if (displaySeqType != "Nuc") {
            return 0
          } else {
            return d.opac
          }
        });
      }
      console.log("sequence blocs OK");
  // ==============================================================
  // Fin modif : blocs de séquences
  // ==============================================================
  firstLoad = 0;
    console.log('Fin updatesvg');
}


// Mise a jour de la couleur des noeuds
// ------------------------------------
function updateNodeColors() {

// Precompute positive rates
g.selectAll('path.symbol')
  .each(function(d) {
    if (d.data.branch_info) {
      d.data.positiveRate = getPositiveRateJSON(
        d.data.branch_info.results,
        leftBoundary,
        rightBoundary
      ).toFixed(3);
//      console.log("DEBUG");
//      console.log(d);

    }
  });

// The color of a node depends on the number of positive sites
g.selectAll('path.symbol')
    .style('fill', function(d) {
      if (d.data.branch_info) {       
//        console.log("DATA ")
//        console.log(d.data)
          var rate = d.data.positiveRate;
        return colorByRate(rate);
      }
    });
  // The color of a branch depends on the number of positive sites in its node
  g.selectAll('path.link')
    .style('stroke', function(d) {
      if (d.target.data.branch_info) {
            var rate = d.target.data.positiveRate;
            return colorByRate(rate);
          }
    });


  // Selected node color and outline
  if (action == svgEvents.showBranchResults && selectedNode) {
    selectedNode
      .style('stroke-width', '2')
      .style('stroke', '#000000')
      .style('fill', '#FFFFFF')
      ;
  }
}
// Mise a jour de la couleur des noeuds apres un click
// ---------------------------------------------------
function updateNodeColorsOnClick(site) {
    // precompute all rates
    g.selectAll('path.symbol')
        .each(function(d) {
            if (d.data.branch_info) {
                d.data.positiveRate = getValueJSON(d.data.branch_info.results, site);
            }
        });
        
    // The color of a node depends on the number of positive sites
    g.selectAll('path.symbol')
        .style('fill', function(d) {
            if (d.data.branch_info) {
                if (d.data.positiveRate > psThresholdMaxUp) {
                    return("#FF0000");
                } else {
                    if (d.data.positiveRate > psThresholdMinUp) {
                        return("#FFFF00");
                    } else {
                        if (d.data.positiveRate < psThresholdMaxDown) {
                            return("#0000FF");
                        } else {
                            if (d.data.positiveRate < psThresholdMinDown) {
                                return("#00FFFF");
                            } else {
                                return("#DDDDDD");
                            }
                        }
                    }
                }
            } else {
                //console.log("debug no data "+d)
            }
        });


    // The color of a branch depends on the value in its node
    g.selectAll('path.link')
        .style('stroke', function(d) {
            if (d.target.data.branch_info) {
                if (d.target.data.positiveRate > psThresholdMaxUp) {
                    return("#FF0000");
                } else {
                    if (d.target.data.positiveRate > psThresholdMinUp) {
                        return("#FFFF00");
                    } else {
                        if (d.target.data.positiveRate < psThresholdMaxDown) {
                            return("#0000FF");
                        } else {
                            if (d.target.data.positiveRate < psThresholdMinDown) {
                                return("#00FFFF");
                            } else {
                                return("#DDDDDD");
                            }
                        }
                    }
                }
            } else {
                //console.log("debug2 no data "+d)
         }
    });

  // Selected node color and outline
  if (action == svgEvents.showBranchResults && selectedNode) {
    selectedNode
      .style('stroke-width', '2')
      .style('stroke', '#000000')
      .style('fill', '#FFFFFF')
      ;
  }
}