svgEvents = {};

svgEvents.switchChildren = function(d) {
  if (d.children[0] && d.children[1]) {
    var temp = d.data.clade[0];
    d.data.clade[0] = d.data.clade[1];
    d.data.clade[1] = temp;
    updateLayout(cladeRoot);
  }
}

svgEvents.collapse = function(d) {
  if (d.data.clade) {
    d.data._clade = d.data.clade;
    d.data.clade = null;
    d.data.nodeinfo = {status : "collapsed"};
  }
  else if (d.data._clade) {
    d.data.clade = d.data._clade;
    d.data._clade = null;
    d.data.nodeinfo = {status : "extended"};
  }
  // var  vertical = d3.select("#svg2");
  // vertical.select("#vertical_rect").remove();
  // displayVerticalLine(d,heightSVG,margin.top + 50 + psGraphHeight  + decalageY);  
  updateLayout(cladeRoot);
}

svgEvents.showBranchResults = function(d) {
  if (d.data.branch_info) {
    resultsShown = d.data.branch_info.branch_results;
    updateLayout(cladeRoot);
  }
}

svgEvents.openall = function(n) {
  var fils = n.descendants();
  fils.forEach(function (d) {
    if   (!d.data.clade) {
      d.data.clade = d.data._clade;;
    }
  })
  // var  vertical = d3.select("#svg2");
  // vertical.select("#vertical_rect").remove();
  updateLayout(cladeRoot);
}

svgEvents.focus = function(d) {
  if(d.data != cladeRoot) {
    var parent = d.parent;
    _cladeRoot = parent.data;
    _Parent = parent;
    cladeRoot = d.data;
  }
  else {
    var parent = _Parent;
    cladeRoot = parent.data;
    _Parent = parent.parent;
  }
  // var  vertical = d3.select("#svg2");
  // vertical.select("#vertical_rect").remove();
  updateLayout(cladeRoot);
}

svgEvents.focus2 = function(d) {
  if (d.data != cladeRoot) {
    if (_cladeRoot == null) {
      _cladeRoot = cladeRoot;
    }
    cladeRoot = d.data;
  }
  else {
    cladeRoot = _cladeRoot;
    _cladeRoot = null;
  }
  // var  vertical = d3.select("#svg2");
  // vertical.select("#vertical_rect").remove();
  updateLayout(cladeRoot);
}

svgEvents.setcolour = function(d) {
  var myStorage = window.localStorage;
  var fils = d.descendants();
  var phylumName = "Undefined";
  if (d.data.name !== undefined){
    phylumName = d.data.name;
  }
  fils.forEach(function (d,i ) {
     if (d.data.name !== undefined)
     {
      myStorage.setItem(d.data.name,colourLeaf);
      myStorage.setItem("phylum_"+d.data.name,phylumName);
     }
  })
  updateLayout(cladeRoot);
}
