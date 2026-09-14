// Fonction sauve le svg de l'arbre
// ---------------------------------
function saveSVGTree(){
  $("a#tree").css('color', '#FF0000');
  $("a#selecp").css('color', '#000000');
  var style = "\n";
  var img = new Image();
  // prepend style to svg
  const treeSvg = d3.select('#svg1');
  treeSvg.insert('defs',":first-child");
  d3.select("#svg1 defs")
      .append('style')
      .attr('type','text/css')
      .html(style);
  var as_text = new XMLSerializer().serializeToString(treeSvg.node());
  img.src = 'data:image/svg+xml;base64,'+window.btoa(unescape(encodeURIComponent(as_text)));
  // window.open().document.write('<p>Please copy or save the the image (larger images may not display properly on the page) <img src="' + img.src + '"/>');
  window.open().document.write('<p>Please copy or save the the image (larger images may not display properly on the page)</p> <img src="' + img.src + '"/>');
};
// Fonction sauve le svg de l'alignement
// -------------------------------------
function saveSVGAlignment(){
  $("a#tree").css('color', '#000000');
  $("a#selecp").css('color', '#FF0000');
  // Application de styles pour l'export
  $('.seqblock text').css({
    'font-family':'Courier',
    'font-size': '16px',
  });
  $('.seqblock g text').css({
    'letter-spacing':(alignmentLetterSpacing)+'px',
    'font-size': fitPoliceSize+'px',
  });
  var style = "\n";
  var img = new Image();
  // prepend style to svg
  const alignmentSvg = d3.select('#svg2');
  alignmentSvg.insert('defs',":first-child");
  d3.select("#svg2 defs")
      .append('style')
      .attr('type','text/css')
      .html(style);
  var as_text = new XMLSerializer().serializeToString(alignmentSvg.node());
  img.src = 'data:image/svg+xml;base64,'+window.btoa(unescape(encodeURIComponent(as_text)));
  window.open().document.write('<p>Please copy or save the the image (larger images may not display properly on the page)</p> <img src="' + img.src + '"/>');
  // var as_text = new XMLSerializer().serializeToString(alignmentSvg.node());
  // img.src = 'data:image/svg+xml;base64,'+window.btoa(unescape(encodeURIComponent(as_text)));
  // window.open().document.write('<p>Please copy or save the the image (it may not display if it is too large) <img src="' + img.src + '"/>');
}
