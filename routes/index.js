// 
// Modified by Genna BEN HASSEN, Camille SIHARATH and Grégoire ALIZADEH NOIRET, 16 December 2021
// Modified by Grégoire ALIZADEH NOIRET, 19 July 2022
// Modified by Laurent GUÉGUEN, 6 March 2023
// 

const bodyParser = require('body-parser');
const {exec} = require('child_process');
const express = require('express');
const router = express.Router();
const fs = require('fs');
const multer = require('multer');
const upload = multer({ dest: 'uploads/' });
const url = require('url');
const upload_dir = 'uploads/';
const xml_dir = 'uploads/';

router.use(bodyParser.urlencoded({ extended: true })); // support encoded bodies
router.use(bodyParser.json()); // support json encoded bodies
if (typeof localStorage === "undefined" || localStorage === null) {
  var LocalStorage = require('node-localstorage').LocalStorage;
  localStorage = new LocalStorage('./scratch');
}

// Home
// ----
router.get('/', function(req, res) {
  res.render('index.ejs', {title: 'Branch Site Interface'});
});

// POST upload_files
// -----------
var fileRemoval = false;
router.post("/upload_files", upload.fields([
  {name: 'file_t'}, 
  {name: 'file_a'}, 
  {name: 'file_r'}
]), 
(req, res) => {
  console.log('\n'+'Uploading files');
  var currentTimestampMs = Date.now(); // Date.now() is the current timestamp in milliseconds
  var fname_t = currentTimestampMs+'-t'+req.files.file_t[0].filename.substring(0, 3);
  var fname_a = currentTimestampMs+'-a'+req.files.file_a[0].filename.substring(0, 3);
  var fname_r = currentTimestampMs+'-r'+req.files.file_r[0].filename.substring(0, 3);
  fs.rename(req.files.file_t[0].path, upload_dir+fname_t, function(err) {
    if (err) {
      console.log('ERROR:', err);
    }
  });
  fs.rename(req.files.file_a[0].path, upload_dir+fname_a, function(err) {
    if (err) {
      console.log('ERROR:', err);
    }
  });
  fs.rename(req.files.file_r[0].path, upload_dir+fname_r, function(err) {
    if (err) {
      console.log('ERROR:', err);
    }
  });

  var fname_xml = currentTimestampMs
    + '-'+fname_t.split('-')[1]
    + '-'+fname_a.split('-')[1]
    + '-'+fname_r.split('-')[1];
  var full_path_xml = xml_dir + fname_xml;
  var statcol = req.body.statcol;
  var nostat = req.body.nostat;
  var resultsType = req.body.resultsType;
  var branchSite;
  var logBranchLength = req.body.logBranchLength;
  var isNuc = req.body.isNuc;

  isNuc = (isNuc != undefined ? true : false);
  logBranchLength = (logBranchLength != undefined ? true : false);
  branchSite = (resultsType == 'branchSiteMode' ? true : false);
  
  // Generate XML file from data
  console.log('Generating PhyloXML tree');
  exec('python3 genere_xml.py'
      +' -t '+upload_dir+fname_t
      +' -a '+upload_dir+fname_a
      +' -r '+upload_dir+fname_r
      +' -o '+xml_dir+fname_xml
      +' -e "\''+statcol+'\'"'
      +' -n '+nostat
      +(branchSite?' -b ':''),
  (error, stdout, stderr) => {
    if (error) {
        console.log(`error: ${error.message}`);
//        error.message = error.stack;
      res.render('error.ejs', {message:"Error during phyloxml generation.",error:error});
    } else {
      if (fileRemoval) {
        console.log('Deleting data files');
        exec('rm'
          +' '+upload_dir+fname_t
          +' '+upload_dir+fname_a
          +' '+upload_dir+fname_r,
          (error, stdout, stderr) => {
            if (error) {
              console.log(`error: ${error.message}`);
            }
            if (stderr) {
              console.log(`error: ${stderr}`);
            }
            console.log(`${stdout}`);
          }
        )
      }
    }
    if (stderr) {
      console.log(`error: ${stderr}`);
    }
    console.log(`${stdout}`);
    
    // Read XML tree as JSON and display data
    const fname = full_path_xml;
    fs.readFile(fname, 'utf8' , (err, data) => {
      if (err) {
        res.render('error.ejs', {message:"Erreur de lecture",error:err});
      }
      var xml_digester = require("xml-digester");
      var handler = new xml_digester.OrderedElementsHandler("eventType");
      var options = {
        "handler": [{
          "path": "eventsRec/*",
          "handler": handler
        }]
      };
      var digester = xml_digester.XmlDigester(options);
      digester.digest(data, function(err, results) {
        if (err) {
          console.log(err);
          return;
        }
        var JSONtree = JSON.stringify(results);
        var JSONpattern = JSON.stringify(""); // Sequence to highlight
        console.log('Rendering view');
        res.render('displaytree.ejs',
        {
          arbre: JSONtree,
          pattern: JSONpattern,
          branchSite: branchSite,
          logBranchLength: logBranchLength,
          isNuc: isNuc
        });
        if (fileRemoval) {
          console.log('Deleting XML file');
          exec('rm'
            +' '+xml_dir+fname_xml,
            (error, stdout, stderr) => {
              if (error) {
                console.log(`error: ${error.message}`);
              }
              if (stderr) {
                console.log(`error: ${stderr}`);
              }
              console.log(`${stdout}`);
            }
          )
        }
      });
    });
  });
});

// GET display_example
// -------------------
router.get('/display_example', function(req, res) {
  // Read XML tree as JSON and display data

  const fname = 'examples/basic_example/4_basic_example_branch_site.xml';

  fs.readFile(fname, 'utf8' , (err, data) => {
    if (err) {
      res.render('error.ejs', {message:"Erreur de lecture",error:err});
    }
    var xml_digester = require("xml-digester");
    var handler = new xml_digester.OrderedElementsHandler("eventType");
    var options = {
      "handler": [{
        "path": "eventsRec/*",
        "handler": handler
      }]
    };
    var digester = xml_digester.XmlDigester(options);
    digester.digest(data, function(err, results) {
      if (err) {
        console.log(err);
        return;
      }
      var JSONtree = JSON.stringify(results);
      var JSONpattern = JSON.stringify(""); // Sequence to highlight
      console.log('Rendering view');
      res.render('displaytree.ejs',
      {
        arbre: JSONtree,
        pattern: JSONpattern,
        branchSite: true,
        logBranchLength: true
      });
    });
  });
});

module.exports = router;
