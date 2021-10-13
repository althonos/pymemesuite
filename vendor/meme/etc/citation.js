//
// citation.js
//
function get_citation_text(doc_type, extra) {
  var html;

  switch (doc_type) {
    case 'AMA':
      return(get_citation_text("GOMo", extra));
    case 'AME':
      return(extra + `
        <span class="citation">
          Robert C. McLeay and Timothy L. Bailey,
          &quot;Motif Enrichment Analysis: a unified framework and an evaluation on ChIP data&quot;,
          <i>BMC Bioinformatics</i>, <b>11</b>:165, 2010.
          <a href="http://www.biomedcentral.com/1471-2105/11/165">[full text]</a>
        </span>
      `);
    case 'CentriMo':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey and Philip Machanick,
          &quot;Inferring direct DNA binding from ChIP-seq&quot;,
          <i>Nucleic Acids Research</i>, <b>40</b>:e128, 2012.
          <a href="http://nar.oxfordjournals.org/content/40/17/e128">[Full Text]</a>
        </span>
      `);
    case 'DREME':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey,
          &quot;DREME: Motif discovery in transcription factor ChIP-seq data&quot;,
          <i>Bioinformatics</i>, <b>27</b>(12):1653-1659, 2011.
          <a href="http://bioinformatics.oxfordjournals.org/content/27/12/1653">[full text]</a>
        </span>
      `);
    case 'FIMO':
      return(extra + `
        <span class="citation">
          Charles E. Grant, Timothy L. Bailey and William Stafford Noble,
          &quot;FIMO: Scanning for occurrences of a given motif&quot;,
          <i>Bioinformatics</i> <b>27</b>(7):1017-1018, 2011.
          <a href="http://bioinformatics.oxfordjournals.org/content/27/7/1017">[full text]</a>
        </span>
      `);
    case 'GLAM2':
    case 'GLAM2SCAN':
      return(extra + `
        <span class="citation">
          Martin C. Frith, Neil F. W. Saunders, Bostjan Kobe and Timothy L. Bailey,
          &quot;Discovering sequence motifs with arbitrary insertions and deletions&quot;,
          <i>PLoS Computational Biology</i>, <b>4</b>(5):e1000071, 2008.
          <a href="https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000071">[full text]</a>
        </span>
      `);
    case 'GOMo':
      return(extra + `
        <span class="citation">
          Fabian A. Buske, Mikael Bod&eacute;n, Denis C. Bauer and Timothy L. Bailey,
          &quot;Assigning roles to DNA regulatory motifs using comparative genomics&quot;,
          <i>Bioinformatics</i>, <b>26</b>(7), 860-866, 2010.
          <a href="http://bioinformatics.oxfordjournals.org/cgi/content/full/26/7/860">[full text]</a>
        </span>
      `);
    case 'MAST':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey and Michael Gribskov,
          &quot;Combining evidence using p-values: application to sequence homology searches&quot;,
          <i>Bioinformatics</i>, <b>14</b>(1):48-54, 1998.
          <a href="http://bioinformatics.oxfordjournals.org/content/14/1/48">[full text]</a>
        </span>
      `);
    case 'MCAST':
      return(extra + `
        <span class="citation">
          Timothy Bailey and William Stafford Noble,
          &quot;Searching for statistically significant regulatory modules&quot;,
          <i>Bioinformatics (Proceedings of the European Conference on Computational Biology)</i>,
          <b>19</b>(Suppl. 2):ii16-ii25, 2003.
          <a href="http://bioinformatics.oxfordjournals.org/cgi/content/abstract/19/suppl_2/ii16">[full text]</a>
        </span>
      `);
    case 'Meta-MEME':
      return(extra + `
        <span class="citation">
	  William N. Grundy, Timothy L. Bailey, Charles P. Elkan and Michael E. Baker.
	  &quot;Meta-MEME: Motif-based Hidden Markov Models of Protein Families&quot;
	  <i>Computer Applications in the Biological Sciences (CABIOS)</i>,
	  <b>13</b>(4):397-406, 1997.
	  <a href="http://bioinformatics.oxfordjournals.org/content/13/4/397">[full text]</a>
        </span>
      `);
    case 'MEME':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey and Charles Elkan,
          &quot;Fitting a mixture model by expectation maximization to
          discover motifs in biopolymers&quot;,
          <em>Proceedings of the Second International Conference on Intelligent Systems
          for Molecular Biology</em>, pp. 28-36, AAAI Press, Menlo Park, California, 1994.
          <a href="http://www.aaai.org/Papers/ISMB/1994/ISMB94-004.pdf">[full text]</a>
        </span>
      `);
    case 'MEME-ChIP':
      return(extra + `
        <span class="citation">
          Philip Machanick and Timothy L. Bailey,
          &quot;MEME-ChIP: motif analysis of large DNA datasets&quot;,
          <i>Bioinformatics</i> <b>27</b>(12):1696-1697, 2011.
        <a href="http://bioinformatics.oxfordjournals.org/content/27/12/1696.full">[full text]</a>
        </span>
      `);
    case 'MEME_SUITE':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey, James Johnson, Charles E. Grant, William S. Noble,
          &quot;The MEME Suite&quot;,
          <i>Nucleic Acids Research</i>, <b>43</b>(W1):W39-W49, 2015.
          <a href="https://academic.oup.com/nar/article/43/W1/W39/2467905">[full text]</a>
        </span>
      `);
    case 'MoMo':
      return(extra + `
        <span class="citation">
          Alice Cheng, Charles Grant, Timothy L. Bailey and William Noble,
          &quot;MoMo: Discovery of statistically significant post-translational modification motifs&quot;, 
          <i>Bioinformatics</i>, <b>35</b>(16):2774-2782, 2018.
          <a href="https://doi.org/10.1093/bioinformatics/bty1058">[full text]</a>
        </span>
      `);
    case 'PSPs':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey, Mikael Bod&eacute;n, Tom Whitington and Philip Machanick,
          &quot;The value of position-specific priors in motif discovery using MEME&quot;,
          <i>BMC Bioinformatics</i>, <b>11</b>(1):179, 2010.
          <a href="http://www.biomedcentral.com/1471-2105/11/179">[full text]</a>
        </span>
      `);
    case 'SEA':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey and Charles E. Grant, &quot;SEA: Simple Enrichment Analysis of motifs&quot;,
          <i>BioRxiv</i>, 2021.
        </span>
      `);
    case 'SpaMo':
      return(extra + `
        <span class="citation">
          Tom Whitington, Martin C. Frith, James Johnson and Timothy L. Bailey
          &quot;Inferring transcription factor complexes from ChIP-seq data&quot;,
          <i>Nucleic Acids Res.</i> <b>39</b>(15):e98, 2011.
          <a href="http://nar.oxfordjournals.org/content/39/15/e98">[full text]</a>
        </span>
      `);
    case 'STREME':
      return(extra + `
        <span class="citation">
          Timothy L. Bailey,
          &quot;STREME: accurate and versatile sequence motif discovery&quot;,
          <i>Bioinformatics</i>, Mar. 24, 2021.
          <a href="https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab203/6184861" >[full text]</a>
        </span>
      `);
    case 'Tomtom':
      return(extra + `
        <span class="citation">
          Shobhit Gupta, JA Stamatoyannopolous, Timothy Bailey and William Stafford Noble,
          &quot;Quantifying similarity between motifs&quot;,
          <i>Genome Biology</i>, <b>8</b>(2):R24, 2007.
          <a href="http://genomebiology.com/2007/8/2/R24">[full text]</a>
        </span>
      `);
    case 'T-Gene':
      return(extra + `
        <span class="citation">
          Timothy O'Connor, Charles E. Grant, Mikael Bod&eacute;n, Timothy L. Bailey,
          &quot;T-Gene: Improved target gene prediction&quot;,
            <i>Bioinformatics</i>, <b>36</b>(12):3902-3904, 2020.
          <a href="https://academic.oup.com/bioinformatics/article/36/12/3902/5815978?guestAccessKey=aa625a49-a2aa-4d7a-858e-8bc82867a534">[Full Text]</a>
        </span>
      `);
    case 'XSTREME':
      return(extra + `
        <span class="citation">
          Charles E. Grant and Timothy L. Bailey, &quot;XSTREME: comprehensive motif analysis of biological sequence datasets&quot;,
          <i>BioRxiv</i>, 2021.
        </span>
      `);
    default:
      return("Unknown program: " + doc_type);
  }
} // get_citation_text

//
// Function to replace the innerHTML of element "id" with an HTML paragraph
// containing the text for 'program', which is known to function get_citation_text.
// If "id" is either "citation" or "reference" some extra text is printed.
//
function print_citation(id, program) {
  var extra;
  switch (id) {
    case 'citation':
      extra = "If you use " + program + " in your research, please cite the following paper:<br>";
      break;
    case 'reference':
      extra = "<h5>Reference</h5>";
      break;
    default:
      extra = "";
      break;
  };
  var html = get_citation_text(program, extra);
  document.getElementById(id).insertAdjacentHTML('beforeend', html);
} // print_citation

// 
// Function to convert a citation for a program to a C #define statement.
//
function print_citation_define(lang, pgm) {
  var citation = get_citation_text(pgm, '');
  citation = citation.replace(/<[^>]*>/g, '');
  citation = citation.replace(/\[.*\]/g, '');
  citation = citation.replace(/\n\s*/g, '\\n');
  citation = citation.replace(/&quot;/g, '\\"');
  citation = citation.replace(/&eacute;/g, 'e');
  citation = citation.replace(/^\\n/, '');
  pgm = pgm.replace(/-/, '');
  citation = "If you use this program in your research, please cite:\\n\\n" + citation;
  if (lang == "C") {
    citation = "#define " + pgm + '_CITE "' + citation + '"';
  } else if (lang == "perl") {
    citation = '"' + pgm + '" => "' + citation + '",';
  }
  return(citation);
} // print_citation_define

//
// Main program (for use with nodejs "node" javascript engine)
// to create citation.js.h and citation.pm from citation.js.
// The command line:
//   node citation.js C > citation.js.h
// will output the C #define statements for each of the
// programs listed below, defining macros <program>_CITE.
// The command line:
//   node citation.js perl > citation.js.pm
// will output perl hash <program> => text
//
//
if (typeof process !== 'undefined') {
  var lang = process.argv[2];
  var programs = ['AMA', 'AME', 'CentriMo', 'DREME', 'FIMO', 'GLAM2', 
    'GLAM2SCAN', 'GOMo', 'MAST', 'MCAST', 'Meta-MEME', 'MEME',
    'MEME-ChIP', 'MEME_SUITE', 'MoMo', 'PSPs', 'SEA', 'SpaMo',
    'STREME', 'Tomtom', 'T-Gene', 'XSTREME'];

  if (lang == "C") {
    console.log("// Do not edit this file.  It is created from etc/citation.js.");
    console.log("#ifndef citation_js_h\n#define citation_js_h\n");
    for (var i=0; i<programs.length; i++) {
      console.log(print_citation_define(lang, programs[i]));
    }
    console.log("\n#endif");
  } else if (lang == "perl") {
    console.log("# Do not edit this file.  It is created from etc/citation.js.");
    console.log("package Citation;");
    console.log("sub cite {\n  my ($pgm) = @_;\n  return $citation{$pgm};\n}");
    console.log("%citation = (");
    for (var i=0; i<programs.length; i++) {
      console.log(print_citation_define(lang, programs[i]));
    }
    console.log(");");
  }
}
