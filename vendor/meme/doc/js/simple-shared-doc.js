//
// simple-shared-doc.js
//

//
// Function to redirect to appropriate doc file.
//
function get_doc_text(pgm, doc_type, extra, extra2) {
  switch (pgm) {
    case 'shared':
      return(get_shared_doc_text(doc_type, extra, extra2));
    case 'ame':
      return(get_ame_doc_text(doc_type, extra, extra2));
    case 'centrimo':
      return(get_centrimo_doc_text(doc_type, extra, extra2));
    case 'discovery':
      return(get_discovery_doc_text(doc_type, extra, extra2));
    case 'fimo':
      return(get_fimo_doc_text(doc_type, extra, extra2));
    case 'gomo':
      return(get_gomo_doc_text(doc_type, extra, extra2));
    case 'mcast':
      return(get_mcast_doc_text(doc_type, extra, extra2));
    case 'meme-chip':
      return(get_meme_chip_doc_text(doc_type, extra, extra2));
    case 'momo':
      return(get_momo_doc_text(doc_type, extra, extra2));
    case 'sea':
      return(get_sea_doc_text(doc_type, extra, extra2));
    case 'spamo':
      return(get_spamo_doc_text(doc_type, extra, extra2));
    case 'streme':
      return(get_streme_doc_text(doc_type, extra, extra2));
    case 'tgene':
      return(get_tgene_doc_text(doc_type, extra, extra2));
    case 'tomtom':
      return(get_tomtom_doc_text(doc_type, extra, extra2));
    case 'xstreme':
      return(get_xstreme_doc_text(doc_type, extra, extra2));
    default:
      return("<b>Unknown program type: <font color=red>" + pgm + "</font></b>");
  }
} // get_doc_text

//
// Function to replace the innerHTML of element "id" with the HTML indicated by "doc_type".
// Easier to read and update than the more flexible approach in shared-doc.js. 
//
function print_doc(id, pgm, doc_type, extra) {
  document.getElementById(id).insertAdjacentHTML('beforeend', get_doc_text(pgm, doc_type, extra));
} // print_doc

//
// Function to replace the innerHTML of element "id" with an HTML paragraph
// containing the text for 'pgm' and 'doc_type'.
// This function can be used in help pop-ups.
//
function print_doc_para(id, pgm, doc_type, extra, extra2) {
  html = "<p>" + get_doc_text(pgm, doc_type, extra, extra2) + "</p>"; 
  document.getElementById(id).insertAdjacentHTML('beforeend', html);
} // print_doc_para

//
// Function to return the Shared HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_shared_doc_text(doc_type, extra, extra2) {
  if (extra == undefined) {extra = ""};
  if (extra2 == undefined) {extra2 = ""};
  switch (doc_type) {
    case 'motif-db':
      return(`
	The name of ` + extra2 + ` a file of motifs ("motif database file") that contains ` + extra + `
      `);
    case 'motif-id':
      return(`
	The name of the ` + extra + ` motif, which is unique ` + extra2 + ` in the motif database file.
      `);
    case 'motif-alt-id':
      return(`
	An alternate name for the ` + extra + ` motif that may be provided ` + extra2 + ` in the motif database file.
      `);
    case 'motif-width':
      return(`
	The width of the motif. No gaps are allowed in motifs supplied to ` + extra + `
        as it only works for motifs of a fixed width.
      `);
    case 'motif-cons':
      return(`
	A consensus sequence computed from the ` + extra + ` motif (as described <a href="#consensus_doc">below</a>).
      `);
    case 'motif-match-score':
     return(`
	` + extra2 + ` The motif match score of a position in a sequence is
	computed by summing the appropriate entry from each column of the
	position-dependent scoring matrix that represents the motif. ` + extra + `
     `);
    case 'motif-match-p-value':
      return(`
	The <i>p</i>-value of a motif match is the probability of a single random
	subsequence of the length of the motif <a href="javascript:help_refine('pop_motif_match_score')">scoring</a>
	at least as well as the observed match.
      `);
    case 'bh-q-value':
      if (extra2 == "") extra2 = "match";
      return(`
	The q-value is the minimum False Discovery Rate (FDR) required to consider this
        ` + extra2 + ` significant.</br>` +
        get_shared_doc_text('bh-q-value-method', extra, extra2) + `
      `);
    case 'bh-q-value-method':
      return(`
        <br>` + extra + ` estimates q-values from all the ` + extra2 + ` <i>p</i>-values 
	using the method proposed by Benjamini & Hochberg (<i>Journal of the Royal Statistical Society B</i>, 57:289-300, 1995).
	See also Storey JD, Tibshirani R. Statistical significance for
	genome-wide studies, <i>Proc. Natl. Acad. Sci. USA</i> (2003) <b>100</b>:9440&ndash;9445.
      `);
    case 'sdb-name':
      return(`
	The name of the (FASTA) sequence database file.
      `);
    case 'sdb-psp':
      return(`
	The name of the position specific priors (PSP) file.
      `);
    case 'sdb-dist':
      return(`
	The name of the binned distribution of priors file.
      `);
    case 'sdb-count':
      return(`
	The number of sequences in the database.
      `);
    case 'sdb-letters':
      return(`
	The number of letters in the sequence database.
      `);
    case 'lastmod':
      return(`
	The date of the last modification to the ` + extra + ` database.
      `);
    case 'sequence-id':
      return(`
        The identifier of the sequence (from the FASTA sequence header line)` + extra + `
      `);
    case 'sequence-name':
      return(`
	` + extra + `name of the sequence extracted from the sequence identifier (in the FASTA sequence header line).<br>
	When you use the <code>--parse-genomic--coord</code> option, the sequence name ends at the
	first colon ':' (if any) present in the sequence\'s FASTA identifier.  Typically this is the
	chromosome or contig name.  With the <code>--parse-genomic--coord</code> option,
	the start and stop positions are in 0-based coordinates relative to the sequence start given 
	in the FASTA sequence identifier (just after the sequence name).</td> </tr>
      `);
    case 'sequence-desc':
      return(`
        The description appearing after the identifier of the sequence in the FASTA header line.
      `);
    case 'sequence-name':
    case 'alph-name':
      return(`
	The name of the alphabet symbol.
      `);
    case 'alph-bg':
      return(`
	The frequency of the alphabet symbol as defined by the background model.
      `);
    case 'match-start-seq':
      return(`
	The start position of the ` + extra + `; 1-based sequence coordinates.
      `);
    case 'match-stop-seq':
      return(`
	The end position of the ` + extra + `; 1-based sequence coordinates.
      `);
    case 'match-start-genomic':
      return(`
	The start position of the ` + extra + `; genomic coordinates.
      `);
    case 'match-stop-genomic':
      return(`
	The end position of the ` + extra + `; genomic coordinates.
      `);
    case 'parse-genomic-coord':
      return(`
	` + extra + ` was run with the <code>--parse-genomic-coord</code> option
	and has split the sequence identifier into sequence name, sequence start and sequence end 
	in genomic coordinates.
      `);
    case 'motif-consensus':
      return(`
        <p id="consensus_doc">
           A <b>consensus sequence</b> is constructed from each column in a
           motif's frequency matrix using the <b>"50% rule"</b>
           as follows:
        </p>
        <ol>
          <li>The letter frequencies in the column are sorted in decreasing order.</li>
          <li>Letters with frequency less 50% of the maximum are discarded.</li>
          <li>The letter used in this position in the consensus sequence is determined
          by the first rule below that applies:</li>
          <ul>
            <li>If there is only one letter left, or if the remaining letters exactly match
            an ambiguous symbol in the alphabet, the <b>letter</b> or <b>ambiguous symbol</b>,
            respectively, is used.</li>
            <li>Otherwise, if the remaining set contains at least 50% of the core
            symbols in the alphabet, the alphabet's <b>wildcard</b>
            (e.g., "N" for DNA or RNA, and "X" for protein) is used.</li>
            <li>Otherwise, the letter with the <b>maximum frequency</b> is used.</li>
          </ul>
        </ol>
      `);
    default:
      return("Error--Unrecognized shared doc_type: " + doc_type);
  }
} // get_shared_doc_text
