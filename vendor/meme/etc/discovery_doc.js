//
// discovery_doc.js
// Documentation common to motif discovery tools.
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_discovery_doc_text(doc_type, extra, extra2) {
  if (extra == undefined) {extra = ""};
  if (extra2 == undefined) {extra2 = ""};

  switch (doc_type) {
    case 'motif_logo':
      return(`
	The sequence logo of the motif.
	  The rules for construction logos are given in
	  the <i>Description</i> section of the documentation for the MEME Suite utility
          <a href="` + extra + `/doc/ceqlogo.html#description">ceqlogo</a>.
      `);
    case 'motif_rc_logo':
      return(`
	The sequence logo of the reverse complement motif.
	  The rules for construction logos are given in
	  the <i>Description</i> section of the documentation for the MEME Suite utility
          <a href="` + extra + `/doc/ceqlogo.html#description">ceqlogo</a>.
      `);
    case 'more':
      return(`
	Click on the blue symbol below to reveal detailed information about the motif.
      `);
    case 'submit_dl':
      return(`
	Click on the blue symbol below to reveal options allowing you
	  to submit this motif to another MEME Suite motif analysis program, to download this
	  motif in various text formats, or to download a sequence "logo" of
	  this motif PNG or EPS format.</p>
	  <h5>Supported Programs</h5>
	  <dl>
	    <dt>Tomtom</dt>
	    <dd>Tomtom is a tool for searching for similar known motifs.
	      [<a href="` + extra + `/doc/tomtom.html?man_type=web">manual</a>]</dd>
	    <dt>MAST</dt>
	    <dd>MAST is a tool for searching biological sequence databases for
	      sequences that contain one or more of a group of known motifs.
	      [<a href="` + extra + `/doc/mast.html?man_type=web">manual</a>]</dd>
	    <dt>FIMO</dt>
	    <dd>FIMO is a tool for searching biological sequence databases for
	      sequences that contain one or more known motifs.
	      [<a href="` + extra + `/doc/fimo.html?man_type=web">manual</a>]</dd>
	    <dt>GOMo</dt>
	    <dd>GOMo is a tool for identifying possible roles (Gene Ontology
	      terms) for DNA binding motifs.
	      [<a href="` + extra + `/doc/gomo.html?man_type=web">manual</a>]</dd>
	    <dt>SpaMo</dt>
	    <dd>SpaMo is a tool for inferring possible transcription factor
	      complexes by finding motifs with enriched spacings.
	      [<a href="` + extra + `/doc/spamo.html?man_type=web">manual</a>]</dd>
	  </dl>
      `);
    case 'site_distr':
      return(`
        This plot shows the positional distribution of the best match to the motif in the ` +
        extra + ` sequences. 
        Only matches with scores at least the ` +
        extra2 + ` score threshold are considered.
        The plot is smoothed with a triangular function whose width is 5% of the maximum ` +
	extra + ` sequence length.
        The position of the dotted vertical line indicates whether the sequences were
	aligned on their left ends, centers, or right ends, respectively.  
      `);
    case 'site_hist':
      return(`
        This histogram shows the distribution of the <b>number</b> of matches to the motif in the ` +
        extra + ` sequences with at least one match.
        Only matches with scores at least the ` +
        extra2 + ` score threshold are considered.
      `);
    default:
      return("Error--Unrecognized discovery doc_type: " + doc_type);
  }
} // get_discovery_doc_text
