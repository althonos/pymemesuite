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
function get_discovery_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};

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
        This graph shows the positional distribution of the best matches to the motif in the ` +
        extra + ` sequences. The vertical line in the center of the graph corresponds to the center
        of the sequences.
      `);
    case 'sequences-tsv':
      return(`
        <p>` + extra + ` outputs a tab-separated values (TSV) file ('sequences.tsv') containing one line for
        each sequence classified as 'positive' by ` + extra + ` for each reported motif.
        The lines are grouped by motif, and groups are separated by a line
        starting with the character "#".
        The first line in the file contains the (tab-separated) names of the fields.
        The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>motif_DB</td> <td>` + get_doc_text('shared', 'motif-db', 'the motif.') + `</td> </tr>
          <tr> <td>2</td> <td>motif_ID</td> <td>` + get_doc_text('shared', 'motif-id') + `</td> </tr>
          <tr> <td>3</td> <td>motif_ALT_ID</td> <td>` + get_doc_text('shared', 'motif-alt-id') + `</td> </tr>
          <tr> <td>4</td> <td>seq_ID</td> <td>The ID of the sequence.</td> </tr>
          <tr> <td>5</td> <td>seq_Score</td> <td>` + get_doc_text('shared', 'motif-match-score', "The seq_Score of the sequence is its maximum match score.") + `</td> </tr>
          <tr> <td>6</td> <td>seq_Class</td> <td>Whether the sequence is a true positive, 'tp', or a false positive, 'fp'.</td> </tr>
          <tr> <td>7</td> <td>is_holdout?</td> <td>Whether the sequence was in the holdout set, '1', or not, '0'.</td> </tr>
        </table>
      `);
    default:
      return("Error--Unrecognized discovery doc_type: " + doc_type);
  }
} // get_discovery_doc_text
