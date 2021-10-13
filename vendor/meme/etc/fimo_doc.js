//
// fimo_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_fimo_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'fimo-results-tsv':
      return(`
        <p>
          FIMO outputs a tab-separated values (TSV) file ('fimo.tsv') that contains one line for each
          significant match to a motif.
          The lines are sorted in order of decreasing statistical significance (increasing <i>p</i>-value).
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.  The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>motif_id</td> <td> ` + get_doc_text('shared', 'motif-id') + `</td> </tr>
          <tr> <td>2</td> <td>motif_alt_id</td> <td> ` + get_doc_text('shared', 'motif-alt-id') + `</td> </tr>
          <tr> <td>3</td> <td>sequence_name</td> <td> ` + get_doc_text('shared', 'sequence-id') + ` --OR-- the `
            + get_doc_text('shared', 'sequence-name') + `</td> </tr>
          <tr> <td>4</td> <td>start</td> <td> ` + get_doc_text('shared', 'match-start-seq', 'motif occurrence') + ` --OR-- `
            + get_doc_text('shared', 'match-start-genomic', 'motif occurrence')
            + get_doc_text('shared', 'parse-genomic-coord', 'The latter case occurs when FIMO') + `</td> </tr>
          <tr> <td>5</td> <td>stop</td> <td> ` + get_doc_text('shared', 'match-stop-seq', 'motif occurrence') + ` --OR-- `
            + get_doc_text('shared', 'match-stop-genomic', 'motif occurrence')
            + get_doc_text('shared', 'parse-genomic-coord', 'The latter case occurs when FIMO') + `</td> </tr>
          <tr> <td>6</td> <td>strand</td> <td>The strand '<code>+</code>' indicates the motif matched the forward
            strand, '<code>-</code>' the reverse strand, and '<code>.</code>'
            indicates strand is not applicable (as for amino acid sequences).</td> </tr>
          <tr> <td>7</td> <td>score</td> <td>The score for the motif occurrence. `
            + get_doc_text('shared', 'motif-match-score') + `</td> </tr>
          <tr> <td>8</td> <td>p-value</td> <td> The <i>p</i>-value of the motif occurrence. `
            + get_doc_text('shared', 'motif-match-p-value') + `
            <div class="active" id="pop_motif_match_score_act"></div></td> </tr>
          <tr> <td>9</td> <td>q-value</td> <td>The q-value of the motif occurrence. `
            + get_doc_text('shared', 'bh-q-value', 'FIMO') + ` <b>Note:</b> This column is empty if
           you used the <code>--text</code> switch.</td> </tr>
          <tr> <td>10</td> <td>sequence</td> <td>The region of the sequence matched to the motif.</td> </tr>
        </table>
      `);
    case 'fimo-results-gff3':
      return(`
        <p>
          FIMO outputs a GFF3</a> file ('fimo.gff') that contains one line for each
          significant match to a motif.
        </p>
        <p>
          The GFF3 format is described <a href="http://gmod.org/wiki/GFF3">here</a>.
          The 'score' reported in the FIMO GFF3 output</a> (in column 6) is<br/>
          &nbsp;&nbsp;&nbsp;&nbsp;<code>min(1000, -10*(log10(pvalue)))</code>, <br/>
          and the 'display name' ('Name=' tag in column 9) is composed of the contents of three
          fields as follows <br/>
          &nbsp;&nbsp;&nbsp;&nbsp;<code>&lt;motif_id&gt;_&lt;sequence_name&gt;&lt;strand&gt;</code>.
        </p>
      `);
    default:
      return("Error--Unrecognized fimo doc_type: " + doc_type);
  }
} // get_fimo_doc_text
