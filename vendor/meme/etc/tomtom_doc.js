//
// tomtom_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_tomtom_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    case 'tomtom-offset':
      return( `
        The offset of the query motif relative to the target motif in the optimal alignment.<br>
        A positive value indicates the query is shifted right.
      `);
    case 'tomtom-p-value':
      return( `
        The probability that a random motif of the same width as the target would have an
        optimal alignment with a match score as good or better than the target's.<br>
        Tomtom estimates the <i>p</i>-value using a null model consisting of sampling
        motif columns from all the columns in the set of target motifs.
      `);
    case 'tomtom-E-value':
      return( `
        The expected number of false positives in the matches up to this point.<br>
        Tomtom estimates the <i>E</i>-value by multiplying the <i>p</i>-value by
        the total number of target motifs in all the target databases.
      `);
    case 'tomtom-overlap':
      return( `
        The number of motif columns that overlap in the optimal alignment.
      `);
    case 'tomtom-orientation':
      return( `
        The orientation of the target motif that gave the optimal alignment. ` + extra + `
      `);
    case 'tomtom-results-tsv':
      return(`
        <p>
          Tomtom outputs a tab-separated values (TSV) file ('tomtom.tsv') that contains one line for each motif
          found to be significantly enriched.
          The lines are grouped by query motif and sorted in order of decreasing statistical significance.
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>Query_ID</td> <td> ` + get_doc_text('shared', 'motif-id', 'query') + `</td> </tr>
          <tr> <td>2</td> <td>Target_ID</td> <td> ` + get_doc_text('shared', 'motif-id', 'target') + `</td> </tr>
          <tr> <td>3</td> <td>Optimal_offset</td> <td> ` + get_doc_text('tomtom', 'tomtom-offset') + `</td> </tr>
          <tr> <td>4</td> <td>p-value</td> <td> ` + get_doc_text('tomtom', 'tomtom-p-value') + `</td> </tr>
          <tr> <td>5</td> <td>E-value</td> <td> ` + get_doc_text('tomtom', 'tomtom-E-value') + `</td> </tr>
          <tr> <td>6</td> <td>q-value</td> <td> ` + get_doc_text('shared', 'bh-q-value', 'Tomtom') + `</td> </tr>
          <tr> <td>7</td> <td>Overlap</td> <td> ` + get_doc_text('tomtom', 'tomtom-overlap') + `</td> </tr>
          <tr> <td>8</td> <td>Query_consensus</td> <td>A consensus sequence
                computed from the letter frequencies in the query motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>9</td> <td>Target_consensus</td> <td>A consensus sequence
                computed from the letter frequencies in the target motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>10</td> <td>Orientation</td> <td> ` + get_doc_text('tomtom', 'tomtom-orientation', "<br>A value of '+' means that the target motif is as it appears in the database. A value of '-' means that the reverse complement of the target motif is shown.") + `</td> </tr>
        </table>
      `);
    default:
      html = "Error--Unrecognized tomtom doc_type: " + doc_type;
  }
} // get_tomtom_doc_text
