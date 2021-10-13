//
// meme-chip_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_meme_chip_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'meme-chip-results-tsv':
      return(`
        <p>
          MEME-ChIP outputs a tab-separated values (TSV) file ('summary.tsv') that
          contains one line for each motif found by MEME-ChIP.
          The lines are sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.
          The names and meanings of the fields in the are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>MOTIF_INDEX</td> <td>The index of the motif in the "Motifs in MEME text format" file ('combined.meme')
                output by MEME-ChIP.</td> </tr>
          <tr> <td>2</td> <td>MOTIF_SOURCE</td> <td> ` + get_doc_text('shared', 'motif-db', 'the motif.', 'of the program that found the <i>de novo</i> motif, or the name of') + `</td> </tr>
          <tr> <td>3</td> <td>MOTIF_ID</td> <td> ` + get_doc_text('shared', 'motif-id', '', 'in the motif discovery program output or ') + `</td> </tr>
          <tr> <td>4</td> <td>ALT_ID</td> <td> ` + get_doc_text('shared', 'motif-alt-id', '', 'in the motif discovery program output or ') + `</td> </tr>
          <tr> <td>5</td> <td>CONSENSUS</td> <td>The ID of the <i>de novo</i> motif, or a consensus sequence
                computed from the letter frequencies in the known motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>6</td> <td>WIDTH</td> <td>The width of the motif.</td> </tr>
          <tr> <td>7</td> <td>SITES</td> <td>The number of sites reported by the <i>de novo</i> program, or the number
                of "Total Matches" reported by CentriMo.</td> </tr>
          <tr> <td>8</td> <td>E-VALUE</td> <td>The statistical significance of the motif.</td> </tr>
          <tr> <td>9</td> <td>E-VALUE_SOURCE</td> <td>The program that reported the <i>E</i>-value.</td> </tr>
          <tr> <td>10</td> <td>MOST_SIMILAR_MOTIF</td> <td>The known motif most similar to this motif according to Tomtom.</td> </tr>
          <tr> <td>11</td> <td>URL</td> <td>A link to a description of the most similar motif, or to the known motif.</td> </tr>
        </table>
      `);
    case 'meme-chip-combined-motifs':
      return(`
        <p>
          MEME-ChIP outputs a text file ('combined.meme') containing all the significant motifs found by MEME-ChIP.
          The motifs are in <a href="` + site_url + `/doc/meme-format.html">Minimal MEME Motif format</a>,
          and their IDs correspond to the motif indices given in the "Summary in TSV Format" file ('summary.tsv').
        </p>
        </p>
          <b>Note:</b> The 'nsites=' and 'E=' fields in the motif headers are only
          relevant for the MEME and DREME motifs.  For known motifs, those values do
          not refer to the number of sites in the input sequences.
        </p>
      `);
    default:
      return("Error--Unrecognized meme_chip doc_type: " + doc_type);
  }
} // get_meme-chip_doc_text
