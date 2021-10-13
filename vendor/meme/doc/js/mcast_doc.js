//
// mcast_doc.js
//

// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_mcast_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'mcast-cluster-score':
      return( `
        The score that the hidden Markov model created by MCAST assigned to the motif cluster.<br>
        This is the sum of the scores of the individual motif matches in the cluster, plus a
        gap penalty, <i>g</i>, multiplied by the total size of the inter-motif gaps
        in the cluster.  Individual motif match scores are log2(<i>P(s)/p</i>), where <i>s</i> is the
        <a href="javascript:help_refine('pop_motif_match_score')">log-odds score</a>
        of the motif match, <i>P(s)</i> is the
        <a href="javascript:help_refine('pop_motif_match_pvalue')"><i>p</i>-value</a>
        of the motif match, and <i>p</i> is the user-specified <i>p</i>-value threshold (default: 0.0005).
        <div class="active" id="pop_motif_match_score_act"></div>
        <div class="active" id="pop_motif_match_pvalue_act"></div>
      `);
    case 'mcast-cluster-p-value':
      return( `
        The <i>p</i>-value of the motif cluster score.<br>
        MCAST estimates <i>p</i>-values by fitting an exponential distribution
        to the observed motif cluster scores.
      `);
    case 'mcast-cluster-E-value':
      return( `
        The <i>E</i>-value of the motif cluster score.<br>
        The <i>E</i>-value is an estimate of the <i>number</i> of false positive matches
        with <i>p</i>-values at least as good as this match\'s.
        MCAST estimates this by multiplying the motif cluster score <i>p</i>-value
        times the (estimated) number of random matches found in the search.
      `);
    case 'mcast-results-tsv':
      return(`
        <p>
          MCAST outputs a tab-separated values (TSV) file ('mcast.tsv') that contains one line for each
          significant match of a cluster of motifs to a sequence region.
          The lines are sorted in order of decreasing statistical significance (increasing <i>p</i>-value).
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.  The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>pattern_name</td> <td>A unique name that MCAST generates for the match. </td> </tr>
          <tr> <td>2</td> <td>sequence_name</td> <td> ` + get_doc_text('shared', 'sequence-id') + ` --OR-- the `
            + get_doc_text('shared', 'sequence-name') + `</td> </tr>
          <tr> <td>3</td> <td>start</td> <td> ` + get_doc_text('shared', 'match-start-seq', 'matched sequence region') + ` --OR-- `
            + get_doc_text('shared', 'match-start-genomic', 'motif occurrence')
            + get_doc_text('shared', 'parse-genomic-coord', 'The latter case occurs when MCAST') + `</td> </tr>
          <tr> <td>4</td> <td>stop</td> <td> ` + get_doc_text('shared', 'match-stop-seq', 'matched sequence region') + ` --OR-- `
            + get_doc_text('shared', 'match-stop-genomic', 'motif occurrence')
            + get_doc_text('shared', 'parse-genomic-coord', 'The latter case occurs when MCAST') + `</td> </tr>
          <tr> <td>5</td> <td>score</td> <td> ` + get_doc_text('mcast', 'mcast-cluster-score') + ` </td> </tr>
          <tr> <td>6</td> <td>p-value</td> <td> ` + get_doc_text('mcast', 'mcast-cluster-p-value') + ` </td> </tr>
          <tr> <td>7</td> <td>E-value</td> <td> ` + get_doc_text('mcast', 'mcast-cluster-E-value') + ` </td> </tr>
          <tr> <td>8</td> <td>q-value</td> <td> ` + get_doc_text('shared', 'bh-q-value', 'MCAST') + `
          <br><br><b>Note:</b> This column is empty if you used the <code>--text</code> switch.</td> </tr>
          <tr> <td>9</td> <td>sequence</td> <td>The region of the sequence matched to the motif cluster.</td> </tr>
        </table>
      `);
    case 'mcast-results-gff3':
      return(`
        <p>
          MCAST outputs a GFF3</a> file ('mcast.gff') that contains one line for each
          significant match to a motif cluster.
        </p>
        <p>
          The GFF3 format is described <a href="http://gmod.org/wiki/GFF3">here</a>.
          The 'score' reported in the MCAST GFF3 output</a> (in column 6) is<br/>
          &nbsp;&nbsp;&nbsp;&nbsp;<code>min(1000, -10*(log10(pvalue)))</code>, <br/>
          and the 'unique identifier' ('ID=' tag in column 9) is the value of the
          &lt;pattern_name&gt; field in the MCAST results TSV format.  Following the
          unique identifier in column 9, the <i>p</i>-value, <i>E</i>-value and q-value
          of the match are given.
        </p>
      `);
    default:
      return("Error--Unrecognized mcast doc_type: " + doc_type);
  }
} // get_mcast_doc_text
