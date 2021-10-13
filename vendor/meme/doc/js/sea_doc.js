//
// sea_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_sea_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    // SEA output fields
    case 'sea-rank':
      return(`
	The rank of the significance (<i>E</i>-value) of the motif in the sorted results.
      `);
    case 'sea-tp':
      return(`
        The number of <b>primary</b> sequences matching the motif / the
        number of primary sequences (the percentage of primary sequences
        matching the motif).
        <br>
        A sequence is said to match the motif if some position
        within it has match score greater than or equal to
        the optimal threshold (Score Threshold).
      `);
    case 'sea-tp%':
      return(`
	The percentage of primary sequences matching the motif with 
	scores greater than or equal to the optimal match score threshold.
      `);
    case 'sea-fp':
      return(`
	The number of <b>control</b> sequences matching the motif / the
        number of control sequences (the percentage of control sequences
        matching the motif).
        <br>
        A sequence is said to match the motif if some position
        within it has match score greater than or equal to
        the optimal threshold (Score Threshold).
      `);
    case 'sea-fp%':
      return(`
	The percentage of control sequences matching the motif with
	scores greater than or equal to the optimal match score threshold.
      `);
    case 'sea-enr_ratio':
      return(`
        The relative enrichment ratio of the motif in the primary vs. control sequences, defined as
        <br>&nbsp;&nbsp;&nbsp;&nbsp;Ratio = ((TP+1)/(NPOS+1)) / ((FP+1)/(NNEG+1)), 
        <br>where NPOS is the number of primary sequences in the input, and NNEG is the
        number of control sequences in the input.
      `);
    case 'sea-score_thr':
      return(`
	The match score threshold giving the optimal <i>p</i>-value.
        This is the score threshold used by SEA to determine the values of "TP" and "FP".
        SEA uses the hold-out sequence set to determine the score threshold unless there
	are too few sequences in the input.
      `);
    case 'sea-pvalue':
      return(`
        The optimal enrichment <i>p</i>-value of the motif according to the statistical test.
        Not adjusted for the number of motifs. 
        <br><br>
	If there are too few sequences in the input to create a hold-out sequence set,
	motifs with the same values of TP and FP can have different <i>p</i>-values.
	This is because a variable correction has been applied to adjust for testing 
	multiple score thresholds.
      `);
    case 'sea-log_pvalue':
      return(`
	The logarithm of PVALUE.
      `);
    case 'sea-evalue':
      return(`
        The expected number of random motifs that would be as enriched in the
        primary sequences as this one.  The <i>E</i>-value is the <i>p</i>-value
        multiplied by the number of motifs in the motif file(s).
      `);
    case 'sea-log_evalue':
      return(`
	The logarithm of EVALUE.
      `);
    case 'sea-qvalue':
      return(
        get_doc_text('shared', 'bh-q-value', 'SEA', 'motif') + `
        <br><br><b>Note: </b>Motif q-values smaller than the smallest positive double precision
        number (about 2.2e-308) have been replaced by the motif <i>E</i>-value.
      `);
    case 'sea-log_qvalue':
      return(`
	The logarithm of QVALUE.
      `);
    case 'sea-tsv-description':
      return(`
        <p>
          SEA outputs a tab-separated values (TSV) file ('sea.tsv') that
          contains one line for each motif found to be significantly enriched,
          sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the TSV file in a comment
          line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr><th>field</th><th>name</th><th>contents</th><tr>
          <tr><td>1</td><td>RANK</td><td> ` + get_doc_text('sea', 'sea-rank') + `</td></tr>
          <tr><td>2</td><td>DB</td><td> ` + get_doc_text('shared', 'motif-db', 'the motif.') + `</td></tr>
          <tr><td>3</td><td>ID</td><td> ` + get_doc_text('shared', 'motif-id') + `</td></tr>
          <tr><td>4</td><td>ALT_ID</td><td> ` + get_doc_text('shared', 'motif-alt-id') + `</td></tr>
          <tr><td>5</td><td>CONSENSUS</td><td> ` + get_doc_text('shared', 'motif-cons') + `</td></tr>
          <tr><td>6</td><td>TP</td><td> ` + get_doc_text('sea', 'sea-tp') + `</td><tr>
          <tr><td>7</td><td>TP%</td><td> ` + get_doc_text('sea', 'sea-tp%') + `</td><tr>
          <tr><td>8</td><td>FP</td><td> ` + get_doc_text('sea', 'sea-fp') + `</td><tr>
          <tr><td>9</td><td>FP%</td><td> ` + get_doc_text('sea', 'sea-fp%') + `</td><tr>
          <tr><td>10</td><td>ENR_RATIO</td><td> ` + get_doc_text('sea', 'sea-enr_ratio') + `</td><tr>
          <tr><td>11</td><td>SCORE_THR</td><td> ` + get_doc_text('sea', 'sea-score_thr') + `</td><tr>
          <tr><td>12</td><td>PVALUE</td><td> ` + get_doc_text('sea', 'sea-pvalue') + `</td></tr>
          <tr><td>13</td><td>LOG_PVALUE</td><td> ` + get_doc_text('sea', 'sea-log_pvalue') + `</td></tr>
          <tr><td>16</td><td>EVALUE</td> <td> ` + get_doc_text('sea', 'sea-evalue') + `</td></tr>
          <tr><td>17</td><td>LOG_EVALUE</td><td> ` + get_doc_text('sea', 'sea-log_evalue') + `</td></tr>
          <tr><td>18</td><td>QVALUE</td> <td> ` + get_doc_text('sea', 'sea-qvalue') + `</td></tr>
          <tr><td>18</td><td>LOG_QVALUE</td> <td> ` + get_doc_text('sea', 'sea-log_qvalue', 'SEA') + `</td></tr>
        </table>
      `);
    case 'sea-sequences-tsv':
      return(`
        <p>SEA outputs a tab-separated values (TSV) file ('sequences.tsv') containing one line for
        each sequence with a site whose score passes the motif's score threshold for each enriched motif 
        reported by SEA.
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
          <tr> <td>5</td> <td>seq_Score</td> <td>` + get_doc_text('shared', 'motif-match-score', '', 'The seq_Score of a sequence is its maximum motif match score over all sequence positions.') + `</td> </tr>
          <tr> <td>6</td> <td>seq_Class</td> <td>Whether the sequence is a true positive, 'tp', or a false positive, 'fp'.</td> </tr>
          <tr> <td>7</td> <td>is_holdout?</td> <td>Whether the sequence was in the holdout set, '1', or not, '0'.</td> </tr>
        </table>
      `);
    default:
      return("Error--Unrecognized sea_doc_type: " + doc_type);
  }
} // get_sea_doc_text
