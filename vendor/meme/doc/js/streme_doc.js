//
// streme_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_streme_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'streme-motif-id':
      return(`
	The name of the motif uses the IUPAC codes for nucleotides or proteins.
	Letters representing multiple nucleotides are used in nucleotide motif
	positions where several nucleotides are favored.
	The name of the motif is &lt;index&gt;-&lt;consensus&gt;,
	where &lt;index&gt; is the rank of the motif according to P-value or Score,
	and &lt;consensus&gt; is an approximation of the motif by an IUPAC sequence.
    `);
    case 'streme-motif-alt-id':
      return(`
	The alternate name of the motif is STREME-&lt;index&gt;,
	where &lt;index&gt; is the rank of the motif according to P-value or Score.
    `);
    case 'streme-pvalue':
      return(`
        The <i>p</i>-value of the motif based on applying the appropriate statistical test
        to the test set sequences. It is not adjusted for the number of motifs reported by STREME.
        <p>
        If STREME reports a <b>single motif</b>, then 
        the <i>p</i>-value is an <b>accurate estimate</b> of
        the statistical significance of the motif as long as the length
        distributions of the positive and negative sequences are essentially the same.
	However, if STREME reports more than one motif, the <i>p</i>-value does <b>NOT</b>
	completely account for multiple testing, and you should use the <i>E</i>-value
	for assessing whether a motif is truly statistically significant.
        </p>
    `);
    case 'streme-evalue':
      return(`
        The <i>E</i>-value is an <b>accurate estimate</b> of
        the statistical significance of the motif as long as the length
        distributions of the positive and negative sequences are essentially the same.
        The <i>E</i>-value is the <i>p</i>-value multiplied by the number of motifs reported by STREME.  
        It is an estimate of the number of motifs that would be found with enrichment
	as high as this motif in shuffled versions of your positive sequences.
    `);
    case 'streme-score':
      return(`
        The Score is the unadjusted <i>p</i>-value of the motif based on the appropriate test
        applied to the training set sequences.
        Since the Score is not adjusted for multiple tests, it <b>cannot</b>
        be used to determine the statistical significance of the motif.
    `);
    case 'streme-sequences-tsv':
      return(`
        <p>STREME outputs a tab-separated values (TSV) file ('sequences.tsv') containing one line for
        each sequence with a site whose score passes the motif's match threshold for each motif discovered by STREME.
        The lines are grouped by motif, and groups are separated by a line
        starting with the character "#".
        The first line in the file contains the (tab-separated) names of the fields.
        The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>motif_ID</td> <td>` + get_doc_text('streme', 'streme-motif-id') + `</td> </tr>
          <tr> <td>2</td> <td>motif_ALT_ID</td> <td>` + get_doc_text('streme', 'streme-motif-alt-id') + `</td> </tr>
          <tr> <td rowspan=2>3</td> <td>motif_P-value</td> <td>` + get_doc_text('streme', 'streme-pvalue') + `</td> </tr>
          <tr>                      <td>motif_Score</td> <td>` + get_doc_text('streme', 'streme-score') + `</td> </tr>
          <tr> <td>4</td> <td>seq_ID</td> <td>The ID of the sequence.</td> </tr>
          <tr> <td>5</td> <td>seq_Score</td> <td>` + get_doc_text('shared', 'motif-match-score', "", "The seq_Score of a sequence is its maximum motif match score over all sequence positions.") + `</td> </tr>
          <tr> <td>5</td> <td>seq_Class</td> <td>Whether the sequence is a true positive, 'tp', or a false positive, 'fp'.</td> </tr>
          <tr> <td>6</td> <td>is_holdout?</td> <td>Whether the sequence was in the holdout set, '1', or not, '0'.</td> </tr>
        </table>
      `);
    default:
      return("Error--Unrecognized streme_doc_type: " + doc_type);
  }
} // get_streme_doc_text
