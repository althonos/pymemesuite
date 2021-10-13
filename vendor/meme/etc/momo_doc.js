//
// momo_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_momo_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_momo_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    // MoMo output fields.
    case 'momo-logo':
      return(`
        MoMo creates a sequence logo for each motif it discovers using
        the <a href="` + site_url + `/doc/ceqlogo.html">ceqlogo</a> utility.
      `);
    case 'momo-mod':
      return(`
        The post-translationally modified residue located
        at the center of the motif.
      `);
    case 'momo-motif':
      return(`
        (A string similar to) a regular expression describing the motif.
        Lower case 'x' represents a match to any residue the non-central
        positions.
        The central peptide and its modification weight (if any) is surrounded
        by '_' characters. If option <code>--single-motif-per-mass</code>
        was specified, the central peptide will be represented by an uppercase 'X',
        which also represents a match to any residue.  With the simple algorithm,
        all non-central positions will contain an 'x', regardless of which residues appeared in
        the input sequences.  (See the Regular Expression column of the output
        for further information.) Motifs are reported in the order
        in which they are found by MoMo, and are not sorted in any way.
      `);
    case 'momo-regexp':
      return(`
        A PERL regular expression suitable for searching for the motif in
        protein sequences using
        <a href="` + site_url + `/doc/fasta-grep.html">fasta-grep</a>
        or PERL scripts. The PERL wildcard character '.' matches any
        resudue. Groups of residues in square brackets '[ ]' match
        any of those residues. For the motif-x and MoDL algorithms,
        this is the regular expression for the motif as output by that
        algorithm.  For the simple algorithm,
        this regular expression includes all of the residues
        in the motif occurrences that match the motif.
      `);
    case 'momo-score':
      return(`
        The algorithm-dependent score of the motif.
        For algorithm motif-x, this is the sum over the significant
        position/residue pairs of -log(p<sub>binomial</sub>);
        for algorithm MoDL, this is the increase in description length in
        bits if the motif were to be removed from the final set of motifs.
      `);
    case 'momo-fg':
      return(`
        fg_match is the number of foreground peptides that match the motif.
        <b>Note: </b>For motif-x, the foreground counts are based on
        <b>all</b> the peptides unless the
        <code>--harvard</code> option is specified, in which
        case the counts are based on the <b>remaining</b> peptides.
      `);
    case 'momo-fg-size':
      return(`
        fg_size is the total number of foreground peptides with the given
        central modification.
      `);
    case 'momo-bg':
      return(`
        bg_match is the number of background peptides that match the motif.
        <b>Note: </b>For motif-x, the background counts are based on
        <b>all</b> the peptides unless the
        <code>--harvard</code> option is specified, in which
        case the counts are based on the <b>remaining</b> peptides.
      `);
    case 'momo-bg-size':
      return(`
        bg_size is the total number of background peptides with the given
        central modification.
      `);
    case 'momo-fold':
      return(`
        The fold enrichment of the foreground matches vs. the background
        matches.  This is equal to: <br>
        &nbsp;&nbsp;(fg_match / fg_size) / (bg_match / bg_size).
      `);
    case 'momo-unadjusted-p':
      return(`
        The <i>p</i>-value of the Fisher Exact test on the
        enrichment of the motif in the foreground vs. the background peptides.
        This value does <b>not</b> accurately represent the statistical
        significance of the motif, and should be interpreted
        as a <b>score</b> only.
      `);
    case 'momo-n-tests':
      return(`
        The number of independent tests performed in the search for the motif,
        which is the number of position/residue pairs algorithm motif-x
        tested for statistical significance.  <b>Note:</b> This field only appears for
        motif-x, and only when the <code>--db-background</code> option is <b>not</b> used.
      `);
    case 'momo-adjusted-p':
      return(`
        The <i>p</i>-value of the Fisher Exact test on the
        enrichment of the motif in the foreground vs. the background
        peptides, adjusted for the number of independent tests
        (1 - (1-pvalue)<sup>tests</sup>).  <b>Note 1:</b> This value <b>does</b>
        accurately represent the statistical significance of the motif.
        <b>Note 2:</b> This field only appears for motif-x, and only when the
        <span style="white-space: nowrap;"><code>--db-background</code></span>
        option is <b>not</b> used.
      `);
    case 'momo-occurrences':
      return(`
        The modified peptides matching the motif.
        <b>Note: </b>For motif-x, only the the
        occurrences that are not covered by a previous
        motif are shown, so there may be fewer occurrences
        than indicated by <code>FG</code> unless you
        specified the <code>--harvard</code> option.
      `);
    case 'momo-modl-log':
      return(`
        The log of steps taken by the MoDL algorithm in discovering the
        group of motifs.  The log is only shown on the line for the
        first motif in a group output by MoDL.
        The description length (DL) given the motifs in the group is
        reported for each step, followed by the motifs in the group.
        The step (e.g., group of motifs)
        achieving the lowest description length is preceded by an aserisk ('*').
        <p>
        The steps are followed by three lines stating the
        final step number, the description length of the motifs in
        the final group, and the decrease in description length
        relative to the empty group (Step 0).
        </p>
      `);
    case 'momo-tsv-description':
      return(`
        <p>
          MoMo outputs a tab-separated values (TSV) file ('momo.tsv') that
          contains one line for each motif found and that is suitable
          for parsing by scripts and viewing with Excel.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the TSV file in a comment line
          starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
          Not all fields are present for all algorithms, and the field numbers are
          indicated for each algorithm in the first three columns of the table.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th colspan=3>Field Number</th> <th rowspan=2>Field<br>Name</th> <th rowspan=2>Field<br>Contents</th> </tr>
          <tr> <th>motif-x</th> <th>MoDL</th> <th>simple</th> </tr>
          <tr> <td>1</td> <td>1</td> <td>1</td> <td>mod</td> <td>` + get_doc_text('momo', 'momo-mod') + `</td> </tr>
          <tr> <td>2</td> <td>2</td> <td>2</td> <td>motif</td> <td>` + get_doc_text('momo', 'momo-motif') + `</td> </tr>
          <tr> <td>3</td> <td>3</td> <td>3</td> <td>regexp</td> <td>` + get_doc_text('momo', 'momo-regexp') + `</td> </tr>
          <tr> <td>4</td> <td>4</td> <td> </td> <td>score</td> <td>` + get_doc_text('momo', 'momo-score') + `</td> </tr>
          <tr> <td>5</td> <td>5</td> <td>4</td> <td>fg_match</td> <td>` + get_doc_text('momo', 'momo-fg') + `</td> </tr>
          <tr> <td>6</td> <td>6</td> <td> </td> <td>fg_size</td> <td>` + get_doc_text('momo', 'momo-fg-size') + `</td> </tr>
          <tr> <td>7</td> <td>7</td> <td> </td> <td>bg_match</td> <td>` + get_doc_text('momo', 'momo-bg') + `</td> </tr>
          <tr> <td>8</td> <td>8</td> <td> </td> <td>bg_size</td> <td>` + get_doc_text('momo', 'momo-bg-size') + `</td> </tr>
          <tr> <td>9</td> <td>9</td> <td> </td> <td>fg/bg</td> <td>` + get_doc_text('momo', 'momo-fold') + `</td> </tr>
          <tr> <td>10</td> <td>10</td> <td> </td> <td>unadjusted_p-value<td>` + get_doc_text('momo', 'momo-unadjusted-p') + `</td> </tr>
          <tr> <td>11</td> <td> </td> <td> </td> <td>tests</td> <td>` + get_doc_text('momo', 'momo-n-tests') + `</td> </tr>
          <tr> <td>12</td> <td> </td> <td> </td> <td>adjusted_p-value</td> <td>` + get_doc_text('momo', 'momo-adjusted-p') + `</td> </tr>
        </table>
      `);
    case 'momo-meme-output':
      return(`
        <p>
          MoMo outputs a file ('momo.txt') that contains each discovered motif
          in <a href="` + site_url + `/doc/meme-format.html">MEME motif format</a>
          that is suitable for use with other MEME Suite programs.
          MoMo creates these position-frequency matrix motifs by aligning the foreground
          peptides matching the motif and computing the position-frequency matrix.
          No pseudo-counts are added.
        <p>
        </p>
          For motif-x motifs where the adjusted
          <i>p</i>-value can be accurately calculated (e.g., you did not specify the
          <code>--db-background</code> option), MoMo writes the value of
          <code>tests</code> times the <code>unadjusted_p-value</code>
          in the <i>E</i>-value field. In all other cases, MoMo writes
          '1' in the <i>E</i>-value field of the motif.
        </p>
      `);
    default:
      return("Error--Unrecognized momo doc_type: " + doc_type);
  }
} // get_momo_doc_text
