//
// centrimo_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_centrimo_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    case 'centrimo-adj-pvalue':
      return(`
        The statistical significance of the enrichment of the motif, adjusted for multiple tests. ` + extra + `
      `);
    case 'centrimo-evalue':
      var evalue_html = `
        at least one region as enriched for best matches to the motif as the reported region
        (or would have optimal average distance to the sequence center as low as observed,
        if you used the <code>--cd</code> option).
      `;
      return(`
        The expected number motifs that would have ` + (extra ? extra : evalue_html) + `
        The <i>E</i>-value is the adjusted <i>p</i>-value multiplied by the number of motifs in the
        input files(s).</td> </tr>
      `);
    case 'centrimo-bin-width':
      return(`
        The width (in sequence positions) of the most enriched region (default),
        <b>or</b> two times the average distance between the center of the best site
        and the sequence center if you used the option <code>--cd</code>.
        A best match to the motif is counted as being in the region if
        the center of the motif falls in the region.
      `);
    case 'centrimo-sites-in-bin':
      return(`
        The number of (positive) sequences whose best match to the motif `
        + (extra ? extra : "falls in the reported region (default) or anywhere in the sequence (if you used the option <code>--cd</code>).") + `
        <br>Note: This number may be less than the number of
        (positive) sequences that have a best match in the region.
        The reason for this is that a sequence may have many matches that score
        equally best.  If <i>n</i> matches have the best score in a sequence, 1/<i>n</i> is added to the
        appropriate bin for each match.</td> </tr>
      `);
    case 'centrimo-mult-tests':
      return(`
        This is the number of multiple tests (<i>n</i>) done for this motif.
        It was used to adjust the <i>p</i>-value of a region for
        multiple tests using the formula:
        <br>&nbsp;&nbsp;
          <i>p</i>' = 1 - (1-<i>p</i>)<sup><i>n</i></sup> where <i>p</i> is the unadjusted <i>p</i>-value.
        <br>
        The number of multiple tests is the number of regions
        considered times the number of score thresholds considered.
        It depends on the motif length, sequence length, and the type of
        optimizations being done (central enrichment, local enrichment, central distance or
        score optimization).</td> </tr>
      `);

    case 'centrimo-results-tsv':
      return(` 
        <p>
          CentriMo outputs a tab-separated values (TSV) file ('centrimo.tsv') that contains one line for each
          region found to be significantly enriched for a motif.
          The lines are sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.
          The names and meanings of each of the fields, which depend on whether or not you provide
          control sequences to CentriMo, are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>db_index</td> <td>The index of the motif file that contains the motif.  Motif
              files are numbered in the order the appeared in the command line.</td> </tr>
          <tr> <td>2</td> <td>motif_id</td> <td> ` + get_doc_text('shared', 'motif-id') + `
              If more than one motif has the same ID, CentriMo uses only the first such motif.
              The name is single-quoted and preceded with '+' or '-' if you scanned separately with
              the reverse complement motif (using the <code>--sep</code> option).</td> </tr>
          <tr> <td>3</td> <td>motif_alt_id</td> <td> ` + get_doc_text('shared', 'motif-alt-id') + `</td> </tr>
          <tr> <td>4</td> <td>consensus</td> <td> ` + get_doc_text('shared', 'motif-cons') + `</td> </tr>
          <tr> <td>5</td> <td>E-value</td> <td> ` + get_doc_text('centrimo', 'centrimo-evalue') + `</td> </tr>
          <tr> <td>6</td> <td>adj_p-value</td> <td> ` + get_doc_text('centrimo', 'centrimo-adj-pvalue') + `
              By default, a <i>p</i>-value is calculated by using the one-tailed binomial
              test on the number of sequences with a match to the motif
              that have their best match in the reported region;
              if you provided control sequences, the <i>p</i>-value of Fisher\'s exact test on the enrichment of
              best matches in the positive sequences relative to the negative sequences is computed instead;
              if you used the <code>--cd</code> option, the <i>p</i>-value is the probability that the average
              distance between the best site and the sequence center would be as low or lower than observed,
              computed using the cumulative Bates distribution, optimized over different score thresholds.
              In all cases, the reported <i>p</i>-value has been adjusted for the number of regions
              and/or score thresholds tested.</td> </tr>
          <tr> <td>7</td> <td>log_adj_p-value</td> <td>Log of adjusted <i>p</i>-value.</td> </tr>
          <tr> <td>8</td> <td>bin_location</td> <td>Location of the center of the most enriched region, or
                0 if you used the <code>--cd</code> option.
          <tr> <td>9</td> <td>bin_width</td> <td> ` + get_doc_text('centrimo', 'centrimo-bin-width') + `</td> </tr>
          <tr> <td>10</td> <td>total_width</td> <td>The maximum number of regions possible for this motif
              <br>&nbsp;&nbsp;
              round(sequence_length - motif_length + 1)/2,<br>
              or the number of places the motif will fit if you used the <code>--cd</code> option.</td> </tr>
          <tr> <td>11</td> <td>sites_in_bin</td> <td> ` + get_doc_text('centrimo', 'centrimo-sites-in-bin') + `</td> </tr>
          <tr> <td>12</td> <td>total_sites</td> <td>The number of sequences containing a match to the motif
              above the score threshold.
          <tr> <td>13</td> <td>p_success</td> <td>The probability of a random match falling into the enriched region:
              <br>&nbsp;&nbsp;
              bin_width / total_width</td> </tr>
          <tr> <td>14</td> <td>p-value</td> <td>The uncorrected <i>p</i>-value before it gets adjusted for the
              number of multiple tests to give the adjusted <i>p</i>-value.</td> </tr>
          <tr> <td>15</td> <td>mult_tests</td> <td> ` + get_doc_text('centrimo', 'centrimo-mult-tests') + `</td> </tr>
          <tr> <th colspan=3>The following additional columns are present when you provide control sequences to CentriMo
          (using the <code>--neg</code> option).</th> </tr>
          <tr> <td>16</td> <td>neg_sites_in_bin</td> <td>The number of negative sequences where the best
              match to the motif falls in the reported region.
              This value is rounded but the underlying value may contain fractional counts.
              Note: This number may be less than the number of negative have a best match in the region.
              The reason for this is that a sequence may have many matches that score equally best.
              If n matches have the best score in a sequence, 1/n is added to the
              appropriate bin for each match.</td> </tr>
          <tr> <td>17</td> <td>neg_sites</td> <td>The number of negative sequences containing a match to the
              motif above the minimum score threshold.
              When score optimization is enabled the score threshold may be raised
              higher than the minimum.</td> </tr>
          <tr> <td>18</td> <td>neg_adj_pvalue</td> <td>The probability that any tested region in the negative
              sequences would be as enriched for best matches to this motif
              according to the Binomial test.</td> </tr>
          <tr> <td>19</td> <td>log_neg_adj_pvalue</td> <td>Log of negative adjusted <i>p</i>-value.</td> </tr>
          <tr> <td>20</td> <td>fisher_adj_pvalue</td> <td>Fisher adjusted <i>p</i>-value before it gets adjusted for the
              number of motifs in the input files(s).</td> </tr>
          <tr> <td>21</td> <td>log_fisher_adj_pvalue</td> <td>Log of Fisher adjusted <i>p</i>-value.</td> </tr>
        </table>
      `);

    case 'centrimo-sites-txt':
      return(`
        <p>
          CentriMo outputs a text file ('site_counts.txt') that contains,
          for each motif, pairs of values (bin_position, site_count),
          or triples of values (bin_position, site_count, neg_site_count) if
          you provided control sequences to CentriMo.
          This data can be used to plot the density of motif best matches (sites)
          along the input sequences.  Fractional counts are possible if multiple (n)
          bins contain the best match for a given sequence, with each bin
          receiving an incremental count of 1/n.
        </p>
        <p>
          The data for each motif begins with a header line with the format:
          <br>&nbsp&nbsp
                DB &lt;db_number&gt; MOTIF &lt;id&gt; &lt;alt&gt;
          </br>
          where &lt;id&gt; and &lt;alt&gt; are as described above.
          The following lines (up to the next header line)
          each contain a single value-pair or value-triple for the motif
          named in the header line.
        </p>
      `);

    default:
      Return("Error--Unrecognized centrimo doc_type: " + doc_type);
  }
} // get_centrimo_doc_text
