//
// ame_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_ame_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    // AME output fields
    case 'ame-pvalue':
      return(`
        The optimal enrichment <i>p</i>-value of the motif according to the statistical test;
        not adjusted for multiple tests.
      `);
    case 'ame-adj-pvalue':
      return(`
        The optimal enrichment <i>p</i>-value of the motif according to the statistical test,
        adjusted for multiple tests using a Bonferroni correction. ` + extra + `
        If the best <i>p</i>-value is <i>p</i> before adjustment,
        and the number of multiple tests is <i>n</i>, then the adjusted
        <i>p</i>-value is <br>&nbsp;&nbsp;&nbsp;&nbsp; 1 - (1-<i>p</i>)<i><sup>n</sup></i>. `);
    case 'ame-evalue':
      return(`
        The expected number of random motifs that would be as enriched in the
        (primary) sequences as this one.  The <i>E</i>-value is the adjusted <i>p</i>-value
        multiplied by the number of motifs in the motif file(s).
      `);
    case 'ame-tsv-description':
      return(`
        <p>
          AME outputs a tab-separated values (TSV) file ('ame.tsv') that
          contains one line for each motif found to be significantly enriched,
          sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the TSV file in a comment
          line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
          Not all fields are present for all types of enrichment analysis,
          and the field numbers are indicated for each type of analysis
          in the first six columns of the table.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th colspan=6>Field Number</th> <th rowspan=2>Field<br>Name</th> <th rowspan=2>Field<br>Contents</th> </tr>
          <tr> <th>fisher</th> <th>ranksum</th> <th>3dmhg</th> <th>4dmhg</th> <th>pearson</th> <th>spearman</th> </tr>
          <tr> <td colspan=6 class='ctr'>1</td> <td>rank</td> <td>The rank of the significance of the motif in the sorted results.</td> </tr>
          <tr> <td colspan=6 class='ctr'>2</td> <td>motif_DB</td> <td> ` + get_doc_text('shared', 'motif-db', 'the motif.') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>2</td> <td>motif_ID</td> <td> ` + get_doc_text('shared', 'motif-id') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>4</td> <td>motif_alt_ID</td> <td> ` + get_doc_text('shared', 'motif-alt-id') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>5</td> <td>consensus</td> <td> ` + get_doc_text('shared', 'motif-cons') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>6</td> <td>p-value</td> <td> ` + get_ame_doc_text('ame-pvalue') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>7</td> <td>adj_p-value</td> <td> ` + get_ame_doc_text('ame-adj-pvalue') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>8</td> <td>E-value</td> <td> ` + get_ame_doc_text('ame-evalue') + `</td> </tr>
          <tr> <td colspan=6 class='ctr'>9</td> <td>tests</td> <td>The number of tests performed; used in correcting the <i>p</i>-value</td> </tr>
          <tr> <td colspan=8></td> </tr>
          <tr> <td colspan=4 class='ctr'>10</td> <td></td> <td></td> <td>FASTA_max</td> <td>The optimal threshold for <b>labeling</b> sequences as positive;
            sequences with FASTA score less than or equal to the threshold are labeled as positive;
            this field will be contain the number of primary sequences if you provided
            a control file (using option <code>--control</code>); this field will contain the size of the
            partition if you specified one (using option <code>--fix-partition</code>).</td> </tr>
          <tr> <td colspan=4 class='ctr'>11</td> <td></td> <td></td> <td>pos</td> <td>The number of sequences <b>labeled</b> as positive.</td> </tr>
          <tr> <td colspan=4 class='ctr'>12</td> <td></td> <td></td> <td>neg</td> <td>The number of sequences <b>labeled</b> as negative.</td> </tr>
          <tr> <td colspan=8></td> </tr>
          <tr> <td class='ctr'>13</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>PWM_min</td> <td>The optimal threshold on PWM score for <b>classifying</b> sequences as positive;
                sequences with PWM score greater than or equal to the threshold are classified as positive.</td> </tr>
          <tr> <td class='ctr'>14</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>TP</td> <td>The number of true positive sequences: sequences both <b>labeled</b> and <b>classified</b> as positive</td> </tr>
          <tr> <td class='ctr'>15</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>%TP</td> <td>The percentage of true positive sequences: percentage of sequences <b>labeled</b> positive and <b>classified</b> as positive.</td> </tr>
          <tr> <td class='ctr'>16</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>FP</td> <td>The number of false positive sequences: sequences <b>labeled</b> negative but <b>classified</b> as positive.</td> </tr>
          <tr> <td class='ctr'>17</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>%TP</td> <td>The percentage of false positive sequences: sequences <b>labeled</b> negative but <b>classified</b> as positive.</td> </tr>
          <tr> <td colspan=8></td> </tr>
          <tr> <td></td> <td class='ctr'>13</td> <td></td> <td></td> <td></td> <td></td> <td>U</td> <td>The value of the Mann-Whitney <i>U</i> statistic.</td> </tr>
          <tr> <td></td> <td class='ctr'>14</td> <td></td> <td></td> <td></td> <td></td> <td>pleft</td> <td>The left-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>15</td> <td></td> <td></td> <td></td> <td></td> <td>pright</td> <td>The right-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>16</td> <td></td> <td></td> <td></td> <td></td> <td>pboth</td> <td>The two-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>17</td> <td></td> <td></td> <td></td> <td></td> <td>adj_pleft</td> <td>The left-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
          <tr> <td></td> <td class='ctr'>18</td> <td></td> <td></td> <td></td> <td></td> <td>adj_pright</td> <td>The right-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
          <tr> <td></td> <td class='ctr'>19</td> <td></td> <td></td> <td></td> <td></td> <td>adj_both</td> <td>The two-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
          <tr> <td colspan=8></td> </tr>
          <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>10</td> <td></td> <td>Pearson_CC</td> <td>The correlation coefficient of the PWM and FASTA scores of positive sequences.</td> </tr>
          <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>11</td> <td></td> <td>mean_squared_error</td> <td>The mean-squared error of the regression line between PWM and FASTA scores.</td> </tr>
          <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>12</td> <td></td> <td>slope</td> <td>The slope of the regression line.</td> </tr>
          <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>13</td> <td></td> <td>intercept</td> <td>The y-intercept of the regression line.</td> </tr>
          <tr> <td colspan=8></td> </tr>
          <tr> <td></td> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>10</td> <td>Spearman_CC</td> <td>The correlation coefficient of the PWM and FASTA ranks of positive sequences.</td> </tr>
        </table>
          </ul>
        </p>
      `);
    case 'ame-sequences-tsv':
      return(`
        <p>AME outputs a tab-separated values (TSV) file ('sequences.tsv') containing one line for
        each sequence classified as 'positive' by AME for each reported motif.
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
          <tr> <td>5</td> <td><i>label_score</i> (either FASTA_score or PWM_score)</td> <td>The value of the score used to label it as positive.</td> </tr>
          <tr> <td>6</td> <td><i>class_score</i> (either PWM_score or FASTA_score)</td> <td>The value of the score used to classify it as positive.</td> </tr>
          <tr> <td>7</td> <td>class</td> <td>Whether the sequence is a true positive, 'tp', or a false positive, 'fp'.</td> </tr>
        </table>
      `);
    default:
      Return("Error--Unrecognized ame doc_type: " + doc_type);
  }
} // get_ame_doc_text
