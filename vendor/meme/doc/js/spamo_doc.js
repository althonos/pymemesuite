//
// spamo_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_spamo_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'spamo-results-tsv':
      return(`
        <p>
          SpaMo outputs a tab-separated values (TSV) file ('spamo.tsv') that contains one line for each motif
          found to be significantly enriched.
          The lines are grouped by secondary motif and sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>prim_db</td> <td> ` + get_doc_text('shared', 'motif-db', 'the primary motif.') + `</td> </tr>
          <tr> <td>2</td> <td>prim_id</td> <td> ` + get_doc_text('shared', 'motif-id', 'primary') + `</td> </tr>
          <tr> <td>3</td> <td>prim_alt</td> <td> ` + get_doc_text('shared', 'motif-alt-id', 'primary') + `</td> </tr>
          <tr> <td>4</td> <td>prim_cons</td> <td> ` + get_doc_text('shared', 'motif-cons', 'primary') + `</td> </tr>
          <tr> <td>5</td> <td>sec_db</td> <td> ` + get_doc_text('shared', 'motif-db', 'the secondary motif.') + `</td> </tr>
          <tr> <td>6</td> <td>sec_id</td> <td> ` + get_doc_text('shared', 'motif-id', 'secondary') + `</td> </tr>
          <tr> <td>7</td> <td>sec_alt</td> <td> ` + get_doc_text('shared', 'motif-alt-id', 'secondary') + `</td> </tr>
          <tr> <td>8</td> <td>sec_cons</td> <td> ` + get_doc_text('shared', 'motif-cons', 'secondary') + `</td> </tr>
          <tr> <td>9</td> <td>trim_left</td> <td>Number of positions trimmed from left of secondary motif.</td> </tr>
          <tr> <td>10</td> <td>trim_right</td> <td>Number of positions trimmed from right of secondary motif.</td> </tr>
          <tr> <th colspan=3>If the next three fields are not blank, the motif is redundant with a more significant ('parent') motif.</th> </tr>
          <tr> <td>11</td> <td>red_db</td> <td> ` + get_doc_text('shared', 'motif-db', 'the parent motif.') + `</td> </tr>
          <tr> <td>12</td> <td>red_id</td> <td> ` + get_doc_text('shared', 'motif-id', 'parent') + `</td> </tr>
          <tr> <td>13</td> <td>red_alt</td> <td> ` + get_doc_text('shared', 'motif-alt-id', 'parent') + `</td> </tr>
          <tr> <td>14</td> <td>E-value</td> <td>The expected number motifs that would have least one spacing as enriched as the best spacing for this secondary.
            The <i>E</i>-value is the best spacing <i>p</i>-value multiplied by the number of motifs in the input database(s).</td> </tr>
          <tr> <td>15</td> <td>gap</td> <td>The distance between the edge of the primary and the (trimmed) secondary motif.</td> </tr>
          <tr> <td>16</td> <td>orient</td> <td>The (combination) of quadrants for which occurrences of this spacing are combined.</td> </tr>
          <tr> <td>17</td> <td>count</td> <td>The number of occurrences of the secondary motif with the given spacing and orientation to the primary motif.</td> </tr>
          <tr> <td>18</td> <td>total</td> <td>The total number of occurrences of the secondary motif within the margins around the best primary motif occurrence.</td> </tr>
          <tr> <td>19</td> <td>adj_p-value</td> <td>The <i>p</i>-value of the gap and orientation, adjusted for nine combinations of quadrants times the number of gaps tested (as controlled by the <code>-range</code> option).</td> </tr>
          <tr> <td>20</td> <td>p-value</td> <td>The <i>p</i>-value of the gap and orientation adjusted only for the number of gaps tested.</td> </tr>
        </table>
        <br>
      `);
    case 'spamo-dumpseqs-tsv':
      return(`
        <p>
          By specifying the options <code>--dumpseqs</code> or <code>--dumpsigs</code>
          you can have SpaMo create TSV (tab-separated values) files
          describing the motif matches SpaMo used to make the histograms in its HTML output.
          The files are named
          '<code>seqs_&lt;primary_motif&gt;_&lt;secondary_db&gt;_&lt;secondary_motif&gt;.txt</code>'.
          The rows in each file are sorted by sequence name.
          The first line in the file contains the (tab-separated) names of the fields.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%">
          <tr><th>field</th><th>name</th><th>contents</th></tr>
          <tr><td>1</td><td>matches</td><td>Trimmed lowercase sequence centered on primary match with matches in uppercase.</td></tr>
          <tr><td>2</td><td>sec_pos</td><td>Position of the secondary match within the whole sequence.</td></tr>
          <tr><td>3</td><td>pri_match</td><td>Sequence fragment that the primary matched.</td></tr>
          <tr><td>4</td><td>pri_strand</td><td>Strand of the primary match (+/-).</td></tr>
          <tr><td>5</td><td>sec_match</td><td>Sequence fragment that the secondary matched.</td></tr>
          <tr><td>6</td><td>sec_strand</td><td>Strand of the secondary match (+/-).</td></tr>
          <tr><td>7</td><td>same_opp</td><td>The primary match on the same (s) or opposite (o) strand as the secondary.</td></tr>
          <tr><td>8</td><td>down_up</td><td>The secondary match is downstream (d) or upstream (u) of the primary.</td></tr>
          <tr><td>9</td><td>gap</td><td>The gap between the primary and secondary matches.</td></tr>
          <tr><td>10</td><td>seq_name</td><td>The name of the sequence.</td></tr>
          <tr><td>11</td><td>adj_p-value</td><td>The <i>p</i>-value of the bin containing the match, adjusted for the number of bins.</td></tr>
          <tr><th colspan="3">If the sequence names are in UCSC Genome Browser position
          format (e.g., "chr5:36715616-36715623"), the following additional fields will be present:</th></tr>
          <tr><td>12</td><td>pri_bed_chr</td><td>Position of primary match in BED coordinates.</td></tr>
          <tr><td>13</td><td>pri_bed_start</td><td>"</td></tr>
          <tr><td>14</td><td>pri_bed_end</td><td>"</td></tr>
          <tr><td>15</td><td>pri_browser</td><td>Position of primary match in UCSC Genome Browser coordinates.</td></tr>
          <tr><td>16</td><td>sec_bed_chr</td><td>Position of secondary match in BED coordinates.</td></tr>
          <tr><td>16</td><td>sec_bed_start</td><td>"</td></tr>
          <tr><td>16</td><td>sec_bed_end</td><td>"</td></tr>
          <tr><td>19</td><td>sec_browser</td><td>Position of secondary match in UCSC Genome Browser coordinates.</td></tr>
        </table>
      `);
    default:
      return("Error--Unrecognized spamo_doc_type: " + doc_type);
  }
} // get_spamo_doc_text
