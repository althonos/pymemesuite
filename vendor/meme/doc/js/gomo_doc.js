//
// gomo_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_gomo_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    case 'gomo-go-term':
      return(`
        The Gene Ontology Consortium term for a specific role or locality.
        Used for annotating genes with their functions.
      `);
    case 'gomo-score':
      return( `
        A score generated as the <a href="https://en.wikipedia.org/wiki/Geometric_mean">
        geometric mean</a> of <a href="https://en.wikipedia.org/wiki/Mann-Whitney_U_test">rank-sum test(s)</a>
        for the particular Gene Ontology term. The two groups compared by the rank-sum test are scores of genes annotated
        with the GO term and scores of genes not annotated with the GO term.</td> </tr>
      `);
    case 'gomo-p-value':
      return(`
        An empirically generated <i>p</i>-value for the enrichment of the GO term.<br>
        The null hypothesis is that by shuffling the labels on gene scores,
        any possible association between the set of genes that a GO term annotates is destroyed.
        A large number of scores are generated using the null hypothesis and the number of null
        hypothesis scores that are better than each of the real scores is summed and then divided
        by the total null hypothesis scores generated to calculate a <i>p</i>-value.</td> </tr>
      `);
    case 'gomo-results-tsv':
      return(`
        <p>
          GOMo outputs a tab-separated values (TSV) file ('gomo.tsv') that contains one line for each
          motif-GO-term pair found to be significantly enriched.
          The lines are grouped by motif and sorted in order of decreasing statistical significance.
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>Motif_Identifier</td> <td> ` + get_doc_text('shared', 'motif-id') + ` </td> </tr>
          <tr> <td>2</td> <td>GO_Term_Identifier</td> <td> ` + get_gomo_doc_text('gomo-go-term') + ` </td> </tr>
          <tr> <td>3</td> <td>GOMo_Score</td> <td> ` + get_gomo_doc_text('gomo-score') + ` </td> </tr>
          <tr> <td>4</td> <td>p-value</td> <td> ` + get_gomo_doc_text('gomo-p-value') + ` </td> </tr>
          <tr> <td>5</td> <td>q-value</td> <td> ` + get_doc_text('shared', 'bh-q-value', 'GOMo') + ` </td> </tr>
        </table>
      `);
    case 'gomo-results-xml':
      return(`
        <p>
        GOMo outputs an XML file ('gomo.xml') with the following format.
        </p>
        <table class="bordertable" border="1">
          <tr>
            <th>Tag</th><th>Child of</th><th>Description</th>
          </tr>
          <tr>
            <td >&lt;gomo&gt;</td><td >Nothing</td>
            <td>
              Information about this run of GOMo.
              <ul>
                <li>version - The version of GOMo that generated the XML file.</li>
                <li>release - The release date of the version that generated the XML.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td >&lt;program&gt;</td>
            <td >&lt;gomo&gt;</td>
            <td>
              Information about the state of the program when it ran.<br />
              <ul>
                <li>name - name of the program.</li>
                <li>cmd - the command line passed to the program.</li>
                <li>gene_url - the url used to lookup further information on the gene ids.
                The url has ampersands (&amp;) converted into <b>&amp;amp;</b> and the place where
                  the gene ID should be replaced by <b>!!GENEID!!</b> .</li>
                <li>outdir - the output directory that the program wrote to.</li>
                <li>clobber - true if GOMo was allowed to overwrite the output directory.</li>
                <li>text_only - true if GOMo wrote to stdout, in which case this file would
                  not exist so it must be false.</li>
                <li>use_e_values - true if GOMo used <i>E</i>-values (converted from <i>p</i>-values) as
                  input scores, false if GOMo used gene scores.</li>
                <li>score_e_thresh - if GOMo used <i>E</i>-values then this is the threshold that
                  GOMo assumed the worst <i>E</i>-value (<i>p</i>-value = 1.0) for the gene to smooth out noise.</li>
                <li>min_gene_count - the minimum number of genes that a GO term was annotated
                  with before GOMo would calculate a score for it.</li>
                <li>motifs - if present then a space delimited list of the motifs that GOMo
                  calculated a score for, otherwise GOMo scored all motifs.</li>
                <li>shuffle_scores - the number of times GOMo generated a shuffled mapping of
                  gene id to gene id to be used to generate scores from the null model.</li>
                <li>q_threshold - GOMo filtered the results to only show those with a better
                (smaller) q-value.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td>&lt;gomapfile&gt;</td>
            <td>&lt;program&gt;</td>
            <td>
              Information about the GO mapping file.<br />
              <ul>
                <li>path - the path to the mapping file.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td>&lt;seqscorefile&gt;</td>
            <td>&lt;program&gt;</td>
            <td>
              Information about the sequence scoring file.<br />
              <ul>
                <li>path - the path to the sequence scoring file.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td>&lt;motif&gt;</td>
            <td>&lt;gomo&gt;</td>
            <td>
              Information about the motif.<br />
              <ul>
                <li>id - the motif identifier.</li>
                <li>genecount - the number of scored sequences that were used to compute the result.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td>&lt;goterm&gt;</td>
            <td>&lt;motif&gt;</td>
            <td>
              Information about the GO term.<br />
              <ul>
                <li>id - the GO identifier.</li>
                <li>score - the geometric mean across all species of the rank-sum test <i>p</i>-value.</li>
                <li>pvalue - the empirically calculated <i>p</i>-value.</li>
                <li>qvalue - the empirically calculated q-value.</li>
                <li>annotated - the number of genes annotated with the go term.</li>
                <li>group - the subgroup that the term belongs to. For the Gene Ontology
                    b = biological process, c = cellular component and m = molecular function.</li>
                <li>nabove - the number of more general terms that link to this one.</li>
                <li>nbelow - the number of more specific terms that link from this one.</li>
                <li>implied - is the go term implied by other significant go terms?
                  Allows values 'y', 'n' or 'u' (default) for yes, no or unknown.</li>
                <li>description - the GO term description.</li>
              </ul>
            </td>
          </tr>
          <tr>
            <td>&lt;gene&gt;</td>
            <td>&lt;goterm&gt;</td>
            <td>
              Information about the GO term's annotated genes for the primary species.<br />
              <ul>
                <li>id - the gene identifier.</li>
                <li>rank - the rank of the scored gene.</li>
              </ul>
            </td>
          </tr>
        </table>
      `);
    default:
      return("Error--Unrecognized doc_type: " + doc_type);
  }
} // get_gomo_doc_text
