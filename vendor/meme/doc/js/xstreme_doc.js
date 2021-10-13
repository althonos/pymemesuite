//
// xstreme_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_xstreme_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};

  switch (doc_type) {
    case 'xstreme-seed':
      return(`
        This column is 1 if this motif is the seed motif of its cluster; otherwise the column is 0.
        The seed motif in a cluster is either the most enriched <b>discovered</b> motif,
	or, if the cluster contains no discovered motifs, the most enriched <b>known</b> motif
	in the cluster. Enrichment is determined by the SEA motif enrichment analysis program.
      `);
    case 'xstreme-cluster':
      return(`
	This column contains the cluster number of the cluster of similar motifs that
	the motif belongs to.
	Clusters are numbered in order of the enrichment of their seed motif.
      `);
    case 'xstreme-width':
      return(`
	The width of the motif.
      `);
    case 'xstreme-sites':
      return(`
        The number of <b>primary</b> sequences matching the motif
	according to the motif discovery program (discovered motifs), or according to 
	the SEA enrichment analysis program (database motifs).
      `);
    case 'xstreme-rank':
      return(`
	The rank of the enrichment of the motif in the primary sequences
	relative to the control sequences according to the SEA motif 
	enrichment analysis program.  The rank is relative to all of the
        <i>ab initio</i> motifs discovered by STREME and MEME, as well as
        all the known motifs you provided to XSTREME. The rank is
	based on the SEA enrichment <b><i>p</i>-values</b> of the motifs.  
        Motifs with tied SEA enrichment <i>p</i>-values are given the same rank.
        Click on the rank to see the motif in the SEA output.  (There will
	be no link to the SEA output if the SEA <b><i>E</i>-value</b> of the motif
        is larger than the <i>E</i>-value threshold.) <b>Caution:</b>The
	SEA statistical significance statistics (<i>p</i>-, <i>E</i>- and q-values)
        are <b>NOT</b> accurate for <b>discovered</b> motifs because they are 
	based on the same sequences as input to STREME and MEME.
      `);
    case 'xstreme-sea_pvalue':
      return(`
        This is the <i>p</i>-value of the motif according
        to the SEA motif enrichment analysis program.
        The RANK of the motif is based on this value.
	This is <b>NOT</b> an accurate measure of the statistical significance
	of the motif because it is not adjusted for multiple tests.
        SEA estimates the <i>p</i>-value of motifs using
	either the Fisher exact test or the Binomial test, as is described
        <a href="` + site_url + `/doc/sea.html#significance">here</a>.
      `);
    case 'xstreme-significance':
      return(`
        The <i>E</i>-value of the motif according
        to the motif discovery program (<i>ab initio</i> motifs),
        or according to the SEA motif enrichment analysis program
        (known motifs).  It is an <b>accurate</b> measure of the statistical 
        significance of the motif <b>unless</b> the ` + extra + `.
	The different programs estimate <i>E</i>-values as follows:
        <ul>
	  <li>SEA estimates the <i>E</i>-value of the motif using
	  either the Fisher exact test or the Binomial test, as is described
	  <a href="` + site_url + `/doc/sea.html#significance">here</a>.
	  </li>
          <li>STREME also estimates the <i>E</i>-value of the motif using
	  either the Fisher exact test or the Binomial test, as is described
	  <a href="` + site_url + `/doc/streme.html#significance">here</a>.
	  </li>
          <li>MEME estimates the <i>E</i>-value of the motif
	  based on its log likelihood ratio,
	  as is described <a href="` + site_url + `/doc/meme.html#objfun">here</a>.
	  <b>Note: </b>MEME <i>E</i>-values are based on a different null model
	  than STREME and SEA <i>E</i>-values and are therefore not directly comparable to them.
          </li>
        </ul>
      `);
    case 'xstreme-sig_acc':
      return(`
	The E-value is an accurate measure of motif significance if this is "1";
        otherwise the EVALUE is not an accurate measure of motif significance.
      `);
    case 'xstreme-sim_motif_source':
      return(`
	Motifs <b>discovered</b> by a motif discovery program (STREME or MEME) are compared
	with known motifs in the motif database(s) you specified. This column and the next two
	columns identify the <b>most similar</b> known motif to the discovered motif. 
	This column gives the similar motif's database name, and the next two columns
	give its ID and its URL, respectively.  
        Only known motifs with Tomtom similarity <i>E</i>-values of less than 1.0 to the 
	discovered motif will be shown here. 
        <br><br>For <b>known</b> motifs, 
	this column and the next two give the database name, ID and URL of the known motif.
      `);
    case 'xstreme-sim_motif':
      return(`
	For <b>discovered</b> motifs this column gives the ID (with any alternate ID in parentheses) of the 
        <b>most similar</b> known motif in the motif database(s) that you provided to XSTREME.  
        For <b>known</b> motifs it gives the known motif's ID.
      `);
    case 'xstreme-motif_url':
      return(`
	For <b>discovered</b> motifs this columns gives the URL of the
        <b>most similar</b> known motif in the motif database(s) that you provided to XSTREME.
        For <b>known</b> motifs it gives the known motif's URL in its online motif database.
      `);
    case 'xstreme-source':
      return(`
	<p>For <b>discovered</b> motifs, this is the name of the motif linked
	  to it in the output of the motif discovery program, followed by the name of that program.
	</p>
	<p>
	  For <b>known</b> motifs, this is the identifier of the motif linked to its entry
	  at its online motif database website, followed by the name of the motif database.
	</p>
      `);
    case 'xstreme-sim_motif_url':
      return(`
      `);
    case 'xstreme-sim_known_motifs':
      return(`
	Motifs <b>discovered</b> by a motif discovery program (STREME or MEME) are compared
	with known motifs in the motif database(s) you specified. This column
	lists the (up to) three most similar motifs. Only known motifs with
	Tomtom similarity <i>E</i>-values of less than 1.0 to the discovered motif will
	be shown here. Clicking any of these links will show the Tomtom alignment
	of the discovered motif with the known motif.</p> 
      `);
    case 'xstreme-show_clustered':
      return(`
	<p>Clicking here will show you all the motifs found by motif discovery or
	motif enrichment analysis that are significantly similar to the reported
	motif.</p>
	<p>The additional motifs are shown aligned with the reported motif,
	sorted in order of significance according to the SEA motif enrichment
	analysis program.</p>
	<p>To cluster the motifs XSTREME does the following:</p>
	<ol>
	  <li>Start with no groups and all significant reported motifs.</li>
	  <li>Run Tomtom with all significant reported motifs to determine
	  pairwise similarity.</li>
	  <li>Group Highly Similar Motifs---
	  <div style="padding-left: 10px;">
	    While ungrouped motifs:
	    <div style="padding-left: 10px;">
	      Select most significant ungrouped motif.
	      <div style="border: 1px solid gray; color: gray">
		This is called the 'seed' motif for the group and we will call the
		<i>E</i>-value of its seed motif the group's "significance".
	      </div>
	      Form a new group from the seed motif and all other motifs that
	      are not yet in a group and who are strongly similar to the seed
	      motif (default: Tomtom <i>E</i>-value &le; 0.05).
	    </div>
	  </div>
	  </li>
	  <li>Merge Groups---
	  <div style="padding-left: 10px">
	    For each group (most significant to least significant), merge it with
	    any less significant group if all of its motifs are weakly similar to
	    the first group's seed motif (default: Tomtom <i>E</i>-value &le; 0.1).
	  </div>
	  </li>
	</ol>
      `);
    case 'xstreme-fimo_sites':
      return(`
      <p>The motif scanning program FIMO is used to predict motif sites in your
        primary sequences.  Predictions are only made for the 'seed' motif
	in each cluster, and only if it is a <b>discovered</b> motif. `
        + get_doc_text('xstreme', 'xstreme-seed') + 
        `The predicted sites are output in <a href=http://gmod.org/wiki/GFF3>GFF3</a>
        format.
      </p>
      <p>
        If your primary sequences have FASTA headers following
        the UCSC style ("chromosome":"start"-"end"),
        and the chromosome names are in UCSC (not ENSEMBL) format,
        the GFF3 output will be suitable for uploading to the UCSC Genome Browser
        as a custom track.
      </p>
      `);
    case 'xstreme-results-tsv':
      return(`
        <p>
	  XSTREME outputs a tab-separated values (TSV) file ('xstreme.tsv') that
          contains one line for each motif found to be significantly enriched,
          sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the TSV file in a comment
          line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr><th>field</th><th>name</th><th>contents</th><tr>
          <tr><td>1</td><td>RANK</td><td> ` + get_doc_text('xstreme', 'xstreme-rank') + `</td></tr>
          <tr><td>2</td><td>SEED_MOTIF</td><td> ` + get_doc_text('xstreme', 'xstreme-seed') + `</td></tr>
          <tr><td>3</td><td>CLUSTER</td><td> ` + get_doc_text('xstreme', 'xstreme-cluster') + `</td></tr>
          <tr><td>4</td><td>SOURCE</td><td> ` + get_doc_text('shared', 'motif-db', 'the motif.', 'of the program that found the <i>de novo</i> motif, or the name of') + `</td></tr>
          <tr><td>5</td><td>ID</td><td> ` + get_doc_text('shared', 'motif-id', '', 'in the motif discovery program output or ') + `</td></tr>
          <tr><td>6</td><td>ALT_ID</td><td> ` + get_doc_text('shared', 'motif-alt-id', '', 'in the motif discovery program output or ') + `</td></tr>
          <tr><td>7</td><td>CONSENSUS</td><td> ` + get_doc_text('shared', 'motif-cons') + `</td></tr>
          <tr><td>8</td><td>WIDTH</td><td> ` + get_doc_text('xstreme', 'xstreme-width') + `</td><tr>
          <tr><td>9</td><td>SITES</td><td> ` + get_doc_text('xstreme', 'xstreme-sites') + `</td><tr>
          <tr><td>10</td><td>SEA_PVALUE</td><td> ` + get_doc_text('xstreme', 'xstreme-sea_pvalue') + `</td></tr>
          <tr><td>11</td><td>EVALUE</td><td> ` + get_doc_text('xstreme', 'xstreme-significance', 'value of EVALUE_ACC is "0"') + `</td></tr>
          <tr><td>12</td><td>EVALUE_ACC</td><td> ` + get_doc_text('xstreme', 'xstreme-sig_acc') + `</td></tr>
          <tr><td>13</td><td>SIM_SOURCE</td><td> ` + get_doc_text('xstreme', 'xstreme-sim_motif_source') + `</td></tr>
          <tr><td>14</td><td>SIM_MOTIF</td><td> ` + get_doc_text('xstreme', 'xstreme-sim_motif') + `</td></tr>
          <tr><td>15</td><td>MOTIF_URL</td><td> ` + get_doc_text('xstreme', 'xstreme-motif_url') + `</td></tr>
        </table>
      `);
    case 'xstreme-nonredundant-motifs':
      return(`
        <p>
          XSTREME outputs a text file ('xstreme.txt') containing the 'seed' motif from
          each cluster of similar (redundant) motifs. `
          + get_doc_text('xstreme', 'xstreme-seed') + 
          `The motifs are in <a href="` + site_url + `/doc/meme-format.html">Minimal MEME Motif format</a>.
        </p>
      `);
    default:
      return("Error--Unrecognized xstreme doc_type: " + doc_type);
  }
} // get_xstreme_doc_text
