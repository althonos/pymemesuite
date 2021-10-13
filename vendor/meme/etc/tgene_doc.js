//
// tgene_doc.js
//

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_tgene_doc_text(doc_type, extra) {
  if (extra == undefined) {extra = ""};
  switch (doc_type) {
    case 'gene-id':
      return(`
	The ID of the gene.
      `);
    case 'gene-name':
      return(`
        The name of the gene. This will be "." if there was no name provided for this 
	gene in your annotation file.
      `);
    case 'tss-id':
      return(`
	The ID of the TSS (transcription start site) of the gene.
      `);
    case 'tss-locus':
      return(`
	The genomic coordinates of the TSS (transcription start site) of the gene.
      `);
    case 'strand':
      return(`
	The chromosomal strand on which the gene is located.
      `);
    case 'max-expr':
      return(`
	The maximum expression level of the TSS across the panel of tissues
	(or 0 if you did not provide a panel).
      `);
    case 're-locus':
      return(`
	The genomic coordinates of a (potential) regulatory element (RE) you 
	provided in your (RE) locus file.
      `);
    case 'max-hist':
      return(`
	The maximum histone level of the RE Locus across the panel of tissues
	(or 0 if you did not provide a panel).
      `);
    case 'distance':
      return(`
	The distance between the TSS and the RE_Locus, taking Strand into account, 
	so that negative distances mean the RE is upstream of the TSS.  The distance
	is measured from the TSS to the closest <i>edge</i> of the RE_Locus.
      `);
    case 'closest-locus':
      return(`
	Equal to 'T' if RE_Locus is the <b>closest locus</b> to the TSS among
	all the loci in your (RE) locus file; equal to 'F' otherwise.
      `);
    case 'closest-tss':
      return(`
	Equal to 'T' if TSS_Locus is the <b>closest TSS</b> to RE_Locus among
	all the TSS loci in your annotation file; equal to 'F' otherwise.  
      `);
    case 'histone':
      return(`
	The name of the histone modification used in calculating the
	Pearson correlation between expression of the TSS and
	the level of the histone modification at the RE_Locus
	(or blank if you did not provide a panel of tissues).
	The correlation is measured after log-transforming both
	variables: <tt>x_new = log(x+1)</tt>.
      `);
    case 'correlation':
      return(`
	The Pearson correlation between the expression of the 
	TSS and the level of the histone modification at the RE_Locus
	(or 0 if you did not provide a panel of tissues).
	The correlation is measured after log-transforming both
	The correlation is computed after log-transforming both
	variables: <tt>x_new = log(x+1)</tt>.
      `);
    case 'correlation-sign':
      return(`
	The sign of the Pearson correlation of the expression of the 
	TSS and the level of the histone modification at the RE_Locus.
	The correlation is computed after log-transforming both
	variables: <tt>x_new = log(x+1)</tt>.
      `);
    case 'corr-pvalue':
      return(`
	The <b>unadjusted</b> correlation <i>p</i>-value of the potential regulatory link between this TSS and RE_Locus
	(or 1.0 if you did not provide a panel of tissues).
	Specifically, this is the empirically estimated,
	two-tailed <i>p</i>-value of the Pearson correlation 
	between the expression of the TSS and the level of the 
	histone modification at the RE_Locus.<br>
        T-Gene estimates correlation <i>p</i>-values by generating a null distribution for
	the Pearson correlation based on multiple permutations
	of the order of the expression values for the TSS.
	Each permutation breaks the relationship between the expression
	values of this TSS and the histone values of this RE_Locus, 
	resulting in a 'random' value for the correlation.
      `);
    case 'dist-pvalue':
      return(`
	The <b>unadjusted</b> distance <i>p</i>-value of the potential regulatory link between this TSS and RE_Locus.
	The distance <i>p</i>-value is based on assuming a uniform distribution of 
	(absolute) distances in the range <br>
        &nbsp;&nbsp;<tt>[0, maximum_link_distance]</tt>.<br>
	It is defined as <br>
	&nbsp;&nbsp;<tt>(2*distance + RE_width) / (2*maximum_link_distance + RE_width)</tt>.
      `);
    case 'cnd-pvalue':
      return(`
	The <b>unadjusted</b> CnD <i>p</i>-value of the potential regulatory link between this TSS and RE_Locus.
        (CnD stands for "Correlation aNd Distance".)
	Specifically, this is the <i>p</i>-value of the product of the correlation <i>p</i>-value
	and the distance <i>p</i>-value.  If you did not provide a panel of tissues, the CnD <i>p</i>-value
	is equal to the Distance <i>p</i>-value.
      `);
    case 'qvalue':
      return(`
	The q-value of the potential regulatory link between this TSS and RE_Locus.
	The q-value is defined as the minimum false discovery rate (FDR) required
	to consider this link statistically significant.
	If you provided a panel of tissues, T-Gene computes the q-value 
        based on the CnD <i>p</i>-value, otherwise it bases it on the 
	Distance <i>p</i>-value.<br><br>
        T-Gene estimates q-values from all the link <i>p</i>-values using the 
        method proposed by Benjamini & Hochberg 
        (<i>Journal of the Royal Statistical Society B</i>, 57:289-300, 1995).
        See also Storey JD, Tibshirani R. Statistical significance for
        genome-wide studies, <i>Proc. Natl. Acad. Sci. USA</i> (2003) <b>100</b>:9440&ndash;9445.
      `);
    case 'links-table':
      var i = 1;
      return(`
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_ID</td> <td>` + get_doc_text('tgene', 'gene-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Gene_Name</td> <td>` + get_doc_text('tgene', 'gene-name') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_ID</td> <td>` + get_doc_text('tgene', 'tss-id') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>TSS_Locus</td> <td>` + get_doc_text('tgene', 'tss-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Strand</td> <td>` + get_doc_text('tgene', 'strand') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Max_Expr</td> <td>` + get_doc_text('tgene', 'max-expr') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>RE_Locus</td> <td>` + get_doc_text('tgene', 're-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Max_Hist</td> <td>` + get_doc_text('tgene', 'max-hist') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Distance</td> <td>` + get_doc_text('tgene', 'distance') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Closest_Locus</td> <td>` + get_doc_text('tgene', 'closest-locus') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Closest_TSS</td> <td>` + get_doc_text('tgene', 'closest-tss') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Histone</td> <td>` + get_doc_text('tgene', 'histone') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Correlation</td> <td>` + get_doc_text('tgene', 'correlation') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Correlation_P_Value</td> <td>` + get_doc_text('tgene', 'corr-pvalue') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Distance_P_Value</td> <td>` + get_doc_text('tgene', 'dist-pvalue') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>CnD_P_Value</td> <td>` + get_doc_text('tgene', 'cnd-pvalue') + `</td> </tr>
	  <tr> <td>` + i++ + `</td> <td>Q_Value</td> <td>` + get_doc_text('tgene', 'qvalue') + `</td> </tr>
      `);
    case 'html-file-short':
      return(`
	T-Gene outputs an HTML file that provides the results in a human-readable format; 
	this file allows interactive selection, filtering and sorting of the potential regulatory links. 
      `);
    case 'html-file':
      return(`<p>` + get_tgene_doc("", 'html-file-short') + `</p>`);
    case 'tsv-general':
      return(`
	The first line in the file contains the (tab-separated) names of the fields.
	Your command line and other program information is given at the end of the file in 
	comment lines starting with the character '#'.
	The names and meanings of each of the fields are described in the table below.
      `);
    case 'links-tsv-short':
      return(`
	T-Gene outputs a comprehensive tab-separated values (TSV) file ('links.tsv') 
	that contains one line for each potential regulatory link that was found.
	The lines are sorted in order of increasing (unadjusted) <i>p</i>-value.
      `);
    case 'links-tsv':
      return(`
	<p>` + 
        get_doc_text('tgene', 'links-tsv-short') + 
	get_doc_text('tgene', 'tsv-general') + 
	get_doc_text('tgene', 'links-table') + `
	</p>`
      );
    default:
      return("Error--Unrecognized tgene doc_type: " + doc_type);
  } // switch
} // get_tgene_doc_text
