//
// submit_or_download_motif.js
//

function make_submit_or_download_motif_form(id, site_url, program) {
  var html = `
    <div class="popup_wrapper">
      <div class="popup" style="display:none; top: -150px;" id="download">
	<div>
	  <div style="float:right; ">
	    <div id="outpop_close" class="close" tabindex="0">x</div>
	  </div>
	  <h2 class="mainh" style="margin:0; padding:0;">Submit or Download</h2>
	  <div style="clear:both"></div>
	</div>
	<div style="height:100px">
	  <div style="float:right; width: 30px;">
	    <div id="outpop_prev" class="navarrow" tabindex="0">
	      <span class="inactive">&#8679;</span><span class="active">&#11014;</span>
	    </div>
	    <div id="outpop_num" class="navnum"></div>
	    <div id="outpop_next" class="navarrow" tabindex="0">
	      <span class="inactive">&#8681;</span><span class="active">&#11015;</span>
	    </div>
	  </div>
	  <div id="logo_box" style="height: 100px; margin-right: 40px;">
	    <canvas id="outpop_logo" height="100" width="250"></canvas>
	    <canvas id="outpop_logo_rc" height="100" width="250"></canvas>
	  </div>
	</div>
	<!-- tabs start -->
	<div class="tabArea top">
	  <span id="outpop_tab_1" class="tab">Submit Motif</span><span
	    id="outpop_tab_2" class="tab middle">Download Motif</span><span
	    id="outpop_tab_3" class="tab middle">Download Logo</span>
	</div>
	<div class="tabMain top">
	  <!-- Submit to another program -->
	  <div id="outpop_pnl_1">
	    <h4 class="compact">Submit to program</h4>
	    <table id="programs" class="programs">
	      <tr>
		<td><input type="radio" name="program" value="tomtom" id="submit_tomtom"></td>
		<td><label for="submit_tomtom">Tomtom</label></td>
		<td><label for="submit_tomtom">Find similar motifs in
		    published libraries or a library you supply.</label></td>
	      </tr>
	      <tr>
		<td><input type="radio" name="program" value="fimo" id="submit_fimo"></td>
		<td><label for="submit_fimo">FIMO</label></td>
		<td><label for="submit_fimo">Find motif occurrences in
		    sequence data.</label></td>
	      </tr>
	      <tr>
		<td><input type="radio" name="program" value="mast" id="submit_mast"></td>
		<td><label for="submit_mast">MAST</label></td>
		<td><label for="submit_mast">Rank sequences by affinity to
		    groups of motifs.</label></td>
	      </tr>
	      <tr class="dna_only">
		<td><input type="radio" name="program" value="gomo" id="submit_gomo"></td>
		<td><label for="submit_gomo">GOMo</label></td>
		<td><label for="submit_gomo">Identify possible roles (Gene
		    Ontology terms) for motifs.</label></td>
	      </tr>
	      <tr class="dna_only">
		<td><input type="radio" name="program" value="spamo" id="submit_spamo"></td>
		<td><label for="submit_spamo">SpaMo</label></td>
		<td><label for="submit_spamo">Find other motifs that are
		    enriched at specific close spacings which might imply the existence of a complex.</label></td>
	      </tr>
	    </table>
	  </div>
	  <!-- download text format  -->
	  <div id="outpop_pnl_2">
	    <div>
	      <label for="text_format">Format:</label>
	      <select id="text_format">
		<option value="0">Count Matrix</option>
		<option value="1">Probability Matrix</option>
		<option value="2">Minimal MEME</option> ` + 
		(program == "MEME" ? `
		<option value="3">FASTA</option>
		<option value="4">Raw</option> ` : ``) + `
	      </select>
	    </div>
	    <textarea id="outpop_text" name="content"
	      style="width:99%; white-space: pre; word-wrap: normal; overflow-x: scroll;"
	      rows="8" readonly="readonly" wrap="off"></textarea>
	    <a id="outpop_text_dl" download="meme.txt" href=""></a>
	  </div>
	  <!-- download logo format -->
	  <div id="outpop_pnl_3">
	    <form id="logo_form" method="post" action="">
	      <input type="hidden" name="program" value=" ` + program + `"/>
	      <input type="hidden" id="logo_motifs" name="motifs" value=""/>
	      <input type="hidden" id="logo_id1" name="id1" value=""/>
	      <table>
		<tr>
		  <td><label for="logo_format">Format:</label></td>
		  <td>
		    <select id="logo_format" name="png">
		      <option value="1">PNG (for web)</option>
		      <option value="0">EPS (for publication)</option>
		    </select>
		  </td>
		</tr>
		<tr>
		  <td><label for="logo_rc">Orientation:</label></td>
		  <td>
		    <select id="logo_rc" name="rc1">
		      <option value="0">Normal</option>
		      <option value="1" id="logo_rc_option">Reverse Complement</option>
		    </select>
		  </td>
		</tr>
		<tr>
		  <td><label for="logo_ssc">Small Sample Correction:</label></td>
		  <td>
		    <input type="hidden" id="logo_err" name="errbars" value="0"/>
		    <select id="logo_ssc" name="ssc">
		      <option value="0">Off</option>
		      <option value="1">On</option>
		    </select>
		  </td>
		</tr>
		<tr>
		  <td><label for="logo_width">Width:</label></td>
		  <td>
		    <input type="text" id="logo_width" size="4" placeholder="default" name="width"/>&nbsp;cm
		  </td>
		</tr>
		<tr>
		  <td><label for="logo_height">Height:</label></td>
		  <td>
		    <input type="text" id="logo_height" size="4" placeholder="default" name="height"/>&nbsp;cm
		  </td>
		</tr>
	      </table>
	    </form>
	  </div>
	  <!-- Buttons -->
	  <div>
	    <div style="float:left;">
	      <input type="button" id="outpop_do" value="Submit" />
	    </div>
	    <div style="float:right;">
	      <input id="outpop_cancel" type="button" value="Cancel" />
	    </div>
	    <div style="clear:both;"></div>
	  </div>
	</div>
      </div>
    </div>
  `;
  document.getElementById(id).insertAdjacentHTML('beforeend', html);
  $("logo_form").action = site_url + "/utilities/generate_logo";
} // make_submit_or_download_motif_form

//
// Functions to update the submit_or_download_motif form.
//

//
// Initialise and display the download popup.
//
function action_show_outpop(e, ordinal) {
  "use strict";
  function init() {
    "use strict";
    var close_btn, next_btn, prev_btn, cancel_btn, do_btn;
    var tab1, tab2, tab3;
    var pnl1, pnl2, pnl3;
    var format_list;
    var tbl_submit, inputs, i, default_prog;
    close_btn = $("outpop_close");
    close_btn.addEventListener("click", action_hide_outpop, false);
    close_btn.addEventListener("keydown", action_hide_outpop, false);
    next_btn = $("outpop_next");
    next_btn.addEventListener("click", action_outpop_next, false);
    next_btn.addEventListener("keydown", action_outpop_next, false);
    prev_btn = $("outpop_prev");
    prev_btn.addEventListener("click", action_outpop_prev, false);
    prev_btn.addEventListener("keydown", action_outpop_prev, false);
    cancel_btn = $("outpop_cancel");
    cancel_btn.addEventListener("click", action_hide_outpop, false);
    do_btn = $("outpop_do");
    do_btn.addEventListener("click", action_outpop_submit, false);
    tab1 = $("outpop_tab_1");
    tab1.tabIndex = 0;
    tab1.addEventListener("click", action_outpop_tab, false);
    tab1.addEventListener("keydown", action_outpop_tab, false);
    tab2 = $("outpop_tab_2");
    tab2.tabIndex = 0;
    tab2.addEventListener("click", action_outpop_tab, false);
    tab2.addEventListener("keydown", action_outpop_tab, false);
    tab3 = $("outpop_tab_3");
    tab3.tabIndex = 0;
    tab3.addEventListener("click", action_outpop_tab, false);
    tab3.addEventListener("keydown", action_outpop_tab, false);
    pnl1 = $("outpop_pnl_1");
    pnl2 = $("outpop_pnl_2");
    pnl3 = $("outpop_pnl_3");
    toggle_class(tab1, "activeTab", true);
    toggle_class(tab2, "activeTab", false);
    toggle_class(tab3, "activeTab", false);
    pnl1.style.display = "block";
    pnl2.style.display = "none";
    pnl3.style.display = "none";
    format_list = $("text_format");
    format_list.addEventListener("change", action_outpop_format, false);
    // setup program selection
    tbl_submit = $("programs");
    // when not dna, hide the inputs for programs that require dna motifs
    toggle_class(tbl_submit, "alphabet_dna", current_alphabet.has_complement());//TODO FIXME alphabet_dna is a bad name for a field when allowing custom alphabets
    // add a click listener for the radio buttons
    inputs = tbl_submit.querySelectorAll("input[type='radio']");
    for (i = 0; i < inputs.length; i++) {
      inputs[i].addEventListener("click", action_outpop_program, false);
    }
    // ensure that a default program option is selected for DNA and Protein
    default_prog = document.getElementById(current_alphabet.has_complement() ? "submit_tomtom" : "submit_fimo"); //TODO FIXME Tomtom might require a more strict definition of DNA
    default_prog.checked = true;
    action_outpop_program.call(default_prog);
    // disable reverse-complement when not DNA
    $("logo_rc_option").disabled = !current_alphabet.has_complement();
    // set errorbars on when ssc is on
    $("logo_ssc").addEventListener("change", action_outpop_ssc, false);
  }
  // store the focused element
  action_hide_outpop.last_active = document.activeElement;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  // hide the help popup
  help_popup();
  // on first load initilize the popup
  if (!action_show_outpop.ready) {
    init();
    action_show_outpop.ready = true;
  }
  update_outpop_motif(ordinal - 1);
  // display the download popup
  $("grey_out_page").style.display = "block";
  $("download").style.display = "block";
  $("outpop_close").focus();
} // action_show_output

//
// Hide the submit or download popup.
//
function action_hide_outpop(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  $("download").style.display = "none";
  $("grey_out_page").style.display = "none";
  if (typeof action_hide_outpop.last_active !== "undefined") {
    action_hide_outpop.last_active.focus();
  }
} // action_hide_outpop

/*
 * Show the next motif in the download popup.
 */
function action_outpop_next(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  update_outpop_motif(current_motif + 1);
} // action_outpop_next

/*
 * Show the previous motif in the download popup.
 */
function action_outpop_prev(e) {
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  update_outpop_motif(current_motif - 1);
} // action_outpop_prev

/*
 * Highlight the selected row in the program list.
 */
function action_outpop_program() {
  "use strict";
  var table, tr, rows, i;
  tr = find_parent_tag(this, "TR");
  table = find_parent_tag(tr, "TABLE");
  rows = table.querySelectorAll("tr");
  for (i = 0; i < rows.length; i++) {
    toggle_class(rows[i], "selected", rows[i] === tr);
  }
} // action_outpop_program

/*
 * Submit the motif to the selected program.
 */
function action_outpop_submit(e) {
  "use strict";
  var form, input, program, motifs;
  // find out which program is selected
  var radios, i;
  radios = document.getElementsByName("program");
  program = "fimo"; // default to fimo, since it works with all alphabet types
  for (i = 0; i < radios.length; i++) {
    if (radios[i].checked) program = radios[i].value;
  }

  motifs = motif_minimal_meme(data.motifs[current_motif]);
  form = document.createElement("form");
  form.setAttribute("method", "post");
  form.setAttribute("action", site_url + "/tools/" + program);

  input = document.createElement("input");
  input.setAttribute("type", "hidden");
  input.setAttribute("name", "motifs_embed");
  input.setAttribute("value", motifs);
  form.appendChild(input);

  document.body.appendChild(form);
  form.submit();
  document.body.removeChild(form);
} // action_outpop_submit(e)

/*
 * Enable error bars when small sample correction is enabled.
 */
function action_outpop_ssc() {
  "use strict";
  $("logo_err").value = $("logo_ssc").value;
} // action_outpop_ssc

//
// Update the motif logos and format download text in the popup.
// This is called whenever the current motif changes.
// 
function update_outpop_motif(index) {
  "use strict";
  var motifs, motif, pspm, logo, canvas, num;
  motifs = data["motifs"];
  if (index < 0 || index >= motifs.length) {return;}
  current_motif = index;
  motif = motifs[index];
  pspm = new Pspm(motif["pwm"]);
  logo = new Logo(current_alphabet, "");
  logo.add_pspm(pspm, 0);
  canvas = $("outpop_logo");
  canvas.width = canvas.width; // clear canvas
  draw_logo_on_canvas(logo, canvas, false);
  canvas = $("outpop_logo_rc");
  canvas.width = canvas.width; // clear rc canvas
  if (data.options.strands === "both" || data.options.revcomp) {
    pspm.reverse_complement(current_alphabet);
    logo = new Logo(current_alphabet, "");
    logo.add_pspm(pspm, 0);
    draw_logo_on_canvas(logo, canvas, false);
  }
  num = $("outpop_num");
  num.innerHTML = "";
  num.appendChild(document.createTextNode("" + (index + 1)));
  update_outpop_format(index);
} // action_outpop_motif

//
// Create the download menu.
//
function update_outpop_format(index) {
  var motif = data.motifs[index];
  var fn = [motif_count_matrix, motif_prob_matrix, motif_minimal_meme, motif_fasta, motif_raw];
  var suffix = ["_counts.txt", "_freqs.txt", ".meme", "_fasta.txt", "_raw.txt"];
  var format = parseInt($("text_format").value);
  var text = fn[format](motif);
  prepare_download(text, "text/plain", motif.id + suffix[format], $("outpop_text_dl"));
  $("outpop_text").value = text;
} // update_outpop_format

/*
 * Update the text in the download format popup.
 */
function action_outpop_format() {
  update_outpop_format(current_motif);
} // action_outpop_format

/*
 * Download the format text.
 * Wire the link containing the data URI text to a download button so it looks
 * the same as the server submit stuff.
 */
function action_outpop_download_motif(e) {
  $("outpop_text_dl").click();
} // action_outpop_download_motif

/*
 * Download the motif logo.
 * The EPS format can be calculated locally in Javascript
 */
function action_outpop_download_logo(e) {
  "use strict";
  var motif = data.motifs[current_motif];
  if ($("logo_format").value === "0") { // EPS
    var pspm, logo, eps;
    var logo_rc, logo_ssc, logo_width, logo_height;
    logo_rc = ($("logo_rc").value === "1");
    logo_ssc = ($("logo_ssc").value === "1");
    logo_width = parseFloat($("logo_width").value);
    if (isNaN(logo_width) || !isFinite(logo_width) || logo_width <= 0) logo_width = null;
    logo_height = parseFloat($("logo_height").value);
    if (isNaN(logo_height) || !isFinite(logo_height) || logo_height <= 0) logo_height = null;
    // create a PSPM from the motif
    pspm = motif_pspm(motif);
    if (logo_rc) pspm.reverse_complement(current_alphabet);
    logo = new Logo(current_alphabet);
    logo.add_pspm(pspm, 0);
    eps = logo.as_eps({"ssc": logo_ssc, "logo_width": logo_width, "logo_height": logo_height});
    prepare_download(eps, "application/postscript", motif.id + (logo_rc ? "_rc" : "") + ".eps");
  } else {
    $("logo_motifs").value = motif_minimal_meme(motif);
    $("logo_id1").value = motif.id;
    $("logo_form").submit();
  }
} // action_outpop_download_logo

/*
 * Change the selected tab in the download popup.
 */
function action_outpop_tab(e) {
  "use strict";
  var tab1, tab2, tab3, pnl1, pnl2, pnl3, do_btn;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  tab1 = $("outpop_tab_1");
  tab2 = $("outpop_tab_2");
  tab3 = $("outpop_tab_3");
  pnl1 = $("outpop_pnl_1");
  pnl2 = $("outpop_pnl_2");
  pnl3 = $("outpop_pnl_3");
  do_btn = $("outpop_do");

  toggle_class(tab1, "activeTab", (this === tab1));
  toggle_class(tab2, "activeTab", (this === tab2));
  toggle_class(tab3, "activeTab", (this === tab3));
  pnl1.style.display = ((this === tab1) ? "block" : "none");
  pnl2.style.display = ((this === tab2) ? "block" : "none");
  pnl3.style.display = ((this === tab3) ? "block" : "none");
  do_btn.value = ((this === tab1) ? "Submit" : "Download");
  do_btn.removeEventListener("click", action_outpop_submit, false);
  do_btn.removeEventListener("click", action_outpop_download_logo, false);
  do_btn.removeEventListener("click", action_outpop_download_motif, false);
  if (this === tab1) {
    do_btn.addEventListener("click", action_outpop_submit, false);
  } else if (this === tab2) {
    do_btn.addEventListener("click", action_outpop_download_motif, false);
  } else {
    do_btn.addEventListener("click", action_outpop_download_logo, false);
  }
} // action_outpop_tab

function motif_fasta(motif) {
  "use strict";
  var sites, site, seq, sequences, sequence, i, num, counter, out;
  counter = {};
  sequences = data["sequence_db"]["sequences"];
  sites = motif["sites"];
  out = "";
  for (i = 0; i < sites.length; i++) {
    site = sites[i];
    seq = site["seq"];
    sequence = sequences[seq];
    counter[seq] = (num = counter[seq]) ? (++num) : (num = 1); // inc counter
    if (i !== 0) {out += "\n";}
    out += ">" + sequence["name"] + "_site_" + num + " offset= " + site["pos"] +
      (site["rc"] ? " RC\n" : "\n");
    out += site["match"];
  }
  return out;
} // motif_fasta

function motif_raw(motif) {
  "use strict";
  var sites, i, out;
  sites = motif["sites"];
  out = "";
  for (i = 0; i < sites.length; i++) {
    if (i !== 0) {out += "\n";}
    out += sites[i]["match"];
  }
  return out;
} // motif_raw

/*
 * Create a pspm for the given motif data
 */
function motif_pspm(motif) {
  var pwm = motif.pwm;
  var name = motif.id;
  var ltrim = 0;
  var rtrim = 0;
  var nsites = motif.nsites;
  var sig = (current_program === "STREME" ? motif.test_pvalue : motif.evalue);
  return new Pspm(pwm, name, ltrim, rtrim, nsites, sig, null, motif.alt, current_program);
} // motif_pspm

/*
 * Create a count matrix from the given motif data
 */
function motif_count_matrix(motif) {
  return motif_pspm(motif).as_count_matrix();
} // motif_count_matrix

/*
 * Create a probablity matrix from the given motif data
 */
function motif_prob_matrix(motif) {
  return motif_pspm(motif).as_probability_matrix();
} // motif_prob_matrix

function motif_minimal_meme(motif) {
  var strands;
  if (current_program === "STREME") {
    strands = (data.options.strands === "both" ? 2 : 1);
  } else {
    strands = (data.options.revcomp ? 2 : 1);
  }
  return motif_pspm(motif).as_meme({
    "with_header": true,
    "with_pspm": true,
    "with_pssm": (current_program === "MEME" ? true : false),
    "version": data["version"],
    "alphabet": current_alphabet,
    "strands": strands
  });
}
