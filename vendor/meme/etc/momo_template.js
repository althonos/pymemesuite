/*
 * pre_load_setup
 *
 *  Sets up initial variables which may be
 *  required for the HTML document creation.
 */
function pre_load_setup() {
}

/*
 * page_loaded
 *
 * Called when the page has loaded for the first time.
 */
function page_loaded() {
  post_load_setup();
}

/*
 * page_loaded
 *
 * Called when a cached page is reshown.
 */
function page_shown(e) {
  if (e.persisted) post_load_setup();
}

/*
 * post_load_setup
 *
 * Setup state that is dependant on everything having been loaded already.
 */
function post_load_setup() {
  update_scroll_pad();
}

/*
 * add_cell
 *
 * Add a cell to the table row.
 */
function add_cell(row, node, cls) {
  var cell = row.insertCell(row.cells.length);
  if (node) cell.appendChild(node);
  if (cls) cell.className = cls;
}

/*
 * add_text_cell
 *
 * Add a text cell to the table row.
 */
function add_text_cell(row, text, cls) {
  var node = null;
  if (text) node = document.createTextNode(text);
  add_cell(row, node, cls);
}

/*
 * add_img_cell
 *
 * Add an IMG cell to the table row.
 */
function add_img_cell(row, src, alt, height, width, cls) {
  var node = null;
  if (src) {
    node = document.createElement("img");
    node.setAttribute("src", src)
    if (alt) { node.setAttribute("alt", alt) }
    if (height > 0) { node.setAttribute("height", height) }
    if (width > 0) { node.setAttribute("width", width) }
  }
  add_cell(row, node, cls);
}

/*
 * add_textarea_cell
 *
 * Add an textarea cell to the table row.
 */
function add_textarea_cell(row, rows, cols, text, readonly, cls) {
  var node = null;
  node = document.createElement("textarea");
  if (rows>0) node.setAttribute("rows", rows);
  if (cols>0) node.setAttribute("cols", cols);
  if (readonly) node.setAttribute("readonly", true);
  //node.setAttribute("data-role", "none");
  node.value = text;
  add_cell(row, node, cls);
} // add_textarea_cell

function toExp(num, prec) {
  if (typeof(num) === "number") return num.toExponential(prec);
  return "";
}

function make_other_settings() {
  $("opt_filter").textContent = data.options.filter ? 
    "'" + data.options.filter.filter_field + "' " + 
    data.options.filter.filter_type + " " + 
    data.options.filter.filter_threshold 
    : "no";
  $("opt_remove_unknowns").textContent = data.options.remove_unknowns ? "yes" : "no";
  $("opt_eliminate_repeats").textContent = data.options.eliminate_repeat_width != 0 ? 
    "yes, if their " + data.options.eliminate_repeat_width + " central residues match" : "no";
  $("opt_single_motif_per_mass").textContent = data.options.single_motif_per_mass ? "yes" : "no";
  $("opt_seed").textContent = data.options.seed;
} // make_other_settings

/*
 * make_results_table
 *
 * Add the results to the results table.
 */
function make_results_table() {
  var width = data['options']['width'];
  var nlines = 9;

  // clear the table and add the items
  var tbl = $("results");

  if (data['options']['algorithm'] == "simple") {
    tbl.className += " no_score no_bg no_fold no_unadjusted_p no_n_tests no_adjusted_p no_modl_log";
  } else if (data['options']['algorithm'] == "MoDL") {
    tbl.className += " no_n_tests no_adjusted_p";
  } else if (data['options']['algorithm'] == "motif-x") {
    tbl.className += " no_modl_log";
    if (! data['pvalues_accurate']) {
      tbl.className += " no_n_tests no_adjusted_p";
    }
  }

  var tbody = tbl.tBodies[0];
  while (tbody.rows.length > 0) {
    tbody.deleteRow(0);
  }

  // add the rows to the table
  var mods = data['mods'];
  for (var i = 0; i < mods.length; i++) {
    var mod = mods[i];

    var motifs = mod['motifs'];
    for (var j = 0; j < motifs.length; j++) {
      var motif = motifs[j];

      // Get text containing the motif occurrences.
      var occurrences = motif['occurrences'];
      var occ_text = "";
      var n_occ = occurrences.length;
      for (var k = 0; k < n_occ; k++) {
	occ_text += occurrences[k] + "\n";
      }

      // Get the text containing the Log.
      var modl_log_text = "";
      var modl_log_text_width = 0;
      var n_steps = 0;
      if (j == 0 && data['options']['algorithm'] == "MoDL") {
	var modl_log = mod['modl_log'];
	var modl_steps = modl_log['modl_steps'];
        var text = "";
	for (var k = 0; k < modl_steps.length; k++) {
          text = "";
	  var step = modl_steps[k];
	  if (k == modl_log['final_step']) text += "*";
	  text += "STEP: " + step['step'] + ", DL: " + step['score'] + "\n";
	  n_steps += 1;
          if (text.length > modl_log_text_width) modl_log_text_width = text.length;
          modl_log_text += text;
	  var reg_exps = step['reg_exps'];
	  for (var m = 0; m < reg_exps.length; m++) {
	    text = "  " + reg_exps[m] + "\n";
	    n_steps += 1;
            if (text.length > modl_log_text_width) modl_log_text_width = text.length;
            modl_log_text += text;
	  }
	}

        modl_log_text += "\n";
	n_steps += 1;
        text = "Final Step: " + modl_log['final_step'] + "\n";
	n_steps += 1;
        if (text.length > modl_log_text_width) modl_log_text_width = text.length;
        modl_log_text += text;
        text = "Final DL: " + modl_log['final_dl'].toFixed(2) + "\n";
	n_steps += 1;
        if (text.length > modl_log_text_width) modl_log_text_width = text.length;
        modl_log_text += text;
        text = "Decrease: " + modl_log['decrease'].toFixed(2) + "\n";
	n_steps += 1;
        if (text.length > modl_log_text_width) modl_log_text_width = text.length;
        modl_log_text += text;
      } // modl_log

      var row = tbody.insertRow(tbody.rows.length);
      row.id = "mod_" + mod['mod_name'] + "_motif_" + motif['motif_name'];
      row['motif'] = motif;
      add_img_cell(row, motif['motif_png'], 'sequence logo of motif', 175, 0, 'col_logo');
      add_text_cell(row, mod['mod_name'], 'col_mod');
      add_text_cell(row, motif['motif_name'], 'col_motif');
      add_text_cell(row, motif['motif_regexp'], 'col_regexp');
      add_text_cell(row, motif['score'].toFixed(2), 'col_score');
      if (data['options']['algorithm'] == "simple") {
	add_text_cell(row, motif['fg_matches'], 'col_fg');
      } else {
	add_text_cell(row, motif['fg_matches'] + "/" + motif['fg_size'], 'col_fg');
      }
      add_text_cell(row, motif['bg_matches'] + "/" + motif['bg_size'], 'col_bg');
      add_text_cell(row, motif['fold'].toFixed(1), 'col_fold');
      add_text_cell(row, motif['m1'].toFixed(1) + "e" + motif['e1'], 'col_unadjusted_p');
      add_text_cell(row, motif['n_tests'], 'col_n_tests');
      add_text_cell(row, motif['m2'].toFixed(1) + "e" + motif['e2'], 'col_adjusted_p');
      add_textarea_cell(row, n_occ < nlines ? n_occ : nlines, width+1, occ_text, true, 'col_motif_occurrences');
      // Only show the MoDL log for first motif with a given mod.
      if (j==0) add_textarea_cell(row, n_steps < nlines ? n_steps : nlines, modl_log_text_width, modl_log_text, true, 'col_modl_log');
    }
  }
} // make_results_table
