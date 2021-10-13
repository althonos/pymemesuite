var current_program = "STREME";
var current_alphabet = new Alphabet(data.alphabet, data.background.freqs);
var current_motif = 0;

/*
 * Make the table header for the discovered motifs.
 */
function make_motif_header() {
  var is_test_set;
  is_test_set = data.test_positives.count > 0;

  var row = document.createElement("tr");
  add_text_header_cell(row, "Motif", "pop_motifs_word", "motif_word");
  add_text_header_cell(row, "Logo", "pop_motifs_logo", "motif_logo");
  if (data.options.strands === "both") {
    add_text_header_cell(row, "RC Logo", "pop_motifs_rc_logo", "motif_logo");
  }
  if (is_test_set) {
    add_text_header_cell(row, "P-value", "pop_motif_test_pvalue", "motif_pvalue");
    add_text_header_cell(row, "E-value", "pop_motif_test_evalue", "motif_evalue", "", "", 1);
  } else {
    add_text_header_cell(row, "Score", "pop_motif_train_pvalue", "motif_score");
  }
  add_text_header_cell(row, "Sites", "pop_motif_sites", "motif_sites");
  add_text_header_cell(row, "More", "pop_more", "motif_more");
  add_text_header_cell(row, "Submit/Download", "pop_submit_dl", "motif_submit");
  row.className = "more";
  add_text_header_cell(row, "Positional Distribution", "pop_site_distr", "col_distribution", "", "", 1);
  add_text_header_cell(row, "Matches per Sequence", "pop_site_hist", "col_histogram", "", "", 1);
  return row;
} // make_motif_header

/*
 * Make a compact motif summary row for the discovered motifs.
 */
function make_motif_row(tbody, ordinal, motif) {
  is_test_set = data.test_positives.count > 0;

  var row = document.createElement("tr");
  row.id = motif["id"];
  add_text_cell(row, motif["id"], "motif_word");
  var pspm = new Pspm(motif["pwm"]);
  add_cell(row, make_logo(current_alphabet, pspm, 50, false, 0, "normal_logo"), "motif_logo");
  if (data.options.strands === "both") {
    add_cell(row, make_logo(current_alphabet, pspm, 50, true, 0, "flipped_logo"), "motif_logo");
  }
  if (is_test_set) {
    add_text_cell(row, motif["test_pvalue"], "motif_pvalue");
    add_text_cell(row, motif["test_evalue"], "motif_evalue");
  } else {
    add_text_cell(row, motif["train_pvalue"], "motif_score");
  }
  var tp = motif["train_pos_count"] + motif["test_pos_count"];
  var tp_percent = (100 * tp / (data["train_positives"]["count"] + data["test_positives"]["count"])).toFixed(1);
  add_text_cell(row, tp + " (" + tp_percent + "%)", "motif_sites");
  add_cell(row, make_sym_btn(text_pair("\u21A7", "less", "\u21A5", "more"), "Show more information.", function(e) { toggle_class(tbody, "collapsed"); }, "\u21A5", ""), "motif_more");
  add_cell(row, make_sym_btn("\u21E2", "Submit the motif to another MEME Suite program or download it.", function(e) { action_show_outpop(e, ordinal); }), "motif_submit");
  // Site Distribution Plot
  add_cell(row, make_distribution(motif), "col_distribution");
  add_cell(row, make_histogram(motif), "col_histogram");
  return row;
} // make_motif_row

/*
 * Make an expanded view of a discovered motif.
 */
function make_motif_exp(tbody, ordinal, motif) {
  "use strict";

  var box = $("tmpl_motif_expanded").cloneNode(true);
  toggle_class(box, "template", false);
  box.id = "";
  var pspm = new Pspm(motif["pwm"]);
  find_child(box, "tvar_logo").appendChild(make_logo(current_alphabet, pspm, 150, false, 0, "normal_logo"));
  if (data.options.strands === "both") {
    find_child(box, "tvar_rclogo").appendChild(make_logo(current_alphabet, pspm, 150, true, 0, "flipped_logo"));
  }

  set_tvar(box, "tvar_train_p", motif["train_pos_count"]);
  set_tvar(box, "tvar_train_p_total", data.train_positives.count);
  set_tvar(box, "tvar_train_pos_ratio", 
    (data.train_positives.count === 0) ? 0+"%" : (100 * motif["train_pos_count"] / data.train_positives.count).toFixed(1) + "%");
  set_tvar(box, "tvar_train_n", motif["train_neg_count"]);
  set_tvar(box, "tvar_train_n_total", data.train_negatives.count);
  set_tvar(box, "tvar_train_neg_ratio", 
    (data.train_negatives.count === 0) ? 0+"%" : (100 * motif["train_neg_count"] / data.train_negatives.count).toFixed(1) + "%");
  set_tvar(box, "tvar_train_p_cd", motif["train_pos_count"]);
  set_tvar(box, "tvar_train_dtc", motif["train_dtc"]);
  set_tvar(box, "tvar_train_pvalue", motif["train_pvalue"]);

  set_tvar(box, "tvar_test_p", motif["test_pos_count"]);
  set_tvar(box, "tvar_test_p_total", data.test_positives.count);
  set_tvar(box, "tvar_test_pos_ratio", 
    (data.test_positives.count === 0) ? 0+"%" : (100 * motif["test_pos_count"] / data.test_positives.count).toFixed(1) + "%");
  set_tvar(box, "tvar_test_n", motif["test_neg_count"]);
  set_tvar(box, "tvar_test_n_total", data.test_negatives.count);
  set_tvar(box, "tvar_test_neg_ratio", 
    (data.test_negatives.count === 0) ? 0+"%" : (100 * motif["test_neg_count"] / data.test_negatives.count).toFixed(1) + "%");
  set_tvar(box, "tvar_test_p_cd", motif["test_pos_count"]);
  set_tvar(box, "tvar_test_dtc", motif["test_dtc"]);
  set_tvar(box, "tvar_test_pvalue", motif["test_pvalue"]);
  set_tvar(box, "tvar_match_threshold", motif["score_threshold"]);

  var cell = document.createElement("td");
  cell.colSpan = 7;
  cell.appendChild(box);
  var row = document.createElement("tr");
  row.className = "more";
  row.appendChild(cell);

  return row;
} // make_motif_exp

/*
 * Create a table of discovered motifs. A fresh table body is used for each
 * motif to make hiding/showing rows with css easier.
 */
function make_motifs() {
  "use strict";
  var i, row, tbody, motif, ordinal;
  // make the motifs table
  var container = $("motifs");
  container.innerHTML = ""; // clear content
  var table = document.createElement("table");
  // add a header that is always shown
  var thead = document.createElement("thead");
  thead.appendChild(make_motif_header());
  table.appendChild(thead);
  for (i = 0; i < data.motifs.length; i++) {
    ordinal = i + 1;
    motif = data.motifs[i];
    tbody = document.createElement("tbody");
    tbody.className = "collapsed";
    tbody.appendChild(make_motif_row(tbody, ordinal, motif));
    tbody.appendChild(make_motif_exp(tbody, ordinal, motif));
    // create a following header for every row except the last one
    if ((i + 1) < data.motifs.length) tbody.appendChild(make_motif_header());
    table.appendChild(tbody);
  }

  // Add stopping reason.
  var tfoot = document.createElement("tfoot");
  table.appendChild(tfoot);
  row = tfoot.insertRow(tfoot.rows.length);
  add_text_header_cell(row, data.stop_reason, "", "stop_reason", "", 6);
  row = tfoot.insertRow(tfoot.rows.length);
  add_text_header_cell(row, "STREME ran for " + data.runtime.cpu + " seconds.", "", "run_time", "", 6);

  container.appendChild(table);
  // Hide unneeded columns.
  var no_test_seq_cols = data.test_positives.count == 0;
  var no_cd_cols = data.options.objfun !== "Central Distance";
  var no_de_cols = data.options.objfun !== "Differential Enrichment";
  var ids, i;
  ids = no_test_seq_cols ? document.getElementsByClassName("test_seq_col") : [];
  for (i=0; i<ids.length; i++) { toggle_class(ids[i], "hide_col", true); }
  ids = no_cd_cols ? document.getElementsByClassName("cd_col") : [];
  for (i=0; i<ids.length; i++) { toggle_class(ids[i], "hide_col", true); }
  ids = no_de_cols ? document.getElementsByClassName("de_col") : [];
  for (i=0; i<ids.length; i++) { toggle_class(ids[i], "hide_col", true); }
} // make_motifs

var DelayLogoTask = function(logo, canvas) {
  this.logo = logo;
  this.canvas = canvas;
};

DelayLogoTask.prototype.run = function () {
  draw_logo_on_canvas(this.logo, this.canvas, false);
};

var DelayCentrimoPlotTask = function(graph, canvas) {
  this.graph = graph;
  this.canvas = canvas;
};

DelayCentrimoPlotTask.prototype.run = function () {
  this.graph.draw_lines(this.canvas.getContext("2d"),
      this.canvas.width, this.canvas.height);
};

function make_distribution(motif) {
  var canvas = document.createElement("canvas");
  canvas.title = "Distribution of best motif sites in sequences.";
  canvas.width = 200;
  canvas.height = 90;
  var color = "blue";
  var seq_len = data["train_positives"]["maxlen"];
  //var smooth = Math.ceil(seq_len*0.05);
  var rset = new CentrimoRSet(seq_len);
  if (motif["site_distr"]) {
    rset.add("", "", color, motif["len"], motif["total_sites"], motif["site_distr"]);
    var graph = new CentrimoGraph(rset, triangular_weights(20), null, null, null, null, data["options"]["align"], motif["len"]);
    add_draw_task(canvas, new DelayCentrimoPlotTask(graph, canvas));
  }
  return canvas;
} // make_distribution

function make_histogram(motif) {
  var canvas = document.createElement("canvas");
  var ctx = canvas.getContext('2d');
  canvas.width = 200;
  canvas.height = 115;
  var x_label = "Matches per Sequence";
  var y_label = "% Sequences";
  var color = "blue";
  var graph = new SiteHistogramGraph(motif["site_hist"]);
  graph.draw_graph(ctx, canvas.width, canvas.height, x_label, y_label, color);
  return canvas;
} // make_histogram
