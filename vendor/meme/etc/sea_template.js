pre_load_setup();
var sea_alphabet = new Alphabet(data.alphabet, data.background['frequencies']);

/*
 * name_from_source
 *
 * Makes a file name more human friendly to read.
 */
function name_from_source(source) {
  if (source == "-") {
    return "-"
  }
  // assume source is a file name
  var file = source.replace(/^.*\/([^\/]+)$/,"$1");
  var noext = file.replace(/\.[^\.]+$/, "");
  return noext.replace(/_/g, " ");
}

/*
 * pre_load_setup
 *
 *  Sets up initial variables which may be
 *  required for the HTML document creation.
 */
function pre_load_setup() {
  var seq_db = data['sequence_db'];
  if (!seq_db['name']) seq_db['name'] = name_from_source(seq_db['source']);
  // get the names of the motif databases
  var dbs = data['motif_dbs'];
  var motif_count = 0;
  for (var i = 0; i < dbs.length; i++) {
    var db = dbs[i];
    if (!db['name']) db['name'] = name_from_source(db['source']);
    motif_count += db['count'];
  }
  dbs['count'] = motif_count;
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

function toExp(num, prec) {
  if (typeof(num) === "number") return num.toExponential(prec);
  return "";
}

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
  var seq_len = data["sequence_db"]["maxlen"];
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

/*
 * make_results_table
 *
 * Add the results to the results table.
 */
function make_results_table() {
  var i;
  var distr = (data['motifs'].length > 0 && data['motifs'][0]['site_distr']);
  var table = $("results");
  var tbody = document.createElement("tbody");
  var motif_dbs = data['motif_dbs'];
  var motifs = data['motifs'];
  for (var i = 0; i < motifs.length; i++) {
    var motif = motifs[i];
    var row = tbody.insertRow(tbody.rows.length);
    row.id = "db_" + motif['db'] + "_motif_" + motif['id'];
    row['motif'] = motif;
    var canvas = document.createElement('canvas');
    var pspm = new Pspm(motif['pwm'], motif['id'], 0, 0, motif['motif_nsites'], motif['motif_evalue']);
    var logo = logo_1(sea_alphabet, '', pspm);
    draw_logo_on_canvas(logo, canvas, false, 0.5);
    var pos = data['sequence_db']['count'] - data['sequence_db']['holdout'];
    var neg = data['control_db']['count'] - data['control_db']['holdout'];
    var tp = motif['tp'];
    var fp = motif['fp'];
    add_cell(row, canvas, 'col_logo');
    add_text_cell(row, motif_dbs[motif['db']]['name'], 'col_db');
    add_cell(row, make_link(motif['id'], motif['url']), 'col_id');
    add_text_cell(row, motif['alt'], 'col_name');
    add_text_cell(row, motif['pvalue'], 'col_pvalue');
    add_text_cell(row, motif['evalue'], 'col_evalue');
    add_text_cell(row, motif['qvalue'], 'col_qvalue');
    add_text_cell(row, tp.toFixed(0) + ' / ' + pos + ' (' + (100*tp/pos).toFixed(1) + "%)", 'col_tp');
    add_text_cell(row, fp.toFixed(0) + ' / ' + neg + ' (' + (100*fp/neg).toFixed(1) + "%)", 'col_fp');
    add_text_cell(row, motif['enr_ratio'].toFixed(2), 'col_ratio');
    add_text_cell(row, (motif['score_thresh']).toFixed(2), 'score_thresh');
    add_cell(row, make_distribution(motif), "col_distribution");
    add_cell(row, make_histogram(motif), "col_histogram");
  }
  table.appendChild(tbody);
}
