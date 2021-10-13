var have_tissue_panel = (data.options["tissues"] != "");

var sort_table = {
  "regulatory_link": [
    {"name": "Closest TSS (CT)", "fn": sort_closest_tss, "priority": 1},
    {"name": "CnD p-value", "fn": sort_cnd_pvalue, "show": have_tissue_panel},
    {"name": "Correlation p-value", "fn": sort_corr_pvalue, "show": have_tissue_panel},
    {"name": "Distance p-value", "fn": sort_dist_pvalue},
    {"name": "Gene ID", "fn": sort_gene_id},
    {"name": "Gene Name", "fn": sort_gene_name},
    {"name": "TSS ID", "fn": sort_tss_id},
    {"name": "RE Locus", "fn": sort_re_locus},
    {"name": "Distance", "fn": sort_distance},
    {"name": "Absolute Distance", "fn": sort_absolute_distance},
    {"name": "Correlation (decreasing)", "fn": sort_correlation_decr, "show": have_tissue_panel},
    {"name": "Correlation (increasing)", "fn": sort_correlation_incr, "show": have_tissue_panel}
  ]
} ;

pre_load_setup();

function pre_load_setup() {
} //pre_load_setup

/*
 * page_loaded
 *
 * Called when the page has loaded for the first time.
 */
function page_loaded() {
  first_load_setup();
  post_load_setup();
} // page_loaded

/*
 * page_loaded
 *
 * Called when a cached page is reshown.
 */
function page_shown(e) {
  if (e.persisted) post_load_setup();
} // page_shown

/*
 * first_load_setup
 *
 * Setup state that is dependent on everything having been loaded already.
 * On browsers which cache state this is only run once.
 */
function first_load_setup() {
  // Turn off these columns by default.
  $('show_tss_locus').checked = false;
  $('show_strand').checked = false;
  $('show_max_expr').checked = false;
  $('show_max_hist').checked = false;
  $('show_closest_locus').checked = false;
  //$('show_closest_tss').checked = false;
  $('show_histone').checked = false;
  $('show_correlation_sign').checked = false;
  $('show_corr_pvalue').checked = false;
  if (have_tissue_panel) {
    $('show_dist_pvalue').checked = false;
  } else {
    //$('show_max_expr').checked = false;
    //$('show_max_hist').checked = false;
    //$('show_histone').checked = false;
    $('show_correlation').checked = false;
    //$('show_correlation_sign').checked = false;
    //$('show_corr_pvalue').checked = false;
    $('show_cnd_pvalue').checked = false;
  }
} // first_load_setup

/*
 * post_load_setup
 *
 * Setup state that is dependent on everything having been loaded already.
 */
function post_load_setup() {
  "use strict";
  var tbl, i;

  if (! data['have_tss_info']) {
    toggle_class($("div_show_tss_id"), "hide", 1);
  }

  $("filter_top").disabled = !($("filter_on_top").checked);
  $("filter_gene_id").disabled = !($("filter_on_gene_id").checked);
  $("filter_gene_name").disabled = !($("filter_on_gene_name").checked);
  $("filter_tss_id").disabled = !($("filter_on_tss_id").checked);
  $("filter_tss_locus").disabled = !($("filter_on_tss_locus").checked);
  $("filter_re_locus").disabled = !($("filter_on_re_locus").checked);
  $("filter_absolute_distance_ge").disabled = !($("filter_on_absolute_distance").checked);
  $("filter_absolute_distance_le").disabled = !($("filter_on_absolute_distance").checked);
  $("filter_distance_ge").disabled = !($("filter_on_distance").checked);
  $("filter_distance_le").disabled = !($("filter_on_distance").checked);
  $("filter_closest_locus").disabled = !($("filter_on_closest_locus").checked);
  $("filter_closest_tss").disabled = !($("filter_on_closest_tss").checked);
  if (have_tissue_panel) {
    $("filter_histone").disabled = !($("filter_on_histone").checked);
    $("filter_correlation_sign").disabled = !($("filter_on_correlation_sign").checked);
    $("filter_corr_pvalue").disabled = !($("filter_on_corr_pvalue").checked);
    $("filter_dist_pvalue").disabled = !($("filter_on_dist_pvalue").checked);
    $("filter_cnd_pvalue").disabled = !($("filter_on_cnd_pvalue").checked);
    $("filter_qvalue").disabled = !($("filter_on_qvalue").checked);
  }

  tbl = $("regulatory_links");
  toggle_class(tbl, "hide_gene_id", !$("show_gene_id").checked);
  toggle_class(tbl, "hide_gene_name", !$("show_gene_name").checked);
  toggle_class(tbl, "hide_tss_id", !$("show_tss_id").checked || !data["have_tss_info"]);
  toggle_class(tbl, "hide_strand", !$("show_strand").checked);
  toggle_class(tbl, "hide_tss_locus", !$("show_tss_locus").checked);
  toggle_class(tbl, "hide_re_locus", !$("show_re_locus").checked);
  toggle_class(tbl, "hide_distance", !$("show_distance").checked);
  toggle_class(tbl, "hide_closest_locus", !$("show_closest_locus").checked);
  toggle_class(tbl, "hide_closest_tss", !$("show_closest_tss").checked);
  if (have_tissue_panel) {
    toggle_class(tbl, "hide_max_expr", !$("show_max_expr").checked);
    toggle_class(tbl, "hide_max_hist", !$("show_max_hist").checked);
    toggle_class(tbl, "hide_histone", !$("show_histone").checked);
    toggle_class(tbl, "hide_correlation", !$("show_correlation").checked);
    toggle_class(tbl, "hide_correlation_sign", !$("show_correlation_sign").checked);
    toggle_class(tbl, "hide_corr_pvalue", !$("show_corr_pvalue").checked);
    toggle_class(tbl, "hide_dist_pvalue", !$("show_dist_pvalue").checked);
    toggle_class(tbl, "hide_cnd_pvalue", !$("show_cnd_pvalue").checked);
    toggle_class(tbl, "hide_qvalue", !$("show_qvalue").checked);
  }

  // Hide the tissue panel-related items if no tissue panel used.
  if (! have_tissue_panel) {
    toggle_class($("regulatory_links"), "hide_tissue_panel", 1);
    toggle_class($("filters"), "hide_tissue_panel", 1);
    toggle_class($("display_columns"), "hide_tissue_panel", 1);
    toggle_class($("tbl_settings"), "hide_tissue_panel", 1);
  }

  make_regulatory_links_table(false);
} // post_load_setup

/*
 * toggle_filter
 *
 * Called when the user clicks a checkbox
 * to enable/disable a filter option.
 */
function toggle_filter(chkbox, filter_id) {
  var filter = $(filter_id);
  filter.disabled = !(chkbox.checked);
  if (!filter.disabled) {
    filter.focus();
    if (filter.select) filter.select();
  }
} // toggle_filter

/*
 * enable_filter
 *
 * Called when the user clicks a filter label.
 * Enables the filter.
 */
function enable_filter(chkbox_id, filter_id) {
  var chkbox = $(chkbox_id);
  if (!chkbox.checked) {
    var filter = $(filter_id);
    $(chkbox_id).checked = true;
    filter.disabled = false;
    filter.focus();
    if (filter.select) filter.select();
  }
} // enable_filter

/*
 * update_filter
 *
 * If the key event is an enter key press then
 * update the filter on the regulatory links table
 */
function update_filter(e) {
  if (!e) var e = window.event;
  var code = (e.keyCode ? e.keyCode : e.which);
  if (code == 13) {
    e.preventDefault();
    make_regulatory_links_table(false);
  }
} // update_filter

function num_keys(e) {
  if (!e) var e = window.event;
  var code = (e.keyCode ? e.keyCode : e.which);
  var keychar = String.fromCharCode(code);
  var numre = /\d/;
  // only allow 0-9 and various control characters (Enter, backspace, delete)
  if (code != 8 && code != 9 && code != 13 && code != 46 && !numre.test(keychar)) {
    e.preventDefault();
  }
}

/*
 * toggle_column
 *
 * Adds or removes a class from the table displaying the
 * predicted regulatory links. This is primary used to set the visibility
 * of columns by using css rules. If the parameter 'show' is not passed
 * then the existence of the class will be toggled, otherwise it will be
 * included if show is false.
 */
function toggle_column(cls) {
  toggle_class($("regulatory_links"), cls);
} // toggle_column

function populate_sort_list (sellist, items) {
  var i, j, item, opt, priority, selected;
  priority = 0;
  selected = 0;
  for (i = 0, j = 0; i < items.length; i++) {
    item = items[i];
    if (typeof item["show"] === 'undefined' || item["show"]) {
      opt = document.createElement("option");
      opt.innerHTML = item["name"];
      opt.value = i;
      sellist.add(opt, null);
      if (typeof item["priority"] !== 'undefined' && item["priority"] > priority) {
        selected = j;
        priority = item["priority"];
      }
      j++;
    }
  }
  sellist.selectedIndex = selected;
} // populate_sort_list

function populate_sort_lists() {
  "use strict";
  var i, motif_sort, regulatory_link_sort_list;
  regulatory_link_sort = $("regulatory_link_sort");
  for (i = regulatory_link_sort.options.length-1; i >= 0; i--) regulatory_link_sort.remove(i);
  regulatory_link_sort_list = sort_table["regulatory_link"];
  populate_sort_list(regulatory_link_sort, regulatory_link_sort_list);
} // populate sort lists

// Sort by increasing Gene ID, then by CnD pvalue etc.
function sort_gene_id(link1, link2) {
  var diff;
  diff = link1['gene_id'].localeCompare(link2['gene_id']);
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // sort_gene_id

// Sort by increasing Gene Name, then by CnD pvalue etc.
function sort_gene_name(link1, link2) {
  var diff;
  diff = link1['gene_name'].localeCompare(link2['gene_name']);
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // sort_gene_name

// Sort by increasing TSS ID, then CnD p-value, then RE_locus then Histone.
function sort_tss_id(link1, link2) {
  var diff;
  diff = link1['tss_id'].localeCompare(link2['tss_id']);
  if (diff != 0) return diff;
  // Cannot call sort_cnd_pvalue() or we'll loop.
  diff = have_tissue_panel ? link1['cnd_pvalue'] - link2['cnd_pvalue'] : link1['dist_pvalue'] - link2['dist_pvalue'];
  if (diff != 0) return diff;
  diff = link1['re_locus'].localeCompare(link2['re_locus']);
  if (diff != 0) return diff;
  return link1['histone'].localeCompare(link2['histone']);
} // sort_tss_id

// Sort by increasing RE_Locus, then by CnD pvalue etc.
function sort_re_locus(link1, link2) {
  var diff;
  diff = link1['re_locus'].localeCompare(link2['re_locus']);
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // 

// Sort by increasing distance, then CnD pvalue etc.
function sort_distance(link1, link2) {
  var diff;
  diff = link1['distance'] - link2['distance'];
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // sort_distance

// Sort by increasing absolute distance, then CnD pvalue etc.
function sort_absolute_distance(link1, link2) {
  var diff;
  diff = Math.abs(link1['distance']) - Math.abs(link2['distance']);
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // sort_absolute_distance

// Sort by closest TSS flag (T less than F), then CnD pvalue etc.
function sort_closest_tss(link1, link2) {
  var diff;
  diff = link2['closest_tss'].localeCompare(link1['closest_tss']);
  if (diff != 0) return diff;
  return have_tissue_panel ? sort_cnd_pvalue(link1, link2) : sort_dist_pvalue(link1, link2);
} // sort_closest_tss

// Sort by correlation decreasing, then CnD pvalue etc.
function sort_correlation_decr(link1, link2) {
  var diff, sign1, sign2;
  diff = link2['correlation'] - link1['correlation'];
  if (diff != 0) return diff;
  return sort_cnd_pvalue(link1, link2);
} // sort_correlation_decr

// Sort by correlation inreasing, then CnD pvalue etc.
function sort_correlation_incr(link1, link2) {
  var diff, sign1, sign2;
  diff = link1['correlation'] - link2['correlation'];
  if (diff != 0) return diff;
  return sort_cnd_pvalue(link1, link2);
} // sort_correlation_incr

// Sort by increasing CnD pvalue, then correlation pvalue
function sort_cnd_pvalue(link1, link2) {
  var diff;
  diff = link1['cnd_pvalue'] - link2['cnd_pvalue'];
  if (diff != 0) return diff;
  return sort_corr_pvalue(link1, link2);
} // sort_cnd_pvalue

// Sort by increasing correlation pvalue, then TSS ID (etc)
function sort_corr_pvalue(link1, link2) {
  var diff;
  diff = link1['corr_pvalue'] - link2['corr_pvalue'];
  if (diff != 0) return diff;
  var corr1 = link1['correlation']; 
  var corr2 = link2['correlation']; 
  if (corr1 < 0) corr1 = -corr1;
  if (corr2 < 0) corr2 = -corr2;
  diff = corr2 - corr1;
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_corr_pvalue

// Sort by increasing distance pvalue, then TSS ID (etc)
function sort_dist_pvalue(link1, link2) {
  var diff;
  diff = link1['dist_pvalue'] - link2['dist_pvalue'];
  if (diff != 0) return diff;
  return sort_tss_id(link1, link2);
} // sort_dist_pvalue

/*
 * regulatory_link_sort_cmp
 *
 * Gets the sorting comparator by index.
 *
 */
function regulatory_link_sort_cmp(index) {
  "use strict";
  if (index < sort_table["regulatory_link"].length) {
    return sort_table["regulatory_link"][index]["fn"];
  }
  return have_tissue_panel ? sort_cnd_pvalue : sort_dist_pvalue;	// default sorting function
}

function get_filter() {
  var filter, pat, value, value1, value2, count;

  filter = {};

  // get the maximum expression filter
  filter["on_max_expr"] = $("filter_on_max_expr").checked;
  if ((value = parseFloat($("filter_max_expr").value)) != null) {
    filter["max_expr"] = value;
    $("filter_max_expr").className = "";
  } else {
    filter["max_expr"] = false;
    $("filter_max_expr").className = "error";
  }

  // get the maximum histone level filter
  filter["on_max_hist"] = $("filter_on_max_hist").checked;
  if ((value = parseFloat($("filter_max_hist").value)) != null) {
    filter["max_hist"] = value;
    $("filter_max_hist").className = "";
  } else {
    filter["max_hist"] = false;
    $("filter_max_hist").className = "error";
  }

  // get the regulatory link count limit
  filter["on_count"] = $("filter_on_top").checked;
  count = parseFloat($("filter_top").value);
  if (isNaN(count) || count < 1) {
    filter["on_count"] = false;
    $("filter_top").className = "error";
  } else {
    filter["count"] = count;
    $("filter_top").className = "";
  }

  // get the filter_on_links filter
  filter["on_links"] = $("filter_on_links").checked;

  // get the filter_on_genes filter
  filter["on_genes"] = $("filter_on_genes").checked;

  // get the filter_on_tsses filter
  filter["on_tsses"] = $("filter_on_tsses").checked;

  // get the filter_on_loci filter
  filter["on_loci"] = $("filter_on_loci").checked;

  // get the gene id filter
  filter["on_gene_id"] = $("filter_on_gene_id").checked;
  pat = $("filter_gene_id").value;
  try {
    filter["gene_id"] = new RegExp(pat, "i");
    $("filter_gene_id").className = "";
  } catch (err) {
    $("filter_gene_id").className = "error";
    filter["on_gene_id"] = false;
  }

  // get the gene name filter
  filter["on_gene_name"] = $("filter_on_gene_name").checked;
  pat = $("filter_gene_name").value;
  try {
    filter["gene_name"] = new RegExp(pat, "i");
    $("filter_gene_name").className = "";
  } catch (err) {
    filter["on_gene_name"] = false;
    $("filter_gene_name").className = "error";
  }

  // get the TSS id filter
  filter["on_tss_id"] = $("filter_on_tss_id").checked;
  pat = $("filter_tss_id").value;
  try {
    filter["tss_id"] = new RegExp(pat, "i");
    $("filter_tss_id").className = "";
  } catch (err) {
    $("filter_tss_id").className = "error";
    filter["on_tss_id"] = false;
  }

  // get the TSS locus filter
  filter["on_tss_locus"] = $("filter_on_tss_locus").checked;
  pat = $("filter_tss_locus").value;
  try {
    filter["tss_locus"] = new RegExp(pat, "i");
    $("filter_tss_locus").className = "";
  } catch (err) {
    $("filter_tss_locus").className = "error";
    filter["on_tss_locus"] = false;
  }

  // get the RE locus filter
  filter["on_re_locus"] = $("filter_on_re_locus").checked;
  pat = $("filter_re_locus").value;
  try {
    filter["re_locus"] = new RegExp(pat, "i");
    $("filter_re_locus").className = "";
  } catch (err) {
    $("filter_re_locus").className = "error";
    filter["on_re_locus"] = false;
  }

  // get the absolute distance filter
  filter["on_absolute_distance"] = $("filter_on_absolute_distance").checked;
  value1 = parseFloat($("filter_absolute_distance_ge").value);
  value2 = parseFloat($("filter_absolute_distance_le").value);
  if (isNaN(value1)) {
    $("filter_absolute_distance_ge").className = "error";
  } else if (isNaN(value2)) {
    $("filter_absolute_distance_le").className = "error";
  } else {
    if (value1 < 0 || value1 > value2) {
      $("filter_absolute_distance_ge").className = "error";
    }
    if (value2 < 0) {
      $("filter_absolute_distance_le").className = "error";
    }
  }
  if (isNaN(value1) || isNaN(value2) || value1 < 0 || value2 < 0 || value1 > value2) {
    filter["on_absolute_distance"] = false;
    filter["on_absolute_distance_ge"] = false;
    filter["on_absolute_distance_le"] = false;
  } else {
    filter["on_absolute_distance_ge"] = true;
    $("filter_absolute_distance_ge").className = "";
    filter["on_absolute_distance_le"] = true;
    $("filter_absolute_distance_le").className = "";
    filter["absolute_distance_ge"] = value1;
    filter["absolute_distance_le"] = value2;
  }

  // get the distance filter
  filter["on_distance"] = $("filter_on_distance").checked;
  value1 = parseFloat($("filter_distance_ge").value);
  value2 = parseFloat($("filter_distance_le").value);
  if (isNaN(value1)) {
    $("filter_distance_ge").className = "error";
  } else if (isNaN(value2)) {
    $("filter_distance_le").className = "error";
  } else {
    if (value1 > value2) {
      $("filter_distance_ge").className = "error";
    }
  }
  if (isNaN(value1) || isNaN(value2) || value1 > value2) {
    filter["on_distance"] = false;
    filter["on_distance_ge"] = false;
    filter["on_distance_le"] = false;
  } else {
    filter["on_distance_ge"] = true;
    $("filter_distance_ge").className = "";
    filter["on_distance_le"] = true;
    $("filter_distance_le").className = "";
    filter["distance_ge"] = value1;
    filter["distance_le"] = value2;
  }

  // get the closest locus filter
  filter["on_closest_locus"] = $("filter_on_closest_locus").checked;
  pat = $("filter_closest_locus").value;
  try {
    filter["closest_locus"] = new RegExp(pat);
    $("filter_closest_locus").className = "";
  } catch (err) {
    $("filter_closest_locus").className = "error";
    filter["on_closest_locus"] = false;
  }

  // get the closest TSS filter
  filter["on_closest_tss"] = $("filter_on_closest_tss").checked;
  pat = $("filter_closest_tss").value;
  try {
    filter["closest_tss"] = new RegExp(pat);
    $("filter_closest_tss").className = "";
  } catch (err) {
    $("filter_closest_tss").className = "error";
    filter["on_closest_tss"] = false;
  }

  // get the histone filter
  filter["on_histone"] = $("filter_on_histone").checked;
  pat = $("filter_histone").value;
  try {
    filter["histone"] = new RegExp(pat, "i");
    $("filter_histone").className = "";
  } catch (err) {
    $("filter_histone").className = "error";
    filter["on_histone"] = false;
  }

  // get the correlation sign filter
  filter["on_correlation_sign"] = $("filter_on_correlation_sign").checked;
  value = $("filter_correlation_sign").value;
  if (value != "+" && value != "-") {
    $("filter_correlation_sign").className = "error";
    filter["on_correlation_sign"] = false;
  } else {
    $("filter_correlation_sign").className = "";
    filter["correlation_sign"] = value[0];
  }

  // get the correlation pvalue filter
  filter["on_corr_pvalue"] = $("filter_on_corr_pvalue").checked;
  if ((value = parseFloat($("filter_corr_pvalue").value)) != null) {
    filter["corr_pvalue"] = value;
    $("filter_corr_pvalue").className = "";
  } else {
    filter["corr_pvalue"] = false;
    $("filter_corr_pvalue").className = "error";
  }

  // get the distance p-value filter
  filter["on_dist_pvalue"] = $("filter_on_dist_pvalue").checked;
  if ((value = parseFloat($("filter_dist_pvalue").value)) != null) {
    filter["dist_pvalue"] = value;
    $("filter_dist_pvalue").className = "";
  } else {
    filter["dist_pvalue"] = false;
    $("filter_dist_pvalue").className = "error";
  }

  // get the CnD pvalue filter
  filter["on_cnd_pvalue"] = $("filter_on_cnd_pvalue").checked;
  if ((value = parseFloat($("filter_cnd_pvalue").value)) != null) {
    filter["cnd_pvalue"] = value;
    $("filter_cnd_pvalue").className = "";
  } else {
    filter["cnd_pvalue"] = false;
    $("filter_cnd_pvalue").className = "error";
  }

  // get the q-value filter
  filter["on_qvalue"] = $("filter_on_qvalue").checked;
  if ((value = parseFloat($("filter_qvalue").value)) != null) {
    filter["qvalue"] = value;
    $("filter_qvalue").className = "";
  } else {
    filter["qvalue"] = false;
    $("filter_qvalue").className = "error";
  }

  return filter;
} //get_filter

function filter_regulatory_link(filter, regulatory_link) {
  if (filter["on_gene_id"] && !filter["gene_id"].test(regulatory_link["gene_id"])) return true;
  if (filter["on_gene_name"] && !filter["gene_name"].test(regulatory_link["gene_name"])) return true;
  if (filter["on_tss_id"] && !filter["tss_id"].test(regulatory_link["tss_id"])) return true;
  if (filter["on_tss_locus"] && !filter["tss_locus"].test(regulatory_link["tss_locus"])) return true;
  if (filter["on_re_locus"] && !filter["re_locus"].test(regulatory_link["re_locus"])) return true;
  if (filter["on_absolute_distance"] && (
    (Math.abs(regulatory_link["distance"]) < filter["absolute_distance_ge"]) ||
    (Math.abs(regulatory_link["distance"]) > filter["absolute_distance_le"]))
  ) return true;
  if (filter["on_distance"] && (
    (regulatory_link["distance"] < filter["distance_ge"]) ||
    (regulatory_link["distance"] > filter["distance_le"]))
  ) return true;
  if (filter["on_closest_locus"] && !filter["closest_locus"].test(regulatory_link["closest_locus"])) return true;
  if (filter["on_closest_tss"] && !filter["closest_tss"].test(regulatory_link["closest_tss"])) return true;
  if (filter["on_histone"] && !filter["histone"].test(regulatory_link["histone"])) return true;
  if (filter["on_correlation_sign"] && 
    ((regulatory_link["correlation"] > 0 && filter["correlation_sign"] != "+") ||
     (regulatory_link["correlation"] <= 0 && filter["correlation_sign"] != "-") )
  ) return true;
  if (filter["on_max_expr"] && regulatory_link["max_expr"] < filter["max_expr"]) return true;
  if (filter["on_max_hist"] && regulatory_link["max_hist"] < filter["max_hist"]) return true;
  if (filter["on_corr_pvalue"] && regulatory_link["corr_pvalue"] > filter["corr_pvalue"]) return true;
  if (filter["on_dist_pvalue"] && regulatory_link["dist_pvalue"] > filter["dist_pvalue"]) return true;
  if (filter["on_cnd_pvalue"] && regulatory_link["cnd_pvalue"] > filter["cnd_pvalue"]) return true;
  if (filter["on_qvalue"] && regulatory_link["qvalue"] > filter["qvalue"]) return true;
  return false;
}

function filter_regulatory_links_table() {
  "use strict";
  var tbl, tbody, tr, regulatory_link, i, filter;
  tbl = $("regulatory_links");
  filter = get_filter();
  for (i = 0; i < tbl.tBodies.length; i++) {
    tbody = tbl.tBodies[i];
    regulatory_link = tbody["data_regulatory_link"];
    toggle_class(tbody, "filtered", filter_motif(filter, regulatory_link));
  }
}

function make_other_settings() {
  $("opt_transcript_types").textContent = data.options.transcript_types;
  $("opt_max_link_distances").textContent = data.options.max_link_distances;
  $("opt_max_pvalue").textContent = data.options.max_pvalue;
  $("opt_histone_root").textContent = data.options.histone_root;
  $("opt_histones").textContent = data.options.histones;
  $("opt_tissues").textContent = data.options.tissues;
  $("opt_rna_source").textContent = data.options.rna_source;
  $("opt_expression_root").textContent = data.options.expression_root;
  $("opt_use_gene_ids").textContent = data.options.use_gene_ids;
  $("opt_lecat").textContent = data.options.lecat;
  $("opt_inc_closest_locus").textContent = data.options.inc_closest_locus;
  $("opt_inc_closest_tss").textContent = data.options.inc_closest_tss;
  $("opt_noise").textContent = data.options.noise;
  $("opt_seed").textContent = data.options.seed;
  $("opt_n_perms").textContent = data.n_perms;
  $("opt_noise_fraction").textContent = data.noise_fraction;
} // make_other_settings

function make_regulatory_links_row(tbody, link) {
  var max_expr = have_tissue_panel ? link['max_expr'].toFixed(2) : '';
  var max_hist = have_tissue_panel ? link['max_hist'].toFixed(2) : '';
  var histone = have_tissue_panel ? link['histone'] : '';
  var correlation = have_tissue_panel ? link['correlation'].toFixed(5) : '';
  var corr_sign = have_tissue_panel ? (link['correlation'] > 0 ? '+' : '-') : '';
  var corr_pvalue = have_tissue_panel ? link['corr_pvalue'].toExponential(2) : '';
  var cnd_pvalue = have_tissue_panel ? link['cnd_pvalue'].toExponential(2) : '';
  var row;
  row = tbody.insertRow(tbody.rows.length);
  add_text_cell(row, link['gene_id'], 'col_gene_id');
  add_text_cell(row, link['gene_name'], 'col_gene_name');
  add_text_cell(row, link['tss_id'], 'col_tss_id');
  add_text_cell(row, link['tss_locus'], 'col_tss_locus');
  add_text_cell(row, link['strand'], 'col_strand');
  add_text_cell(row, max_expr, 'col_max_expr tissue_panel');
  add_text_cell(row, link['re_locus'], 'col_re_locus');
  add_text_cell(row, max_hist, 'col_max_hist tissue_panel');
  add_text_cell(row, link['distance'], 'col_distance');
  add_text_cell(row, link['closest_locus'], 'col_closest_locus');
  add_text_cell(row, link['closest_tss'], 'col_closest_tss');
  add_text_cell(row, histone, 'col_histone tissue_panel');
  add_text_cell(row, correlation, 'col_correlation tissue_panel');
  add_text_cell(row, corr_sign, 'col_correlation_sign tissue_panel');
  add_text_cell(row, corr_pvalue, 'col_corr_pvalue tissue_panel');
  add_text_cell(row, link['dist_pvalue'].toExponential(2), 'col_dist_pvalue');
  add_text_cell(row, cnd_pvalue, 'col_cnd_pvalue tissue_panel');
  add_text_cell(row, link['qvalue'].toExponential(2), 'col_qvalue');
} // make_regulatory_links_row

function get_regulatory_link_line(link) {
  var line = "";
  var tab = "";
  var max_expr = have_tissue_panel ? link['max_expr'].toFixed(2) : '';
  var max_hist = have_tissue_panel ? link['max_hist'].toFixed(2) : '';
  var histone = have_tissue_panel ? link['histone'] : '';
  var correlation = have_tissue_panel ? link['correlation'].toFixed(5) : '';
  var corr_sign = have_tissue_panel ? (link['correlation'] > 0 ? '+' : '-') : '';
  var corr_pvalue = have_tissue_panel ? link['corr_pvalue'].toExponential(2) : '';
  var cnd_pvalue = have_tissue_panel ? link['cnd_pvalue'].toExponential(2) : '';
  if ($('show_gene_id').checked) { line = line + tab + link['gene_id']; tab = '\t'; }
  if ($('show_gene_name').checked) { line = line + tab + link['gene_name']; tab = '\t'; }
  if ($('show_tss_id').checked) { line = line + tab + link['tss_id']; tab = '\t'; }
  if ($('show_tss_locus').checked) { line = line + tab + link['tss_locus']; tab = '\t'; }
  if ($('show_strand').checked) { line = line + tab + link['strand']; tab = '\t'; }
  if ($('show_max_expr').checked) { line = line + tab + max_expr; tab = '\t'; }
  if ($('show_re_locus').checked) { line = line + tab + link['re_locus']; tab = '\t'; }
  if ($('show_max_hist').checked) { line = line + tab + max_hist; tab = '\t'; }
  if ($('show_distance').checked) { line = line + tab + link['distance']; tab = '\t'; }
  if ($('show_closest_locus').checked) { line = line + tab + link['closest_locus']; tab = '\t'; }
  if ($('show_closest_tss').checked) { line = line + tab + link['closest_tss']; tab = '\t'; }
  if ($('show_histone').checked) { line = line + tab + histone; tab = '\t'; }
  if ($('show_correlation').checked) { line = line + tab + correlation; tab = '\t'; }
  if ($('show_correlation_sign').checked) { line = line + tab + corr_sign; tab = '\t'; }
  if ($('show_corr_pvalue').checked) { line = line + tab + corr_pvalue; tab = '\t'; }
  if ($('show_dist_pvalue').checked) { line = line + tab + link['dist_pvalue'].toExponential(2); tab = '\t'; }
  if ($('show_cnd_pvalue').checked) { line = line + tab + cnd_pvalue; tab = '\t'; }
  if ($('show_qvalue').checked) { line = line + tab + link['qvalue'].toExponential(2); tab = '\t'; }
  line = line + '\n';
  return(line);
} // get_regulatory_link_line

function get_regulatory_link_hdr() {
  var line = '';
  var tab = '';
  var field_names = {
    "gene_id" : "Gene_ID",
    "gene_name" : "Gene_Name",
    "tss_id" : "TSS_ID",
    "tss_locus" : "TSS_Locus",
    "strand" : "Strand",
    "max_expr" : "Max_Expr",
    "re_locus" : "RE_Locus",
    "max_hist" : "Max_Hist",
    "distance" : "Distance",
    "closest_locus" : "Closest_Locus",
    "closest_tss" : "Closest_TSS",
    "histone" : "Histone",
    "correlation" : "Correlation",
    "correlation_sign" : "Correlation_Sign",
    "corr_pvalue" : "Correlation_P_Value",
    "dist_pvalue" : "Distance_P_Value",
    "cnd_pvalue" : "CnD_P_Value",
    "qvalue" : "Q_Value"
  };
  
  if ($('show_gene_id').checked) { line = line + tab + field_names['gene_id']; tab = '\t'; }
  if ($('show_gene_name').checked) { line = line + tab + field_names['gene_name']; tab = '\t'; }
  if ($('show_tss_id').checked) { line = line + tab + field_names['tss_id']; tab = '\t'; }
  if ($('show_tss_locus').checked) { line = line + tab + field_names['tss_locus']; tab = '\t'; }
  if ($('show_strand').checked) { line = line + tab + field_names['strand']; tab = '\t'; }
  if ($('show_max_expr').checked) { line = line + tab + field_names['max_expr']; tab = '\t'; }
  if ($('show_re_locus').checked) { line = line + tab + field_names['re_locus']; tab = '\t'; }
  if ($('show_max_hist').checked) { line = line + tab + field_names['max_hist']; tab = '\t'; }
  if ($('show_distance').checked) { line = line + tab + field_names['distance']; tab = '\t'; }
  if ($('show_closest_locus').checked) { line = line + tab + field_names['closest_locus']; tab = '\t'; }
  if ($('show_closest_tss').checked) { line = line + tab + field_names['closest_tss']; tab = '\t'; }
  if ($('show_histone').checked) { line = line + tab + field_names['histone']; tab = '\t'; }
  if ($('show_correlation').checked) { line = line + tab + field_names['correlation']; tab = '\t'; }
  if ($('show_correlation_sign').checked) { line = line + tab + field_names['correlation_sign']; tab = '\t'; }
  if ($('show_corr_pvalue').checked) { line = line + tab + field_names['corr_pvalue']; tab = '\t'; }
  if ($('show_dist_pvalue').checked) { line = line + tab + field_names['dist_pvalue']; tab = '\t'; }
  if ($('show_cnd_pvalue').checked) { line = line + tab + field_names['cnd_pvalue']; tab = '\t'; }
  if ($('show_qvalue').checked) { line = line + tab + field_names['qvalue']; tab = '\t'; }
  line = line + "\n";
  return(line);
} // get_regulatory_link_hdr

function make_regulatory_links_table(download) {
  var regulatory_link_sort;
  var filter, regulatory_links, filtered, regulatory_link;
  var tbl, tbody, row, cell, i, skipped;

  // get the filter and sort comparator
  filter = get_filter();
  regulatory_link_sort = regulatory_link_sort_cmp(parseInt($('regulatory_link_sort').value));
 
  // Apply the pre-sort filters.
  regulatory_links = data['regulatory_links']; 
  pre_sort_filtered = [];
  for (i = 0; i < regulatory_links.length; i++) {
    if (filter_regulatory_link(filter, regulatory_links[i])) continue;
    if (regulatory_links[i]['filtered']) continue;
    pre_sort_filtered.push(regulatory_links[i]);
  }

  // Sort the filtered links.
  pre_sort_filtered.sort(regulatory_link_sort);

  // Apply the post-sort filters.
  filtered = [];
  best_tss_for_gene = {};
  best_re_for_tss = {};
  best_tss_for_re = {};
  for (i = 0; i < pre_sort_filtered.length; i++) {
    link = pre_sort_filtered[i];
    gene_id = link['gene_id'];
    tss_id = link['tss_id'];
    re_locus = link['re_locus'];
    // Show only best TSS for each Gene (TSS may have multiple RE loci).
    if (filter['on_genes']) {
      if (best_tss_for_gene[gene_id] && best_tss_for_gene[gene_id] != tss_id) continue;
    }
    // Show only best RE locus for each TSS.
    if (filter['on_tsses'] || filter['on_genes']) {
      if (best_re_for_tss[tss_id]) continue;
    }
    // Show only best TSS for each RE locus.
    if (filter['on_loci']) {
      if (best_tss_for_re[re_locus]) continue;
    }
    best_tss_for_gene[gene_id] = tss_id;
    best_re_for_tss[tss_id] = true;
    best_tss_for_re[re_locus] = true;
    filtered.push(link);
  }

  // limit to top N regulatory links
  if (filter['on_count']) {
    if (filtered.length > filter['count']) filtered.length = filter['count'];
  }

  // clear the table
  tbl = $('regulatory_links');
  for (i = tbl.tBodies.length - 1; i >= 0; i--) {
    tbody = tbl.tBodies[i];
    tbody.parentNode.removeChild(tbody);
  }

  // add the new rows to the table
  for (i = 0; i < filtered.length; i++) {
    regulatory_link = filtered[i];
    tbody = document.createElement('tbody');
    tbody.className = "regulatory_link_group";
    tbody['data_regulatory_link'] = regulatory_link;
    make_regulatory_links_row(tbody, regulatory_link);
    tbl.appendChild(tbody);
  }

  // download (displayed columns of) the (filtered and sorted) links
  if (download) {
    var links_text = get_regulatory_link_hdr();
    for (i = 0; i < filtered.length; i++) {
      links_text = links_text + get_regulatory_link_line(filtered[i]);
    }
    var links_file_name = "tgene_links.tsv";
    var mimetype = "text/plain";
    //console.log(links_text);
    prepare_download(links_text, mimetype, links_file_name);
  }

  // note the count of filtered regulatory links
  if (filtered.length != regulatory_links.length) {
    skipped =  regulatory_links.length - filtered.length;
    tbody = document.createElement('tbody');
    row = tbody.insertRow(tbody.rows.length);

    if (skipped === 1) {
      desc = "1 potential regulatory link hidden due to filters";
    } else {
      desc = skipped + " potential regulatory links hidden due to filters";
    }
    cell = row.insertCell(row.cells.length);
    cell.colSpan = 19;
    cell.style.textAlign = "center";
    cell.style.fontWeight = "bold";
    cell.appendChild(document.createTextNode(desc));
    tbl.appendChild(tbody);
  }
} // make_regulatory_links_table
