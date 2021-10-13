var current_program = "MEME";
var current_alphabet = new Alphabet(data.alphabet, data.background.freqs);
var current_motif = 0;
//var new_icon_src = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABwAAAAQCAMAAAAyEe/dAAAACVBMVEX//wAAAAD////s2cV/AAAAdUlEQVQYlXVRBxLAIAhL+P+jC2HZhXcBZEWEldDsZcLIcAhHWWnK8SDcWQhMFUHdAQ1CqQ5+CWPmlHojl+nCJNRtzu4qRc3IUzmTVpXYK0nox0z0PI1stgchdK7lEv7ekhvalw8WW547Gyzedt/2/gLx8WXjXF/1AYFriNAWAAAAAElFTkSuQmCC";

var DelayLogoTask = function(logo, canvas) {
  this.logo = logo;
  this.canvas = canvas;
};

DelayLogoTask.prototype.run = function () {
  draw_logo_on_canvas(this.logo, this.canvas, false);
};

function clone_template(template) {
  "use strict";
  var node, help_btns, i, button;
  node = $(template).cloneNode(true);
  toggle_class(node, "template", false);
  node.id = "";
  help_btns = node.querySelectorAll(".help");
  for (i = 0; i < help_btns.length; i++) {
    button = help_btns[i];
    if (button.hasAttribute("data-topic")) {
      button.tabIndex = "0";
      button.addEventListener("click", __toggle_help, false);
      button.addEventListener("keydown", __toggle_help, false);
    }
  }
  return node;
}

function make_small_logo(alphabet, pspm, options) {
  if (typeof options === "undefined") options = {};
  if (options.rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, {x_axis: false, y_axis: false});
  logo.add_pspm(pspm, (typeof options.offset === "number" ? options.offset : 0));
  var canvas = document.createElement('canvas');
  if (typeof options.className === "string") canvas.className = options.className;
  if (typeof options.width === "number" && options.width > 0) {
    canvas.height = 0;
    canvas.width = options.width;
    draw_logo_on_canvas(logo, canvas, false);
  } else {
    draw_logo_on_canvas(logo, canvas, false, 1/3);
  }
  return canvas;
}

function make_large_logo(alphabet, pspm, rc, offset, className) {
  if (rc) pspm = pspm.copy().reverse_complement(alphabet);
  var logo = new Logo(alphabet, "");
  logo.add_pspm(pspm, offset);
  var canvas = document.createElement('canvas');
  canvas.height = 200;
  canvas.width = 0;
  canvas.className = className;
  size_logo_on_canvas(logo, canvas, false);
  add_draw_task(canvas, new DelayLogoTask(logo, canvas));
  return canvas;
}

function make_sym_btn(symbol, title, action) {
  var box;
  box = document.createElement("div");
  box.tabIndex = 0;
  box.className = "sym_btn";
  box.appendChild(document.createTextNode(symbol));
  box.title = title;
  box.addEventListener('click', action, false);
  box.addEventListener('keydown', action, false);
  return box;
}

function make_seq(alphabet, seq) {
  var i, j, letter, lbox, sbox;
  sbox = document.createElement("span");
  for (i = 0; i < seq.length; i = j) {
    letter = seq.charAt(i);
    for (j = i+1; j < seq.length; j++) {
      if (seq.charAt(j) !== letter) {
        break;
      }
    }
    lbox = document.createElement("span");
    lbox.style.color = alphabet.get_colour(alphabet.get_index(letter));
    lbox.appendChild(document.createTextNode(seq.substring(i, j)));
    sbox.appendChild(lbox);
  }
  return sbox;
}

//
// make_pv_text
//
// Returns the string p-value, with the p italicised.
///
function make_pv_text() {
  var pv_text = document.createElement("span");
  var pv_italic_text = document.createElement("span");
  pv_italic_text.appendChild(document.createTextNode("p"));
  pv_italic_text.style.fontStyle = "italic";
  pv_text.appendChild(pv_italic_text);
  pv_text.appendChild(document.createTextNode("-value"));
  return pv_text;
}

function append_site_entries(tbody, motif, site_index, count) {
  "use strict";
  var i, end;
  var sites, site, sequences, sequence;
  var rbody;
  if (typeof count !== "number") {
    count = 20;
  }
  sequences = data["sequence_db"]["sequences"];
  sites = motif["sites"];
  end = Math.min(site_index + count, sites.length);
  for (i = site_index; i < end; i++) {
    site = sites[i];
    sequence = sequences[site["seq"]];

    rbody = tbody.insertRow(tbody.rows.length);
    add_text_cell(rbody, "" + (site["seq"] + 1) + ".", "site_num");
    add_text_cell(rbody, sequence["name"], "site_name");
    add_text_cell(rbody, site["rc"] ? "-" : "+", "site_strand");
    add_text_cell(rbody, site["pos"] + 1, "site_start");
    add_text_cell(rbody, site["pvalue"].toExponential(2), "site_pvalue");
    add_text_cell(rbody, site["lflank"], "site lflank");
    add_cell(rbody, make_seq(current_alphabet, site["match"]), "site match");
    add_text_cell(rbody, site["rflank"], "site rflank");
  }
  return i;
}

function make_site_entries() {
  "use strict";
  var region;
  region = this;
  if (region.data_site_index >= region.data_motif["sites"].length) {
    // all sites created
    region.removeEventListener('scroll', make_site_entries, false);
    return;
  }
  // if there's still 100 pixels to scroll than don't do anything yet
  if (region.scrollHeight - (region.scrollTop + region.offsetHeight) > 100) {
    return;
  }

  region.data_site_index = append_site_entries(
      find_child(region, "sites_tbl").tBodies[0], 
      region.data_motif, region.data_site_index, 20
    ); 
}

function make_sites(motif) {
  "use strict";
  function add_site_header(row, title, nopad, help_topic, tag_class) {
    var div, divcp, th;
    th = document.createElement("th");
    div = document.createElement("div");
    div.className = "sites_th_inner";
    if (typeof title !== "object") {
      title = document.createTextNode("" + title);
    }
    div.appendChild(title);
    if (help_topic) {
      div.appendChild(document.createTextNode("\xA0"));
      div.appendChild(help_button(help_topic));
    }
    divcp = div.cloneNode(true);
    divcp.className = "sites_th_hidden";
    th.appendChild(div);
    th.appendChild(divcp);
    if (nopad) {
      th.className = "nopad";
    }
    if (tag_class) {
      th.className += " " + tag_class;
    }
    row.appendChild(th);
  }
  var outer_tbl, inner_tbl, tbl, thead, tbody, rhead;

  outer_tbl = document.createElement("div");
  outer_tbl.className = "sites_outer";

  inner_tbl = document.createElement("div");
  inner_tbl.className = "sites_inner";
  outer_tbl.appendChild(inner_tbl);

  tbl = document.createElement("table");
  tbl.className = "sites_tbl";
  inner_tbl.appendChild(tbl);

  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);

  rhead = thead.insertRow(thead.rows.length);
  add_site_header(rhead, "", true);
  add_site_header(rhead, "Name", false, "pop_seq_name");
  add_site_header(rhead, "Strand", false, "pop_site_strand", "site_strand_title");
  add_site_header(rhead, "Start", false, "pop_site_start");
  add_site_header(rhead, make_pv_text(), false, "pop_site_pvalue");
  add_site_header(rhead, "", false);
  add_site_header(rhead, "Sites", true, "pop_site_match");
  add_site_header(rhead, "", false);

  inner_tbl.data_motif = motif;
  inner_tbl.data_site_index = append_site_entries(tbody, motif, 0, 20);
  if (inner_tbl.data_site_index < motif["sites"].length) {
    inner_tbl.addEventListener('scroll', make_site_entries, false);
  }
  return outer_tbl;
}

function make_motif_table_entry(row, alphabet, ordinal, motif, colw) {
  "use strict";
  function ev_sig(evalue_str) {
    "use strict";
    var ev_re, match, sig, exp, num;
    ev_re = /^(.*)e(.*)$/;
    if (match = ev_re.exec(evalue_str)) {
      sig = parseFloat(match[1]);
      exp = parseInt(match[2]);
      if (exp >= 0) {
        return false;
      } else if (exp <= -3) {
        return true;
      } else {
        return sig * Math.pow(10, exp) <= 0.05;
      }
    }
    return true;
  }
  function make_preview(alphabet, motif) {
    "use strict";
    var pspm, preview, preview_rc;
    var box, btn_box, logo_box, btn_plus, btn_minus;
    if (motif["preview_logo"]) {
      preview = motif["preview_logo"];
      preview_rc = motif["preview_logo_rc"];
    } else {
      pspm = new Pspm(motif["pwm"]);
      preview = make_logo(alphabet, pspm, 50, false, 0);
      motif["preview_logo"] = preview;
      if (alphabet.has_complement()) {
        preview_rc = make_logo(alphabet, pspm, 50, true, 0, "logo_rc");
        motif["preview_logo_rc"] = preview_rc;
      }
    }
    if (preview_rc) {
      btn_plus = document.createElement("div");
      btn_plus.appendChild(document.createTextNode("+"));
      btn_plus.className = "preview_btn plus";
      btn_plus.tabIndex = "0";
      btn_plus.addEventListener("click", action_btn_rc, false);
      btn_plus.addEventListener("keydown", action_btn_rc, false);
      btn_minus = document.createElement("div");
      btn_minus.appendChild(document.createTextNode("-"));
      btn_minus.className = "preview_btn minus";
      btn_minus.tabIndex = "0";
      btn_minus.addEventListener("click", action_btn_rc, false);
      btn_minus.addEventListener("keydown", action_btn_rc, false);
      btn_box = document.createElement("div");
      btn_box.className = "preview_btn_box";
      btn_box.appendChild(btn_plus);
      btn_box.appendChild(btn_minus);
    }
    logo_box = document.createElement("div");
    logo_box.className = "preview_logo_box";
    logo_box.appendChild(preview);
    if (preview_rc) logo_box.appendChild(preview_rc);
    box = document.createElement("div");
    box.className = "preview_box";
    if (preview_rc) box.appendChild(btn_box);
    box.appendChild(logo_box);
    if (preview_rc) {
      if (motif["rc"]) {
        btn_minus.className += " active";
        logo_box.className += " show_rc_logo";
      } else {
        btn_plus.className += " active";
      }
    }
    return box;
  }
  var pspm, preview, preview_rc, c;
  row.data_motif = motif;
  row.id = motif["alt"];
  row.data_ordinal = ordinal;
  if (!ev_sig(motif["evalue"])) {
    row.style.opacity = 0.4;
  }
  add_text_cell(row, "" + ordinal + ".", "motif_ordinal");
  add_cell(row, make_preview(alphabet, motif), "motif_logo");
  add_text_cell(row, motif["evalue"], "motif_evalue");
  add_text_cell(row, motif["nsites"], "motif_nsites");
  add_text_cell(row, motif["len"], "motif_width");
  add_cell(row, make_sym_btn("\u21A7", "Show more information.", 
    action_show_more), "motif_more");
  add_cell(row, make_sym_btn("\u21E2", "Submit the motif to another MEME Suite program or download it.", function(e) { action_show_outpop(e, ordinal); }), "motif_submit");
  if (colw) {
    for (c = 0; c < row.cells.length; c++) {
      row.cells[c].style.minWidth = colw[c] + "px";
    }
  }
}

function make_motifs_table(alphabet, start_ordinal, motifs, colw, stop_reason) {
  var i, j;
  var tbl, thead, tbody, tfoot, row, preview;
  var motif, pspm;

  tbl = document.createElement("table");
  
  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);
  tfoot = document.createElement("tfoot");
  tbl.appendChild(tfoot);

  row = thead.insertRow(thead.rows.length);
  add_text_header_cell(row, "", "", "motif_ordinal");
  add_text_header_cell(row, "Logo", "pop_logo", "motif_logo");
  add_text_header_cell(row, "E-value", "pop_ev", "motif_evalue");
  add_text_header_cell(row, "Sites", "pop_sites", "motif_nsites");
  add_text_header_cell(row, "Width", "pop_width", "motif_width");
  add_text_header_cell(row, "More", "pop_more", "motif_more");
  add_text_header_cell(row, "Submit/Download", "pop_submit_dl", "motif_submit");

  for (i = 0; i < motifs.length; i++) {
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, alphabet, start_ordinal + i, motifs[i], colw);
  }

  row = tfoot.insertRow(tfoot.rows.length);
  add_text_header_cell(row, stop_reason, "", "stop_reason", "", 6);

  return tbl;
}

function make_expanded_motif(alphabet, ordinal, motif, less_x, submit_x) {
  "use strict";
  var box, pspm, logo_box, large_logo, large_logo_rc, tab_logo, tab_logo_rc;
  var btn, offset, norc;

  box = clone_template("tmpl_motif_expanded");
  box.data_motif = motif;
  box.data_ordinal = ordinal;

  pspm = new Pspm(motif["pwm"]);
  if (typeof motif["rc"] !== "boolean") {
    motif["rc"] = false;
  }
  if (motif["large_logo"]) {
    large_logo = motif["large_logo"];
    large_logo_rc = motif["large_logo_rc"];
  } else {
    large_logo = make_large_logo(alphabet, pspm, false, 0);
    motif["large_logo"] = large_logo;
    if (alphabet.has_complement()) {
      large_logo_rc = make_large_logo(alphabet, pspm, true, 0, "logo_rc");
      motif["large_logo_rc"] = large_logo_rc;
    }
  }
  norc = (large_logo_rc == null);
  toggle_class(box, "norc", norc);

  logo_box = find_child(box, "tvar_logo");
  logo_box.appendChild(large_logo);
  if (large_logo_rc) logo_box.appendChild(large_logo_rc);
  toggle_class(logo_box, "show_rc_logo", motif["rc"]);

  tab_logo = find_child(box, "tvar_tab");
  tab_logo_rc = find_child(box, "tvar_tab_rc");

  toggle_class(tab_logo, "activeTab", !motif["rc"]);
  toggle_class(tab_logo_rc, "activeTab", motif["rc"]);

  tab_logo.addEventListener('click', action_rc_tab, false);
  tab_logo.addEventListener('keydown', action_rc_tab, false);
  tab_logo_rc.addEventListener('click', action_rc_tab, false);
  tab_logo_rc.addEventListener('keydown', action_rc_tab, false);

  set_tvar(box, "tvar_ordinal", ordinal); 
  set_tvar(box, "tvar_evalue", motif["evalue"]);
  set_tvar(box, "tvar_width", motif["len"]);
  set_tvar(box, "tvar_site_count", motif["nsites"]);
  set_tvar(box, "tvar_llr", motif["llr"]);
  set_tvar(box, "tvar_ic", motif["ic"]);
  set_tvar(box, "tvar_re", motif["re"]);
  set_tvar(box, "tvar_bt", motif["bt"]);
  if (data.sequence_db.primary_count > data.options.brief) {
    if (data.options.brief == 1000) {
      set_tvar(box, "tvar_sites", "Output of sites suppressed because there were more than 1000 (primary) sequences.");
    } else {
      set_tvar(box, "tvar_sites", "Output of sites suppressed by -brief option.");
    }
  } else {
    set_tvar(box, "tvar_sites", make_sites(motif));
  }

  offset = 32; // 1* 5px padding + 2 * 10px padding + 2 * 2px border + 3px ??

  btn = find_child(box, "tvar_less");
  btn.style.left = (less_x - offset) + "px";
  btn.addEventListener('click', action_show_less, false);
  btn.addEventListener('keydown', action_show_less, false);
  btn = find_child(box, "tvar_submit");
  btn.style.left = (submit_x - offset) + "px";
  btn.addEventListener('click', action_show_outpop, false);
  btn.addEventListener('keydown', action_show_outpop, false);
  return box;
}

//
//
//
function make_motifs() {
  "use strict";
  function pixel_value(str_in) {
    "use strict";
    var px_re, match;
    px_re = /^(\d+)px$/;
    if (match = px_re.exec(str_in)) {
      return parseInt(match[1], 10);
    }
    return 0;
  }
  var container, tbl;
  var colw, r, row, c, cell, cell_style, pad_left, pad_right;

  // make the motifs table
  container = $("motifs");
  container.innerHTML = ""; // clear content

  tbl = make_motifs_table(current_alphabet, 1, data["motifs"], colw, data["stop_reason"]);
  container.appendChild(tbl);

  // measure table column widths
  colw = [];
  row = tbl.tBodies[0].rows[0];
  for (c = 0; c < row.cells.length; c++) {
    var padLeft, padRight;
    cell = row.cells[c];
    cell_style = window.getComputedStyle(cell, null);
    pad_left = pixel_value(cell_style.getPropertyValue("padding-left"));
    pad_right = pixel_value(cell_style.getPropertyValue("padding-right"));
    colw[c] = cell.clientWidth - pad_left - pad_right;
    if (typeof colw[c] !== "number" || colw[c] < 0) {
      colw[c] = 1;
    }
  }

  // set minimum table column widths on each row so later when we remove rows it still aligns
  for (r = 0; r < tbl.tBodies[0].rows.length; r++) {
    row = tbl.tBodies[0].rows[r];
    for (c = 0; c < row.cells.length; c++) {
      row.cells[c].style.minWidth = colw[c] + "px";
    }
  }

  // store the table column widths so we can create rows latter with the same minimums
  container.data_colw = colw;

  // calculate the x offset for the buttons
  row = tbl.tBodies[0].rows[0];
  container.data_more_x = coords(find_child(find_child(row, "motif_more"), "sym_btn"))[0];
  container.data_submit_x = coords(find_child(find_child(row, "motif_submit"), "sym_btn"))[0];

  draw_on_screen();
}

function make_meme_block(container, max_seq_len, is_scan, site) {
  "use strict";
  var motif = data.motifs[site.motif];
  var block = make_block(container, max_seq_len, site.pos, motif.len,
      site.pvalue, site.rc, site.motif, is_scan);
  var handler = (is_scan ?
      make_scan_popup(site, motif, block) :
      make_block_popup(site, motif, block));
  block.addEventListener("mouseover", handler, false);
  block.addEventListener("mouseout", handler, false);
}

function append_blocks_entries(tbody, seq_index, count) {
  "use strict";
  var i, end, j;
  var max_pvalue, max_block_height, max_seq_len, sequences;
  var sequence, sites, scans, scan;
  var container, plus, minus, rule, row;
  // define some constants
  max_seq_len = data.sequence_db.max_length;
  // determine how many to load
  end = Math.min(seq_index + count, data.sequence_db.sequences.length);
  for (i = seq_index; i < end; i++) {
    // get the sequence
    sequence = data.sequence_db.sequences[i];
    // make the containers for the block diagram
    container = make_block_container(current_alphabet.has_complement(),
        data.options.revcomp, max_seq_len, sequence.length);
    // create blocks for the motif sites
    sites = sequence["sites"];
    for (j = 0; j < sites.length; j++)
      make_meme_block(container, max_seq_len, false, sites[j]);
    // create blocks for the scanned sites
    scan = data.scan[i];
    for (j = 0; j < scan.sites.length; j++)
      make_meme_block(container, max_seq_len, true, scan.sites[j]);
    // create a row for the sequence
    row = tbody.insertRow(tbody.rows.length);
    toggle_class(row, "empty_seq", sites.length == 0 && scan.sites.length == 0);
    toggle_class(row, "only_scan", sites.length == 0 && scan.sites.length > 0);
    add_text_cell(row, (i + 1) + ".", "blockdiag_num");
    add_text_cell(row, sequence["name"], "blockdiag_name");
    add_text_cell(row, scan["pvalue"].toExponential(2), "blockdiag_pvalue");
    add_cell(row, container, "block_td"); 
  }
  return end;
}

function make_blocks_entries() {
  "use strict";
  var region;
  region = this;
  if (region.data_blocks_index >= data["sequence_db"]["sequences"].length) {
    // all sites created
    region.removeEventListener('scroll', make_blocks_entries, false);
    return;
  }
  // if there's still 100 pixels to scroll then don't do anything yet
  if (region.scrollHeight - (region.scrollTop + region.offsetHeight) > 100) {
    return;
  }

  region.data_blocks_index = append_blocks_entries(
      find_child(region, "blocks_tbl").tBodies[0], 
      region.data_blocks_index, 20
    ); 
}

// Apply opacity alpha to color rgb with backrgound bkg.
function RGBAtoRGB(rgb, bkg, opacity) {
  var i;
  var rgb_new = [];
  for (i=0; i<3; i++) {
    rgb_new[i] = Math.round(((1-opacity) * bkg[i]) + (opacity * rgb[i]));
  }
  return rgb_new;
}

// Function to measure the size of text on a canvas.
var MeasureText = function(font, text) {
  var image = document.createElement("canvas");
  var image_ctx = image.getContext('2d');
  image_ctx.save();
  image_ctx.font = font;
  var text_length = image_ctx.measureText(text).width;
  image.remove();
  return text_length;
} // MeasureText

// Functions to download the motif block diagram as a PDF or SVG file.
function download_PDF_block_diagram() {
  downloadBlockDiagram(true);
}
function download_SVG_block_diagram() {
  downloadBlockDiagram(false);
}

// Helper function to download the motif block diagram as a PDF or SVG file.
var downloadBlockDiagram = function(
  make_pdf
) {

  // Check that the required javascript was loaded.
  if ( (make_pdf && typeof jsPDF === 'undefined') ||
    (!make_pdf && typeof d3 === 'undefined') ) {
    var id = make_pdf ? $("pdfButton") : $("svgButton");
    help_popup(id, "pop_offline");
    return;
  }

  // Determine which lines are visible in the HTML inner scroll window.
  var inner_tbl = $("blocks_scroll");
  var pix_per_sequence = 27;		// (vertical) pixels per sequence diagram line
  var first = Math.round(inner_tbl.scrollTop / pix_per_sequence) + 1;
  var last = first + Math.round(inner_tbl.offsetHeight / pix_per_sequence) - 1;

  // Get the contents of the HTML inner scroll window while saving the sequences to be printed.
  var numbers = document.getElementsByClassName("blockdiag_num");
  var bars = {};
  var visible_motifs = {};
  var seq_index = 0;
  var rgb;
  for (var i=0; i<numbers.length && seq_index < last; i++) {
    var row_node = numbers[i].parentNode;
    // Check if the sequence is displayed in the outer scrolling window.
    var seq_name = numbers[i].nextSibling.innerHTML;
    if (
      ($("rdo_sites_only").checked && row_node.getAttribute("class").includes("only_scan"))
      || 
      (! $("rdo_all_seqs").checked && row_node.getAttribute("class").includes("empty_seq"))
    ) { continue; }
    seq_index++;

    if (seq_index < first) { continue; }		// sequence not in HTML inner scrolling window 

    var pvalue = numbers[i].nextSibling.nextSibling.innerHTML;
    var far = numbers[i].nextSibling.nextSibling.nextSibling.children[0].children;
    var seq_length = data.sequence_db.sequences[i].length;

    var seqObj = [];
    seqObj["length"] = seq_length;
    seqObj["pvalue"] = pvalue;
    seqObj["pn"] = [];
    seqObj["width"] = [];
    seqObj["left"] = [];
    seqObj["height"] = [];
    seqObj["color"] = [];
    seqObj["opacity"] = [];
    for (var x = 0; x < far.length; x++) {
      if ((far[x].getAttribute("style") != null)
        && (
          ( $("rdo_sites_only").checked && ! far[x].getAttribute("class").includes("scanned"))
          || ( $("rdo_sites_and_scan").checked || $("rdo_all_seqs").checked )
        )
      ) {
        if (far[x].getAttribute("style").includes("rgb")) {
          var compStyles = far[x].style;
	  // Make scanned sites get displayed first so they will not "cover" regular sites.
	  var site_pn = far[x].getAttribute("class").includes("top") ? "+" : "-";
	  var site_width = parseFloat(compStyles.width.slice(0, -1));
	  var site_left = parseFloat(compStyles.left.slice(0, -1));
	  var site_height = parseFloat(compStyles.height.slice(0, -2));
	  var site_color = compStyles.backgroundColor.slice(4, -1).replace(/ /g, "");
          if (far[x].getAttribute("class").includes("scanned")) {
	    seqObj["pn"].unshift(site_pn);
	    seqObj["width"].unshift(site_width);
	    seqObj["left"].unshift(site_left);
	    seqObj["height"].unshift(site_height);
	    seqObj["color"].unshift(site_color);
            seqObj["opacity"].unshift(0.3);
          } else {
	    seqObj["pn"].push(site_pn);
	    seqObj["width"].push(site_width);
	    seqObj["left"].push(site_left);
	    seqObj["height"].push(site_height);
	    seqObj["color"].push(site_color);
            seqObj["opacity"].push(1);
          }
	  visible_motifs[far[x].getAttribute("data-colour-index")] = site_color;
        }
      }
    }
    // Save the sequence data if it has motifs (or rdo_all_seqs is checked)
    if ($("rdo_all_seqs").checked || seqObj["width"].length > 0) { 
      bars[seq_name] = seqObj; 
    }
  }

  // jsPDF coordinates are always in points.
  var font_size = 13;
  var nbars = Object.keys(bars).length;
  var legend_font_size = 0.8 * font_size;

  // Initialize field widths in points by measuring header text.
  var font = "bold " + font_size + "pt Helvetica, sans-serif";
  var max_name_width = MeasureText(font, "Name");
  var max_pvalue_width = MeasureText(font, "p-value");

  var max_seq_length = 0;		// in characters
  var has_complement = current_alphabet.has_complement();
  var revcomp = data.options["revcomp"];

  // Measure text of numbers, names and p-values, convert to points and save the max.
  font = font_size + "pt Helvetica, sans-serif";
  var seq_name;
  for (seq_name in bars) {
    var seq_name_width = MeasureText(font, seq_name);
    var pvalue_width = MeasureText(font, pvalue);
    var seq_length = bars[seq_name]["length"];
    if (seq_length > max_seq_length) { max_seq_length = seq_length; }
    if (seq_name_width > max_name_width) { max_name_width = seq_name_width; }
    if (pvalue_width > max_pvalue_width) { max_pvalue_width = pvalue_width; }
  }

  // Get the length in characters of the longest visible motif.
  var max_motif_length = 0;
  var motif_index, motif_length;
  for (motif_index in visible_motifs) {
    motif_length = data.motifs[motif_index].len;
    if (motif_length > max_motif_length) { max_motif_length = motif_length; }
  }

  // Sort the motif indices.
  var motif_indices = [];
  var sorted_motif_indices = [];
  for (motif_index in visible_motifs) {
    motif_indices.push(Number(motif_index));
  }
  sorted_motif_indices = motif_indices.sort(function(a, b){return a-b;});

  // Set up values for main section.
  var height = (nbars+1) * (2.6*font_size);
  var nmotifs = Object.keys(visible_motifs).length;
  var name_field_width = max_name_width + font_size;
  var pvalue_field_width = max_pvalue_width + font_size;
  var plus_minus_field_width = has_complement ? 2*font_size : font_size;
  var non_diagram_width = name_field_width + pvalue_field_width + plus_minus_field_width;
  var diagram_width = 47 * font_size;
  var pix_per_char = diagram_width/max_seq_length;
  var x_scale_factor = data.sequence_db.max_length/100;	// Scale factor comes from function make_block().
  var diagram_line_height = (height-2*font_size)/nbars; 
  var doc_width = diagram_width + non_diagram_width + 2*font_size;
  var doc_height = height + 0.5*font_size;

  // Set up values for the legend.
  var tmp_font = legend_font_size + "pt Courier, normal";
  var courier_width = MeasureText(tmp_font, "A");

  var legend_line_height = 1.2 * legend_font_size;
  var index_field_width = 3 * legend_font_size;
  var symbol_field_width = 5 * legend_font_size;
  var legend_non_consensus_width = index_field_width + symbol_field_width + 3*legend_font_size;
  var legend_hdr_font = legend_font_size + "pt Helvetica, sans-serif";
  var consensus_hdr_width = MeasureText(legend_hdr_font, "Motif Consensus");
  var consensus_field_width = doc_width - legend_non_consensus_width - legend_font_size;
  // Get number of characters that will fit in legend consensus field.
  var legend_split_length = Math.floor(consensus_field_width/courier_width);
  // Get number of lines in legend.
  var n_legend_lines = 0;
  for (motif_index in visible_motifs) {
    motif_length = data.motifs[motif_index].len;
    n_legend_lines += Math.ceil(motif_length/legend_split_length);
  }
  if (n_legend_lines > 0) { n_legend_lines += 3; }	// header line + 2*space
  var legend_width = legend_non_consensus_width + Math.min(legend_split_length, max_motif_length)*courier_width;
  var legend_height = n_legend_lines * legend_line_height;
  doc_height += legend_height + 1*font_size;

  if (make_pdf) {
    // Now create the PDF document.
    // This next line is necessary because jsPDF silently swaps width and height.
    var orient = doc_width > doc_height ? 'landscape' : 'portrait';
    doc = new jsPDF(
      {
	orientation: orient,
	unit: 'pt',
	format: [doc_width, doc_height]
      }
    );

    // Set the font size for the PDF.
    doc.setFontSize(1.33*font_size);

    // Create the header.
    var offset = font_size;
    var liney = 1.5*font_size;
    // .. Name hdr ..
    doc.setFont("Helvetica", "bold");
    doc.text("Name", offset, liney);
    offset += name_field_width;

    // p-value hdr
    doc.setFont("Helvetica", "bolditalic");
    doc.text("p", offset + font_size, liney);
    doc.setFont("Helvetica", "bold");
    doc.text("-value", offset + 2*font_size, liney);
    offset += pvalue_field_width + plus_minus_field_width;

    // Motif Location hdr
    doc.text("Motif Locations", offset, liney);

    // Generate the data object for the PDF.
    liney -= 0.5*font_size;
    var dy = font_size/3.5;
    for (var seq_name in bars) {
      liney += diagram_line_height;
      offset = font_size;

      //
      // Generate the text fields.
      //
      doc.setFont("Helvetica", "normal");

      // Sequence name text
      doc.text(seq_name, offset, liney + dy);
      offset += name_field_width;

      // p-value text
      doc.text(bars[seq_name]["pvalue"], offset + pvalue_field_width, liney + dy, {align: "right"});
      offset += pvalue_field_width;

      // +/- text (optional)
      if (has_complement) {
	doc.text("+", offset+font_size, liney + dy - font_size/2);
	if (revcomp) {
	  doc.text("-", offset+1.15*font_size, liney + dy + font_size/2);
	}
      }
      offset += plus_minus_field_width;

      // Generate the base line.
      doc.setLineWidth(0.35);
      doc.line(offset, liney, offset + (bars[seq_name]["length"] * pix_per_char), liney);

      // Generate the blocks.
      for (var i = 0; i < bars[seq_name]["width"].length; i++) {
	if (bars[seq_name]["pn"][i] == undefined) { continue; }
	rgb = bars[seq_name]["color"][i].split(",").map(Number);
	var opacity =  bars[seq_name]["opacity"][i]; 
	if (opacity != 1) { rgb = RGBAtoRGB(rgb, [255,255,255], opacity); }
	var bar_x = offset + (bars[seq_name]["left"][i] * x_scale_factor * pix_per_char);
	var bar_y = (bars[seq_name]["pn"][i] == "+") ? (liney - 0.1*font_size*bars[seq_name]["height"][i]) : liney;
	doc.setFillColor(rgb[0], rgb[1], rgb[2]);
	doc.rect(bar_x, bar_y, bars[seq_name]["width"][i] * x_scale_factor * pix_per_char, 0.1*font_size*bars[seq_name]["height"][i], 'FD');
      }
    }

    //
    // Generate the legend.
    //
    if (n_legend_lines > 0) {
      doc.setFontSize(1.33*legend_font_size);
      dy = 0.8 * legend_font_size;

      // The legend header.
      var legend_top = liney + 2*legend_font_size;
      liney += 4.5*legend_font_size;
      offset = legend_font_size;
      doc.setFont("Helvetica", "bold");
      doc.text("Motif", offset, liney);
      offset += index_field_width + legend_font_size;
      doc.text("Symbol", offset, liney);
      offset += symbol_field_width + legend_font_size;
      doc.text("Motif Consensus", offset, liney);
      liney -= 0.5*legend_font_size;
      liney += legend_line_height;

      for (var i=0; i<motif_indices.length; i++) {
	motif_index = sorted_motif_indices[i];
	offset = legend_font_size;

	// Motif Name
	doc.setFont("Helvetica", "normal");
	var motif_index_string = (motif_index+1).toString();
	motif_index_string = motif_index_string + ".";
	var dx = 3 * legend_font_size;
	doc.text(motif_index_string, offset+dx, liney+dy, {align: "right"});
	offset += index_field_width + legend_font_size;

	// Motif Symbol
	motif_length = data.motifs[motif_index].len;
	rgb = visible_motifs[motif_index].split(",").map(Number);
	var bar_x = offset;
	var bar_y = liney;
	doc.setFillColor(rgb[0], rgb[1], rgb[2]);
	doc.rect(bar_x, bar_y, symbol_field_width*(motif_length/max_motif_length), legend_font_size, 'FD');
	offset += symbol_field_width + legend_font_size;

	// Motif Consensus Sequence
	doc.setFont("Courier", "normal");
	var motif_consensus = data.motifs[motif_index].id;
	doc.text(motif_consensus, offset, liney+dy, {maxWidth: legend_split_length*courier_width});
	liney += Math.ceil(motif_length/legend_split_length) * legend_line_height;
      }

      // Draw box around legend.
      doc.rect(
	0.5*legend_font_size, 
	legend_top + 0.5*legend_font_size, 
	Math.max(legend_width, legend_non_consensus_width + consensus_hdr_width + courier_width),
	legend_height
      );
    } // legend

    doc.save('motif_locations.pdf');

  } else {
    // Download an SVG document.
    var body = d3.select("#blocks").append("svg")
      .attr("width", (diagram_width + non_diagram_width) + "pt")
      .attr("height", (doc_height+legend_font_size).toString())
      .attr("background-color", "lightgrey")
      .attr("id", "memeSVG")
      .attr("xmlns", "http://www.w3.org/2000/svg");

    // Create the header.
    var x = 0;
    var offset = font_size;
    var liney = 1.5*font_size;

    // .. Name hdr ..
    body.append("text")
      .attr("x", offset)
      .attr("y", liney)
      .attr("font-family", "Helvetica, sans-serif")
      .attr("font-size", font_size+"pt")
      .attr("font-weight", "bold")
      .text("Name");
      offset += name_field_width;

    // p-value hdr
    body.append("text")
      .attr("x", offset + 2*font_size)
      .attr("y", liney)
      .attr("font-family", "Helvetica, sans-serif")
      .attr("font-size", font_size+"pt")
      .attr("font-weight", "bold")
      .attr("font-style", "italic")
      .text("p");
    body.append("text")
      .attr("x", offset + 3*font_size)
      .attr("y", liney)
      .attr("font-family", "Helvetica, sans-serif")
      .attr("font-size", font_size+"pt")
      .attr("font-weight", "bold")
      .text("-value");
    offset += pvalue_field_width + plus_minus_field_width + font_size;

    // Motif Location hdr
    body.append("text")
      .attr("x", offset)
      .attr("y", liney)
      .attr("font-family", "Helvetica, sans-serif")
      .attr("font-size", font_size+"pt")
      .attr("font-weight", "bold")
      .text("Motif Locations");

    // Generate the data for the SVG.
    liney -= 0.5*font_size;
    var dy = font_size/3.5;
    for (var seq_name in bars) {
      liney += diagram_line_height;
      offset = font_size;

      //
      // Generate the text fields.
      //

      // Sequence name text
      body.append("text")
	.attr("x", offset)
	.attr("y", liney + dy)
	.attr("font-family", "Helvetica, sans-serif")
	.attr("font-size", font_size+"pt")
	.text(seq_name);
      offset += name_field_width;

      // p-value text
      body.append("text")
	.attr("x", offset + font_size+ pvalue_field_width)
	.attr("y", liney + dy)
	.attr("font-family", "Helvetica, sans-serif")
	.attr("font-size", font_size+"pt")
	.attr("text-anchor", "end")
	.text(bars[seq_name]["pvalue"]);
      offset += pvalue_field_width + font_size;

      // +/- text (optional)
      if (has_complement) {
	body.append("text")
	  .attr("x", offset+font_size)
	  .attr("y", liney + dy - font_size/2)
	  .attr("font-family", "Helvetica, sans-serif")
	  .attr("font-size", font_size+"pt")
	  .text("+");
	if (revcomp) {
	  body.append("text")
	    .attr("x", offset+1.15*font_size)
	    .attr("y", liney + dy + font_size/2)
	    .attr("font-family", "Helvetica, sans-serif")
	    .attr("font-size", font_size+"pt")
	    .text("-");
	}
      }
      offset += plus_minus_field_width;

      // Generate the base line.
      body.append("line")
	.attr("x1", offset)
	.attr("x2", offset + (bars[seq_name]["length"] * pix_per_char))
	.attr("y1", liney)
	.attr("y2", liney)
	.attr("stroke-width", 0.5)
	.attr("stroke","black");

      // Generate the blocks.
      for (var i = 0; i < bars[seq_name]["width"].length; i++) {
	if (bars[seq_name]["pn"][i] == undefined) { continue; }
	body.append("rect")
	  .attr("x", offset + (bars[seq_name]["left"][i] * x_scale_factor * pix_per_char) )
	  .attr("y", (bars[seq_name]["pn"][i] == "+") ? (liney - 0.1*font_size*bars[seq_name]["height"][i]) : liney)
	  .attr("width", bars[seq_name]["width"][i] * x_scale_factor * pix_per_char) 
	  .attr("height", 0.1*font_size*bars[seq_name]["height"][i])
	  .attr("fill", "rgb("+bars[seq_name]["color"][i] + ")")
	  .attr("fill-opacity", bars[seq_name]["opacity"][i])
	  .attr("stroke-width", 0.5)
	  .attr("stroke","black");
      }
    }
    
    //
    // Generate the legend.
    //
    if (n_legend_lines > 0) {
      dy = 0.8 * legend_font_size;

      // The legend header.
      var legend_top = liney + 2*legend_font_size;
      liney += 4.5*legend_font_size;
      offset = legend_font_size;
      body.append("text")
	.attr("x", offset)
	.attr("y", liney)
	.attr("font-family", "Helvetica, sans-serif")
	.attr("font-size", legend_font_size+"pt")
        .attr("font-weight", "bold")
	.text("Motif");
      offset += index_field_width + legend_font_size;
      body.append("text")
	.attr("x", offset)
	.attr("y", liney)
	.attr("font-family", "Helvetica, sans-serif")
	.attr("font-size", legend_font_size+"pt")
        .attr("font-weight", "bold")
	.text("Symbol");
      offset += symbol_field_width + legend_font_size;
      body.append("text")
	.attr("x", offset)
	.attr("y", liney)
	.attr("font-family", "Helvetica, sans-serif")
	.attr("font-size", legend_font_size+"pt")
        .attr("font-weight", "bold")
	.text("Motif Consensus");
      liney -= 0.5*legend_font_size;
      liney += legend_line_height;

      for (var i=0; i<motif_indices.length; i++) {
	motif_index = sorted_motif_indices[i];
	offset = legend_font_size;

	// Motif Name
	var motif_index_string = (motif_index+1).toString();
	motif_index_string = motif_index_string + ".";
	var dx = 3.3 * legend_font_size;
	body.append("text")
	  .attr("x", offset+dx)
	  .attr("y", liney+dy)
	  .attr("font-family", "Helvetica, sans-serif")
	  .attr("font-size", legend_font_size+"pt")
	  .attr("text-anchor", "end")
	  .text(motif_index_string);
	offset += index_field_width + legend_font_size;

	// Motif Symbol
	motif_length = data.motifs[motif_index].len;
	var bar_x = offset;
	var bar_y = liney;
	body.append("rect")
	  .attr("x", bar_x )
	  .attr("y", liney)
	  .attr("width", symbol_field_width*(motif_length/max_motif_length)) 
	  .attr("height", legend_font_size)
	  .attr("fill", "rgb("+ visible_motifs[motif_index] + ")")
	  .attr("stroke-width", 0.5)
	  .attr("stroke","black");
	offset += symbol_field_width + legend_font_size;

	// Motif Consensus Sequence
	var motif_consensus = data.motifs[motif_index].id;
	var cons_length = motif_consensus.length;
        var start_index = 0;
        while (start_index < cons_length) {
	  body.append("text")
	    .attr("x", offset)
	    .attr("y", liney+dy)
	    .attr("font-family", "Courier")
	    .attr("font-size", legend_font_size+"pt")
	    .text(motif_consensus.slice(start_index, Math.min(cons_length, start_index+legend_split_length)))
	  liney += legend_line_height;
	  start_index += legend_split_length;
        }
      }

      // Draw box around legend.
      body.append("rect")
	.attr("x", 0.5*legend_font_size)
	.attr("y", legend_top + 0.5*legend_font_size)
	.attr("width", Math.max(legend_width, legend_non_consensus_width + consensus_hdr_width + courier_width))
	.attr("height", legend_height)
	.attr("fill", "none")
	.attr("stroke-width", 1)
	.attr("stroke", "black");

    } // legend

    var svg = document.getElementsByTagName("svg")[0].outerHTML;
    var svgBlob = new Blob([svg], {type:"image/svg+xml;charset=utf-8"});
    var svgUrl = URL.createObjectURL(svgBlob);
    var downloadLink = document.createElement("a");
    downloadLink.href = svgUrl;
    downloadLink.download = "meme-motif-locations.svg";
    document.getElementById("sites_sec").appendChild(downloadLink);
    downloadLink.click();
    downloadLink.remove();
    document.getElementById("memeSVG").remove();

  } // SVG

};

function make_blocks() {
  "use strict";
  function add_seqs_filter(container, id, checked, label_text, help_topic) {
    "use strict";
    var label, radio;
    radio = document.createElement("input");
    radio.type = "radio";
    radio.name = "seqs_display";
    radio.id = id;
    radio.checked = checked;
    radio.addEventListener('click', action_seqs_filter, false);
    label = document.createElement("label");
    label.appendChild(document.createTextNode(label_text));
    label.htmlFor = id;
    container.appendChild(radio);
    container.appendChild(label);
    if (help_topic) {
      container.appendChild(document.createTextNode("\xA0"));
      container.appendChild(help_button(help_topic));
    }
  }
  function add_block_diagram_button(container, id, buttonText, help_topic) {
    var button, label;
    button = document.createElement("button");
    button.id = id;
    label = document.createTextNode(buttonText);
    button.appendChild(label);
    button.onclick = (id === "pdfButton") ? download_PDF_block_diagram : download_SVG_block_diagram;
    container.appendChild(document.createTextNode("  "));
    container.appendChild(button);
    if (help_topic) {
      container.appendChild(document.createTextNode("\xA0"));
      container.appendChild(help_button(help_topic));
    }
    //var new_icon = document.createElement("img");
    //new_icon.src = new_icon_src;
    //new_icon.alt = "NEW";
    //container.appendChild(document.createTextNode("  "));
    //container.appendChild(new_icon);
  }
  function add_blocks_header(row, title, nopad, help_topic) {
    "use strict";
    var div, divcp, th;
    th = document.createElement("th");
    div = document.createElement("div");
    div.className = "blocks_th_inner";
    if (typeof title !== "object") {
      title = document.createTextNode("" + title);
    }
    div.appendChild(title);
    if (help_topic) {
      div.appendChild(document.createTextNode("\xA0"));
      div.appendChild(help_button(help_topic));
    }
    divcp = div.cloneNode(true);
    divcp.className = "blocks_th_hidden";
    th.appendChild(div);
    th.appendChild(divcp);
    if (nopad) {
      th.className = "nopad";
    }
    row.appendChild(th);
  }
  var container;
  var page, view_height, outer_tbl, inner_tbl, tbl, thead, tbody, rhead;
  var in_view, i, seq_count;
  
  page = (document.compatMode === "CSS1Compat") ? document.documentElement : document.body;
  view_height = Math.max(page.clientHeight - 300, 300);

  container = $("blocks");
  toggle_class(container, "hide_empty_seqs", true);
  toggle_class(container, "hide_only_scan", true);
  container.innerHTML = "";
  add_seqs_filter(container, "rdo_sites_only", true, "Only Motif Sites", "pop_motif_sites");
  add_seqs_filter(container, "rdo_sites_and_scan", false, "Motif Sites+Scanned Sites", "pop_scanned_sites");
  add_seqs_filter(container, "rdo_all_seqs", false, "All Sequences", "pop_all_sequences");
  add_block_diagram_button(container, "pdfButton", "Download PDF", "pop_download_pdf_motif_locations");
  add_block_diagram_button(container, "svgButton", "Download SVG", "pop_download_svg_motif_locations");

  outer_tbl = document.createElement("div");
  outer_tbl.className = "blocks_outer";

  inner_tbl = document.createElement("div");
  inner_tbl.id = "blocks_scroll";
  inner_tbl.className = "blocks_inner";
  inner_tbl.style.maxHeight = view_height + "px";
  outer_tbl.appendChild(inner_tbl);

  tbl = document.createElement("table");
  tbl.className = "blocks_tbl";
  inner_tbl.appendChild(tbl);

  thead = document.createElement("thead");
  tbl.appendChild(thead);
  tbody = document.createElement("tbody");
  tbl.appendChild(tbody);

  rhead = thead.insertRow(thead.rows.length);
  add_blocks_header(rhead, "", true);
  add_blocks_header(rhead, "Name", false, "pop_seq_name");
  add_blocks_header(rhead, make_pv_text(), false, "pop_seq_pvalue");
  add_blocks_header(rhead, "Motif Locations", false, "pop_motif_location");

  container.appendChild(outer_tbl);

  seq_count = data["sequence_db"]["sequences"].length;
  in_view = Math.max(Math.ceil(view_height / 25), 1);
  i = append_blocks_entries(tbody, 0, in_view);

  while (i < seq_count && inner_tbl.scrollHeight - (inner_tbl.scrollTop + inner_tbl.offsetHeight) < 400) {
    i = append_blocks_entries(tbody, i, 20);
  }
  inner_tbl.data_blocks_index = i;
  if (i < seq_count) {
    inner_tbl.addEventListener('scroll', make_blocks_entries, false);
  }
}

function make_scan_popup(site, motif) {
  return function (e) {
    "use strict";
    var pop, xy, padding, edge_padding, pop_left, pop_top, page_width;
    var lflank, match, rflank, pspm;
    if (!e) var e = window.event;
    pop = make_scan_popup.pop;
    if (e.type === "mouseover") {
      if (pop) return;
      pop = clone_template("tmpl_scan_info");
      pspm = new Pspm(motif.pwm);
      if (site.rc) pspm.reverse_complement(current_alphabet);
      set_tvar(pop, "tvar_logo", make_small_logo(current_alphabet, pspm, {"className": "scan_logo"}));
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", site.pvalue.toExponential(2));
      set_tvar(pop, "tvar_start", site.pos + 1);
      set_tvar(pop, "tvar_end", site.pos + motif.len);

      document.body.appendChild(pop);
      position_popup(this, pop);
      make_scan_popup.pop = pop;
    } else if (e.type === "mouseout") {
      if (pop) {
        pop.parentNode.removeChild(pop);
        make_scan_popup.pop = null;
      }
    }
  };
}

function make_block_popup(site, motif, block) {
  return function (e) {
    "use strict";
    var pop;
    var lflank, match, rflank, pspm, ruler, match_seq, match_width;
    if (!e) var e = window.event;
    pop = make_block_popup.pop;
    if (e.type === "mouseover") {
      if (pop) return;
      pop = clone_template("tmpl_block_info");
      pspm = new Pspm(motif.pwm);
      if (site.rc) { // must be dna
        pspm.reverse_complement(current_alphabet);
        lflank = current_alphabet.invcomp_seq(site.rflank);
        match = current_alphabet.invcomp_seq(site.match);
        rflank = current_alphabet.invcomp_seq(site.lflank);
      } else {
        lflank = site.lflank;
        match = site.match;
        rflank = site.rflank;
      }
      ruler = document.getElementById("measure_match");
      match_seq = make_seq(current_alphabet, match);
      ruler.innerHTML = "";
      ruler.appendChild(match_seq);
      match_width = ruler.clientWidth;
      ruler.removeChild(match_seq);
      set_tvar(pop, "tvar_lflank", lflank);
      set_tvar(pop, "tvar_match", match_seq);
      set_tvar(pop, "tvar_rflank", rflank);
      set_tvar(pop, "tvar_logo_pad", lflank);
      set_tvar(pop, "tvar_logo", make_small_logo(current_alphabet, pspm, {"width": match_width}));
      set_tvar(pop, "tvar_motif", motif.id);
      set_tvar(pop, "tvar_pvalue", site.pvalue.toExponential(2));
      set_tvar(pop, "tvar_start", site.pos + 1);
      set_tvar(pop, "tvar_end", site.pos + motif.len);

      document.body.appendChild(pop);
      position_popup(block, pop);
      make_block_popup.pop = pop;
    } else if (e.type === "mouseout") {
      if (pop) {
        pop.parentNode.removeChild(pop);
        make_block_popup.pop = null;
      }
    }
  };
}

//
// action_show_more
//
// Show more information on the motif.
///
function action_show_more(e) {
  var node, tr, tbody, table, container, motif, ordinal;
  var expanded_motif;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  // find the row that contains the cell
  node = this;
  do {
    if (node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find row!?");
  tr = node;
  // get info
  motif = tr.data_motif;
  ordinal = tr.data_ordinal;
  // find tbody
  do {
    if (node.tagName === "TBODY") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find tbody!?");
  tbody = node;
  // find table
  do {
    if (node.tagName === "TABLE") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find table!?");
  table = node;
  // find container
  container = node.parentNode;
  // make a expanded motif
  motif["expanded"] = true;
  expanded_motif = make_expanded_motif(current_alphabet, ordinal, motif, 
      container.data_more_x, container.data_submit_x);
  // now determine how to place it
  if (tbody.rows.length === 1) {
    // only us in the table so the table can be replaced
    container.replaceChild(expanded_motif, table);
  } else if (tbody.rows[0] === tr) {
    // first row, so remove and insert an expanded motif before
    table.deleteRow(tr.rowIndex);
    container.insertBefore(expanded_motif, table);
  } else if (tbody.rows[tbody.rows.length - 1] === tr) {
    // last row, so remove and insert an expanded motif after
    table.deleteRow(tr.rowIndex);
    container.insertBefore(expanded_motif, table.nextSibling);
  } else {
    var table2, tbody2;
    table2 = table.cloneNode(false);
    table2.appendChild(table.tHead.cloneNode(true));
    tbody2 = table.tBodies[0].cloneNode(false);
    table2.appendChild(tbody2);
    container.insertBefore(table2, table.nextSibling);
    for (i = tbody.rows.length - 1; i >= 0; i--) {
      row = tbody.rows[i];
      row.parentNode.removeChild(row);
      if (row === tr) {
        break;
      }
      tbody2.insertBefore(row, tbody2.rows[0]);
    }
    container.insertBefore(expanded_motif, table2);
  }
  find_child(expanded_motif, "tvar_less").focus();
}

//
// action_show_less
//
// Show less information on the motif.
///
function action_show_less(e) {
  var btn;
  var expanded_motif, container, motif, ordinal, colw, focus_target;
  var table, tbody, tbody2, row, table_before, table_after;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  btn = this;
  // find expanded motif
  expanded_motif = find_parent(btn, "expanded_motif");
  if (!expanded_motif) throw new Error("Expected expanded motif.");
  // find the container
  container = expanded_motif.parentNode;
  // get data
  motif = expanded_motif.data_motif;
  ordinal = expanded_motif.data_ordinal;
  colw = container.data_colw;
  // get the table before
  table_before = expanded_motif.previousSibling;
  if (table_before && table_before.tagName !== "TABLE") {
    table_before = null;
  }
  // get the table after
  table_after = expanded_motif.nextSibling;
  if (table_after && table_after.tagName !== "TABLE") {
    table_after = null;
  }
  // see if there is a table below or above that we can put this in.
  // if there is a table both below and above then add this motif and
  // all ones below to the above table
  motif["expanded"] = false;
  if (table_before && table_after) {
    tbody = table_before.tBodies[0];
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, current_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
    tbody2 = table_after.tBodies[0];
    while (tbody2.rows.length > 0) {
      row = tbody2.rows[0];
      row.parentNode.removeChild(row);
      tbody.appendChild(row);
    }
    container.removeChild(table_after);
  } else if (table_before) {
    tbody = table_before.tBodies[0];
    row = tbody.insertRow(tbody.rows.length);
    make_motif_table_entry(row, current_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
  } else if (table_after) {
    tbody = table_after.tBodies[0];
    row = tbody.insertRow(0);
    make_motif_table_entry(row, current_alphabet, ordinal, motif, colw);
    focus_target = find_child(row.cells[5], "sym_btn");
    container.removeChild(expanded_motif);
  } else {
    //no table above or below!
    // make a new table
    table = make_motifs_table(current_alphabet, ordinal, [motif], colw, data["stop_reason"]);
    focus_target = find_child(table.tBodies[0].rows[0].cells[5], "sym_btn");
    container.replaceChild(table, expanded_motif);
  }
  focus_target.focus();
}

//FIXME -- can we delete this junk?
//function action_show_outpop(e) {
function fred_action_show_outpop(e) {
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
  var node;
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
  // load the motif logo
  node = this;
  do {
    if (/\bexpanded_motif\b/.test(node.className) || node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (node === null) throw new Error("Expected node!");
  update_outpop_motif(node.data_ordinal - 1);
  // display the download popup
  $("grey_out_page").style.display = "block";
  $("download").style.display = "block";
  $("outpop_close").focus();
} // fred_action_show_outpop

function action_btn_rc(e) {
  "use strict";
  var node, tr, motif, box, logo_box, tab_st, tab_rc, rc;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  node = this;
  do {
    if (node.tagName === "TR") break;
  } while (node = node.parentNode);
  if (!node) throw new Error("Expected to find row!?");
  tr = node;
  // get info
  motif = tr.data_motif;
  box = find_parent(this, "preview_box");
  logo_box = find_child(box, "preview_logo_box");
  tab_st = find_child(box, "plus");
  tab_rc = find_child(box, "minus");
  rc = (this === tab_rc);
  motif["rc"] = rc;
  toggle_class(logo_box, "show_rc_logo", rc);
  toggle_class(tab_st, "active", !rc);
  toggle_class(tab_rc, "active", rc);
}

function action_rc_tab(e) {
  "use strict";
  var box, logo_box, tab_st, tab_rc, rc;
  if (!e) e = window.event;
  if (e.type === "keydown") {
    if (e.keyCode !== 13 && e.keyCode !== 32) {
      return;
    }
    // stop a submit or something like that
    e.preventDefault();
  }
  box = find_parent(this, "expanded_motif");
  logo_box = find_child(box, "tvar_logo");
  tab_st = find_child(box, "tvar_tab");
  tab_rc = find_child(box, "tvar_tab_rc");
  rc = (this === tab_rc);
  box.data_motif["rc"] = rc;
  toggle_class(logo_box, "show_rc_logo", rc);
  toggle_class(tab_st, "activeTab", !rc);
  toggle_class(tab_rc, "activeTab", rc);
}

function action_seqs_filter() {
  "use strict";
  var block_container;
  block_container = $("blocks");
  if ($("rdo_all_seqs").checked) {
    toggle_class(block_container, "hide_empty_seqs", false);
    toggle_class(block_container, "hide_only_scan", false);
  } else if ($("rdo_sites_and_scan").checked) {
    toggle_class(block_container, "hide_empty_seqs", true);
    toggle_class(block_container, "hide_only_scan", false);
  } else if ($("rdo_sites_only").checked) {
    toggle_class(block_container, "hide_empty_seqs", true);
    toggle_class(block_container, "hide_only_scan", true);
  }
}

//
// page_loaded
//
// Called when the page has loaded for the first time.
///
function page_loaded() {
  post_load_setup();
}

//
// page_loaded
//
// Called when a cached page is reshown.
///
function page_shown(e) {
  if (e.persisted) post_load_setup();
}

//
// page_loaded
//
// Called when the page is resized
///
function page_resized() {
  var page, blocks_scroll;
  update_scroll_pad();
  page = (document.compatMode === "CSS1Compat") ? document.documentElement : document.body;
  blocks_scroll = $("blocks_scroll");
  if (blocks_scroll) {
    blocks_scroll.style.maxHeight = Math.max(page.clientHeight - 300, 300) + "px";
  }
}

//
// pre_load_setup
//
// Run before the page is displayed
///
function pre_load_setup() {
  var start, hue, sat, light, divisions;
  var i, j, motifs, motif, sites, site, sequences, sequence;
  var max_seq_len;
  motifs = data["motifs"];
  sequences = data["sequence_db"]["sequences"];
  max_seq_len = 1;
  if (sequences) {		// no sequences if -brief
    for (i = 0; i < sequences.length; i++) {
      sequence = sequences[i];
      sequence["sites"] = [];
      if (sequence["length"] > max_seq_len) {
	max_seq_len = sequence["length"];
      }
   }
  }
  data["sequence_db"]["max_length"] = max_seq_len;
  // use hsl colours
  start = 0; //red
  sat = 100;
  light = 50;
  for (i = 0; i < motifs.length; i++) {
    motif = motifs[i];
    // give the motif a colour
    divisions = 1 << Math.ceil(Math.log(i + 1) / Math.LN2);
    hue = start + (360 / divisions) * ((i - (divisions >> 1)) * 2 + 1);
    motif["colour"] = "hsl(" + hue + ", " + sat + "%, " + light + "%)";
    // associate sites with sequences as well 
    // to make generating the block diagram easier
    sites = motif["sites"];
    for (j = 0; j < sites.length; j++) {
      site = sites[j];
      sequence = sequences[site["seq"]];
      // record the motif index
      site["motif"] = i;
      // add the site to the sequence
      sequence["sites"].push(site);
    }
  }
}

//
// post_load_setup
//
// Run when the page has loaded, or been reloaded.
//
function post_load_setup() {
  update_scroll_pad();
  if (data["motifs"].length > 0) {
    make_motifs();
    if (data.sequence_db.primary_count > data.options.brief) {
      if (data.options.brief == 1000) {
        $("blocks").innerHTML = "<p>Output of sites suppressed because there were more than 1000 (primary) sequences.</p>";
      } else {
        $("blocks").innerHTML = "<p>Output of motif locations suppressed by -brief option.</p>";
      }
    } else {
      make_blocks();
    }
  } else {
    $("motifs").innerHTML = "<p>No significant motifs found!</p>"; // clear content
    $("motifs").innerHTML += "<p><b>" + data["stop_reason"] + "</b></p>";
    $("blocks").innerHTML = "<p>No significant motifs found!</p>";
  }
}

pre_load_setup();
