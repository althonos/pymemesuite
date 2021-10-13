var current_program = "DREME";
var current_motif = 0;
var current_alphabet = new Alphabet(data.alphabet);

/*
 * Make a colorized sequence.
 */
function make_seq(seq) {
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
    lbox.style.color = current_alphabet.get_colour(current_alphabet.get_index(letter));
    lbox.appendChild(document.createTextNode(seq.substring(i, j)));
    sbox.appendChild(lbox);
  }
  return sbox;
}

/*
 * Make the table header for the discovered motifs.
 */
function make_motif_header() {
  var row = document.createElement("tr");
  add_text_header_cell(row, "", "", "motif_ordinal");
  add_text_header_cell(row, "Motif", "pop_motifs_word", "motif_word");
  add_text_header_cell(row, "Logo", "pop_motifs_logo", "motif_logo");
  if (data.options.revcomp) {
    add_text_header_cell(row, "RC Logo", "pop_motifs_rc_logo", "motif_logo");
  }
  add_text_header_cell(row, "E-value", "pop_motifs_evalue", "motif_evalue");
  add_text_header_cell(row, "Unerased E-value", "pop_motifs_uevalue", "motif_evalue");
  add_text_header_cell(row, "More", "pop_more", "motif_more");
  add_text_header_cell(row, "Submit/Download", "pop_submit_dl", "motif_submit");
  row.className = "more";
  return row;
}

/*
 * Make a compact motif summary row for the discovered motifs.
 */
function make_motif_row(tbody, ordinal, motif) {
  var row = document.createElement("tr");
  add_text_cell(row, "" + ordinal + ".", "motif_ordinal");
  add_text_cell(row, motif["id"], "motif_word");
  var pspm = new Pspm(motif["pwm"]);
  add_cell(row, make_logo(current_alphabet, pspm, 50, false, 0, "normal_logo"), "motif_logo");
  if (data.options.revcomp) {
    add_cell(row, make_logo(current_alphabet, pspm, 50, true, 0, "flipped_logo"), "motif_logo");
  }
  add_text_cell(row, motif["evalue"], "motif_evalue");
  add_text_cell(row, motif["unerased_evalue"], "motif_evalue");
  add_cell(row, make_sym_btn(text_pair("\u21A7", "less", "\u21A5", "more"), "Show more information.", function(e) { toggle_class(tbody, "collapsed"); }, "\u21A5", ""), "motif_more");
  add_cell(row, make_sym_btn("\u21E2", "Submit the motif to another MEME Suite program or download it.", function(e) { action_show_outpop(e, ordinal); }), "motif_submit");
  return row;
}

/*
 * Make a sortable table of enriched matching rows.
 */
function make_motif_words(motif) {
  var row, i, match;
  var table = document.createElement("table");
  var thead = document.createElement("thead");
  row = document.createElement("tr");
  add_text_header_cell(row, "Word", "pop_match_word", "match_word", function(e) {sort_table(this, compare_words);});
  add_text_header_cell(row, "Positives", "pop_match_pos", "match_count", function(e) {sort_table(this, compare_counts);});
  add_text_header_cell(row, "Negatives", "pop_match_neg", "match_count", function(e) {sort_table(this, compare_counts);});
  add_text_header_cell(row, "P-value", "pop_match_pval", "match_evalue", function(e) {sort_table(this, compare_evalues);});
  add_text_header_cell(row, "E-value", "pop_match_eval", "match_evalue", function(e) {sort_table(this, compare_evalues);});
  thead.appendChild(row);
  table.appendChild(thead);
  var tbody = document.createElement("tbody");
  for (i = 0; i < motif.matches.length; i++) {
    match = motif.matches[i];
    row = document.createElement("tr");
    add_cell(row, make_seq(match.seq), "match_word");
    add_text_cell(row, match.p + " / " + data.sequence_db.count, "match_count");
    add_text_cell(row, match.n + " / " + data.control_db.count, "match_count");
    add_text_cell(row, match.pvalue, "match_evalue");
    add_text_cell(row, match.evalue, "match_evalue");
    tbody.appendChild(row);
  }
  table.appendChild(tbody);
  return table;
}

/*
 * Make an expanded view of a discovered motif.
 */
function make_motif_exp(tbody, ordinal, motif) {
  "use strict";
  var box, pspm, logo_box;
  box = $("tmpl_motif_expanded").cloneNode(true);
  toggle_class(box, "template", false);
  box.id = "";
  var pspm = new Pspm(motif["pwm"]);
  find_child(box, "tvar_logo").appendChild(make_logo(current_alphabet, pspm, 150, false, 0, "normal_logo"));
  if (data.options.revcomp) {
    find_child(box, "tvar_rclogo").appendChild(make_logo(current_alphabet, pspm, 150, true, 0, "flipped_logo"));
  }
  set_tvar(box, "tvar_p", motif["p"]);
  set_tvar(box, "tvar_p_total", data.sequence_db.count);
  set_tvar(box, "tvar_n", motif["n"]);
  set_tvar(box, "tvar_n_total", data.control_db.count);
  set_tvar(box, "tvar_pvalue", motif["pvalue"]);
  set_tvar(box, "tvar_evalue", motif["evalue"]);
  set_tvar(box, "tvar_uevalue", motif["unerased_evalue"]);
  set_tvar(box, "tvar_words", make_motif_words(motif));
  var cell = document.createElement("td");
  cell.colSpan = 8;
  cell.appendChild(box);
  var row = document.createElement("tr");
  row.className = "more";
  row.appendChild(cell);
  return row;
}

/*
 * Convert a string containing a scientific number into the log of that number
 * without having an intermediate representation of the number.
 * This is intended to avoid underflow problems with the tiny evalues that
 * MEME and DREME can create.
 */
function sci2log(scinum) {
  "use strict";
  var ev_re, match, sig, exp;
  ev_re = /^(.*)e(.*)$/;
  if (match = ev_re.exec(scinum)) {
    sig = parseFloat(match[1]);
    exp = parseInt(match[2]);
    return Math.log(sig) + (exp * Math.log(10));
  }
  return 0;
}

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
  container.appendChild(table);
}

/*
 * Sort all rows in the first table body based on the column of the given element and the comparison function.
 * The sort is not very fast and is intended for small tables only.
 */
function sort_table(colEle, compare_function) {
  //find the parent of colEle that is either a td or th
  var i, j;
  var cell = colEle;
  while (true) {
    if (cell == null) return;
    if (cell.nodeType == Node.ELEMENT_NODE && 
        (cell.tagName.toLowerCase() == "td" || cell.tagName.toLowerCase() == "th")) {
      break;
    }
    cell = cell.parentNode;
  }
  //find the parent of cell that is a tr
  var row = cell;
  while (true) {
    if (row == null) return;
    if (row.nodeType == Node.ELEMENT_NODE && row.tagName.toLowerCase() == "tr") {
      break;
    }
    row = row.parentNode;
  }
  //find the parent of row that is a table
  var table = row;
  while (true) {
    if (table == null) return;
    if (table.nodeType == Node.ELEMENT_NODE && table.tagName.toLowerCase() == "table") {
      break;
    }
    table = table.parentNode;
  }
  var column_index = cell.cellIndex;
  // do a bubble sort, because the tables are so small it doesn't matter
  var change;
  var trs = table.tBodies[0].getElementsByTagName('tr');
  var already_sorted = true;
  var reverse = false;
  while (true) {
    do {
      change = false;
      for (i = 0; i < trs.length -1; i++) {
        var v1 = elem_text(trs[i].cells[column_index]);
        var v2 = elem_text(trs[i+1].cells[column_index]);
        if (reverse) {
          var tmp = v1;
          v1 = v2;
          v2 = tmp;
        }
        if (compare_function(v1, v2) > 0) {
          exchange(trs[i], trs[i+1], table);
          change = true;
          already_sorted = false;
          trs = table.tBodies[0].getElementsByTagName('tr');
        }
      }
    } while (change);
    if (reverse) break;// we've sorted twice so exit
    if (!already_sorted) break;// sort did something so exit
    // when it's sorted one way already then sort the opposite way
    reverse = true;
  }
  // put arrows on the headers
  var headers = table.tHead.getElementsByTagName('tr');
  for (i = 0; i < headers.length; i++) {
    for (j = 0; j < headers[i].cells.length; j++) {
      var cell = headers[i].cells[j];
      var arrows = cell.getElementsByClassName("sort_arrow");
      var arrow;
      if (arrows.length == 0) {
        arrow = document.createElement("span");
        arrow.className = "sort_arrow";
        cell.insertBefore(arrow, cell.firstChild);
      } else {
        arrow = arrows[0];
      }
      arrow.innerHTML = "";
      if (j == column_index) {
        arrow.appendChild(document.createTextNode(reverse ? "\u25B2" : "\u25BC"));
      }
    }
  }
}

/*
 * Swap two rows in a table.
 */
function exchange(oRowI, oRowJ, oTable) {
  var i = oRowI.rowIndex;
  var j = oRowJ.rowIndex;
   if (i == j+1) {
    oTable.tBodies[0].insertBefore(oRowI, oRowJ);
  } if (j == i+1) {
    oTable.tBodies[0].insertBefore(oRowJ, oRowI);
  } else {
    var tmpNode = oTable.tBodies[0].replaceChild(oRowI, oRowJ);
    if(typeof(oRowI) != "undefined") {
      oTable.tBodies[0].insertBefore(tmpNode, oRowI);
    } else {
      oTable.appendChild(tmpNode);
    }
  }
}

/*
 * Compare two E-values which may be very small.
 */
function compare_evalues(v1, v2) {
  var e1 = sci2log(v1);
  var e2 = sci2log(v2);
  if (e1 < e2) return -1;
  else if (e1 > e2) return 1;
  return 0;
}

/*
 * Compare two counts.
 */
function compare_counts(v1, v2) {
  var re = /(\d+)\s*\/\s*\d+/;
  var m1 = re.exec(v1);
  var m2 = re.exec(v2);
  if (m1 == null && m2 == null) return 0;
  if (m1 == null) return -1;
  if (m2 == null) return 1;
  return parseInt(m2[1]) - parseInt(m1[1]);
}

/*
 * Compare two sequence words.
 */
function compare_words(v1, v2) {
  return v1.localeCompare(v2);
}

