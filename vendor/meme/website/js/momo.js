var sequences = null;

window.onpageshow = function(event) {
  if (event.persisted) {
    window.location.reload() 
  }
}

window.onload = function() {
  var i;
  var num_psms = sessionStorage.getItem("num_psms");
  if (! num_psms || num_psms==0) {
    sessionStorage.setItem("num_psms", 1);
    sessionStorage.setItem("psm_1_format", "unknown");
    num_psms = 1;
  }
  for (i = 1; i <= num_psms; i++) {
    var psm_row = "psm_row_" + i;
    $(psm_row).style.display = '';
    var psm_format = "psm_" + i + "_format";
    $(psm_format).innerHTML = sessionStorage.getItem(psm_format);
    if (i==1) { updateColumnName($(psm_format).innerHTML); }
  }
  check_column_name($("psm_column_name").value);
} // window.onload

function register_component(id, element, controler) {
  "use strict";
  if (id == "sequences") {
    sequences = controler;
  }
}

function on_sequences_change() {
  if (sequences.source.value == "noseq") $("db_background").checked = false;
  on_alg_change();
}

function on_alg_change(partial) {
  $("flanking_sequences").style.display = (($('alg_mtfx').checked || $('alg_modl').checked) ? 'block' : 'none');
  $("db_bkg").style.display = ((sequences.source.value != "noseq" && ($('alg_mtfx').checked || $('alg_modl').checked))
    ? 'block' : 'none');
  $("motifx_thresholds").style.display = ($('alg_mtfx').checked ? 'block' : 'none');
  $("modl_thresholds").style.display = ($('alg_modl').checked ? 'block' : 'none');
  if (partial != 1) {$("occurs").value = ($('alg_mtfx').checked ? 20 : 5);}
}

function add_ptm() {
  "use strict";
  var num_psm_rows = $("psm_rows").rows.length;
  var old_num = sessionStorage.getItem("num_psms");

  // Check that another row is available in the table.
  if (old_num == num_psm_rows) {
    alert("Sorry-- you can only enter up to " + num_psm_rows + " PTM files.");
    return false;
  }

  // Check that previous row has an uploaded PSM.
  var previousFile = $("psm_" + old_num);
  if (! previousFile.files[0]) {
    alert("Please select a file using the previous Browse button first.");
    return(false);
  }

  // Display the new row and record the new number of psms.
  var new_num = Number(old_num) + 1;
  var new_row = "psm_row_" + new_num;
  $(new_row).style.display = '';
  sessionStorage.setItem("num_psms", new_num);

  return(true);
} // add_ptm

// Guess PTM file type based on header and return it.
function guess_ptm_file_type(ptm_file) {
  var header = sessionStorage.getItem(ptm_file + "_line0").split('\t');
  var nfields = header.length;
  var old_psm_column_name = sessionStorage.getItem("psm_column_name");
  var filetype = {};
  if (/^>/.test(header[0])) {
    filetype.format = "FASTA";
    filetype.psm_column_name= "";
  } else if (nfields == 1) {
    filetype.format = "Raw"
    filetype.psm_column_name= "";
  } else {
    var ms_found = false;
    var peptide_found = false;
    var sequence_found = false;
    var old_name_found = false;
    for (i=0; i<nfields; i++) {
      // See if column name is one of the known ones.
      hdr = trim_whitespace(header[i]);
      if (hdr == 'modified sequence') {
	ms_found = true;
      } else if (hdr == 'Peptide') {
	peptide_found = true;
      } else if (hdr == 'sequence') {
	sequence_found = true;
      } else if (old_psm_column_name && hdr == old_psm_column_name) {
        old_name_found = true;
      }
    }
    if (ms_found == true) {			// comet
      filetype.format = 'comet';
      filetype.psm_column_name = 'modified sequence';
    } else if (peptide_found == true) {		// ms-gf+
      filetype.format = 'ms-gf+';
      filetype.psm_column_name = 'Peptide';
    } else if (sequence_found == true) {	// tide/percolator
      filetype.format = 'tide/percolator';
      filetype.psm_column_name = 'sequence';
    } else {					// unknown PSM format
      filetype.format = 'unknown PSM format';
      if (old_name_found) { 
        filetype.psm_column_name = old_psm_column_name;
      } else {
        filetype.psm_column_name = '';
      }
    }
  }
  return filetype;
} // guess_ptm_file_type

function trim_whitespace(string) {
  // Trim leading and trailing whitespace.
  string = string.replace(/^\s+/g, '');	// trim leading spaces
  string = string.replace(/\s+$/g, '');	// trim trailing spaces
  return string;
} // trim_whitespace

function check_column_name(psm_column_name) {

  // Trim whitespace from the column name before checking.
  psm_column_name = trim_whitespace(psm_column_name);

  // Check that all the PTMs contain the named column (or are prealigned if column name is blank).
  var n_bad = 0;
  var bad_files = '';
  var num_rows = sessionStorage.getItem("num_psms");
  for (i=1; i<=num_rows; i++) {
    var psm_name = "psm_" + i;
    var files = $(psm_name).files;
    if (files == null || files[0] == null) { 			// No PTM file loaded on this line. 
      if (i==1) {
        $("psm_1_format").innerHTML = '';
        $("psm_column_name").style.display = 'none';
      }
      break;
    }
    var file_name = files[0].name;
    if (psm_column_name == '') {				// name is blank
      var format = sessionStorage.getItem(psm_name + "_format");
      if (! (format == "FASTA" || format == "Raw")) {
        if (i==1) {		// First file must be prealigned if column name is blank.
	 alert("You cannot specify a blank 'Modified Peptide Column Name' with '" + format + "' PTM files.");
	 return false;
        } else {
	  bad_files = bad_files + '\n\t\t' + file_name;
	  n_bad++;
        }
      }
    } else {							// name is not blank
      var line0 = sessionStorage.getItem(psm_name + "_line0");
      var header = line0.split('\t');
      var num_fields = header.length;
      for (j=0; j<num_fields; j++) {
	if (trim_whitespace(header[j]) == psm_column_name) break;
      }
      if (j == num_fields) {
	bad_files = bad_files + '\n\t\t' + file_name;
	n_bad++;
      }
    }
  }

  // Report any bad files.
  if (n_bad > 0) {
    var file = (n_bad == 1) ? 'file' : 'files';
    var does = (n_bad == 1) ? 'does' : 'do';
    var is = (n_bad == 1) ? 'is' : 'are';
    if (psm_column_name != '') {
      alert("The PTM " + file + ":" + bad_files + 
	"\n" + does + " not contain a column named:\n\t\t" +
	"'" + psm_column_name + "'" +
	"\nPlease enter a new Modified Peptide Column Name or remove/change the PTM " + file + ".");
    } else {
      alert("The PTM " + file + ":" + bad_files + 
	"\n" + is + " not in 'FASTA' or 'Raw' format but the first file is." +
	"\nPlease enter a new first file or remove/change the listed PTM " + file + ".");
    }
    return false;
  } else {
    return true;
  }
} // check_column_name

// Check that all files have the named column.
function on_column_name_change() {
  var psm_column_name_old = $("psm_column_name").oldvalue;	// local copy of old value
  var psm_column_name_new = $("psm_column_name").value;		// local copy of new value
  var success = false;

  // Check that all files have the specified column.
  if (! check_column_name(psm_column_name_new)) {
    // Restore the column name to the old name.
    $("psm_column_name").value = psm_column_name_old;
  } else {
    // Replace column name with its trimmed value.
    $("psm_column_name").value = trim_whitespace(psm_column_name_new);
    success = true;
  }

  // Save the column name.
  sessionStorage.setItem("psm_column_name", $("psm_column_name").value);

  // Display the filter if the column name was set for the first time.
  if ($("psm_column_name").oldvalue == "" && $("psm_column_name").value != "") {
    updateFilterMenu();
    updateFilterDisplay();
  }

  return success;
} // on_column_name_change

function on_filter_field_change() {
  var filter_field = $("filter_field");
  sessionStorage.setItem("filter_field", filter_field.value);
} // on_filter_field_change

function on_file_change(psm_name) {
  var filter_enable = $("filter_enable");
  var filter_field = $("filter_field");
  
  // Retrieve the first file from the FileList object
  var f = $(psm_name).files[0];

  if (f) {
    var r = new FileReader();
    r.onload = function(e) {
      var contents = e.target.result;
      var lines = contents.split('\n');
      var line0 = '';
      var line1 = '';

      // Find the first non-empty line (header).
      for (i=0; i<lines.length; i++) {
        var tmp = lines[i].replace(/\s+/g, '');
        if (tmp.length != 0) break;
      }
      if (i < lines.length) {		// in case file is empty
        line0 = lines[i];
      }

      // Find the next non-empty line (fields).
      for (++i; i<lines.length; i++) {
        var tmp = lines[i].replace(/\s+/g, '');
        if (tmp.length != 0) break;
      }
      if (i < lines.length) {		// in case file only has one line
        line1 = lines[i];
      }
      
      // Save the first two (non-blank) lines of the file in session storage.
      sessionStorage.setItem(psm_name + "_line0", line0);
      sessionStorage.setItem(psm_name + "_line1", line1);

      // Guess the file format.
      var filetype = guess_ptm_file_type(psm_name);
      sessionStorage.setItem(psm_name + "_format", filetype.format);
      $(psm_name + "_format").innerHTML = filetype.format;
   
      // Set the probable modified peptide column name and update the
      // filter menu if this is the first PTM file.
      if (psm_name == "psm_1") { 
        //if (filetype.format == "unknown PSM format" && filetype.psm_column_name == "") {
        if (filetype.format == "unknown PSM format") {
          // new first file of unknown PSM format.
	  updateColumnName(filetype.format);	// Update the modified peptide column name
	  $("psm_column_name").value = filetype.psm_column_name;
	  //$("psm_column_name").value = "";
	  updateFilterMenu();			// Update the filter menu.
        } else {
	  if (check_column_name(filetype.psm_column_name)) {
	    updateColumnName(filetype.format);	// Update the modified peptide column name
	    $("psm_column_name").value = filetype.psm_column_name;
	    updateFilterMenu();			// Update the filter menu.
	  } else {
	    if (sessionStorage.getItem("num_psms") > 1) {
	      $("psm_column_name").value = filetype.psm_column_name;
	    }
	    return false;
	  }
        }
      } else {
        if (! check_column_name($("psm_column_name").value)) {
          remove_ptm(true);			// Remove the PTM line.
          return false;
        }
      }
    }
    r.readAsText(f);
  } else {
    alert("Failed to load file");
    filter_field.options.length = 0;
    filter_enable.checked = false;
    filter_enable.disabled = (filter_field.options.length == 0);
    updateFilterDisplay();
  }
} // on_file_change

function updateColumnName(format) {
  if (format == "FASTA" || format == "Raw" || format == '') {
    $("psm_column_name").style.display = 'none';
  } else {
    $("psm_column_name").style.display = '';
  }
} // update_column_name

function updateFilterMenu() {
  // Get the first two non-blank lines from the first PTM file.
  var line0 = sessionStorage.getItem("psm_1_line0");
  var line1 = sessionStorage.getItem("psm_1_line1");

  // Update the menu if the first PTM has been read.
  if (line0) {
    var header = line0.split('\t');
    var fields = line1.split('\t');

    // Add each field as an option in the table and save in sessionStorage
    // so the data will persist after refresh/reload of the page.
    var num_fields = 0;
    var old_filter_field = sessionStorage.getItem("filter_field");
    var new_filter_field;
    filter_field.options.length = 0;
    for (i = 0; i < header.length && i < fields.length; i++) {
      if (!isNaN(fields[i])) {	// numeric field
	filter_field.options.add(new Option(header[i], header[i], false, false));
        // Use first numeric as filter field unless we find the saved field.
        if (num_fields == 0 || (old_filter_field && old_filter_field == header[i])) {
          new_filter_field = header[i];
        } 
	num_fields++;
      }
    }
    $("filter_field").value = new_filter_field;
    sessionStorage.setItem("filter_field", new_filter_field);
  }
  
  if (filter_field.options.length == 0) {
    filter_enable.checked = false;
  }
  filter_enable.disabled = (filter_field.options.length == 0);

  updateFilterDisplay();

} // updateFilterMenu

// Clear all or just the last PTM input rows (if evt != null).
function remove_ptm(evt) {
  "use strict";
  var num_rows = sessionStorage.getItem("num_psms");
  while (num_rows >= 1) {
    // Clear the filename and format, and turn off row display.
    $("psm_" + num_rows).value = $("psm_" + num_rows).defaultValue;
    var psm_format = "psm_" + num_rows + "_format";
    $(psm_format).innerHTML = '';
    sessionStorage.setItem(psm_format, '');
    if (num_rows == 1) { 
      $("psm_column_name").style.display = 'none';	// always display row 1
    } else {
      $("psm_row_" + num_rows).style.display = "none"
    }
    num_rows--;
    if (evt) break;	// clear last row only
  }
  if (num_rows < 1) num_rows = 1;
  sessionStorage.setItem("num_psms", num_rows);
} // remove_ptm

// Display or suppress the filter checkbox
// based on the type (or absence) of (first) PTM.
// Display/suppress the input fields based on checkbox.
function updateFilterDisplay() {
  var files = $("psm_1").files;
  var ptm_file_given = (files && files[0]);
  // Display checkbox if a PTM file has been read and 
  // modified peptide column name has been defined.
  if (ptm_file_given && $("psm_column_name").value != "") {
    $("filter").style.display = 'block';
    // Display filter input fields if the filter box is checked.
    var do_show = $("filter_enable").checked;
    var colheader = $("filter_table").rows[0].cells;
    var colbody = $("filter_table").rows[1].cells;
    colheader[1].style.display = do_show ? '' : 'none';
    colheader[2].style.display = do_show ? '' : 'none';
    colheader[3].style.display = do_show ? '' : 'none';
    colbody[1].style.display = do_show ? '' : 'none';
    colbody[2].style.display = do_show ? '' : 'none';
    colbody[3].style.display = do_show ? '' : 'none';
  } else {
    $("filter").style.display = 'none';
  }
}

function check() {
  "use strict";
  // Check that a PTM file is given.
  if (! $("psm_1").files[0]) {
    alert("Please input at least one PTM file.");
    return false;
  }
  // Check that the Modified Peptide Column Name is valid for all PTM files.
  if (! check_column_name($("psm_column_name").value)) return false;
  // Check the context sequences.
  if (sequences != null) {
    if (!sequences.check()) return false;
  }
  if ((! $('alg_mtfx').checked) && $('harvard').checked) {
    $('harvard').checked = false;
  }
  if (sequences.source.value == "noseq" && $('db_background').checked) {
    alert("Please input context sequences or uncheck 'Get background peptides' under 'Advanced options'.");
    return false;
  }
  // Check the job details.
  if (!check_job_details()) return false;
  if (!check_int_value("motif width", "width", 1, 51, 13)) return false;
  return true;
}

function options_changed() {
  if ($("filter_enable").checked) return true;
  if (!/^\s*5\s*$/.test($("occurs").value)) return true;
  if (!/^\s*13\s*$/.test($("width").value)) return true;
  if ($("single_per_mass").checked) return true;
  if ($("remove_unknowns").checked) return true;
  if (!$("eliminate_enable").checked) return true;
  if ($("db_background").checked) return true;
  if (!/^\s*13\s*$/.test($("eliminate_width").value)) return true;
  if (!/^\s*0.000001\s*$/.test($("score_threshold").value)) return true;
  if ($("harvard").checked) return true;
  if (!/^\s*100\s*$/.test($("max_motifs").value)) return true;
  if (!/^\s*50\s*$/.test($("max_iterations").value)) return true;
  if (!/^\s*10\s*$/.test($("max_no_decrease").value)) return true;
  return false;
}

function options_reset(evt) {
  $("occurs").value = 5;
  $("width").value = 13;
  $("single_per_mass").checked = false;
  $("remove_unknowns").checked = false;
  $("eliminate_enable").checked = true;
  $("eliminate_width").value = 13;
  $("db_background").checked = false;
  $("score_threshold").value = 0.000001;
  $("harvard").checked = false;
  $("max_motifs").value = 100;
  $("max_iterations").value = 50;
  $("max_no_decrease").value = 10;
}

function fix_reset() {
  on_alg_change(1);
  updateFilterMenu(); 		// update the filter menu from the first PTM file
  // Make sure "Hidden Modifications" gets turned off on form reset.
  var i, more_opts = document.getElementsByClassName("more_opts");
  for (i=0; i<more_opts.length; i++) { toggle_class(more_opts[i], 'modified', false); }
}

function reloadpage() {
  location.reload();
}

function on_form_submit(evt) {
  if (!check()) {
    evt.preventDefault();
  }
}

function on_form_reset(evt) {
  window.setTimeout(function(evt) {
    remove_ptm(null);
    fix_reset();
    sessionStorage.clear();     // reset all storage
  }, 50);
}

function on_pageshow() {
  sequences._source_update();
}

function on_load() {
  $("sequences").addEventListener("click", on_sequences_change, false);
  $("alg_simp").addEventListener("click", on_alg_change, false);
  $("alg_mtfx").addEventListener("click", on_alg_change, false);
  $("alg_modl").addEventListener("click", on_alg_change, false);
  // add on_file_change listeners
  var num_psm_rows = $("psm_rows").rows.length;
  for (var i=1; i<=num_psm_rows; i++) {
    (function() { // Need function scope for this to work.
      var psm_id = "psm_" + i;
      $(psm_id).addEventListener('change', function(){on_file_change(psm_id);}, false);
    }());
  }
  $("psm_column_name").addEventListener('change', on_column_name_change, false);
  $("more_psms").addEventListener("click", add_ptm, false);
  $("less_psms").addEventListener("click", remove_ptm, false);
  $("momo_form").addEventListener("submit", on_form_submit, false);
  $("momo_form").addEventListener("reset", on_form_reset, false);
  $("filter_enable").addEventListener('change', updateFilterDisplay, false);
  $("filter_field").addEventListener('change', on_filter_field_change, false);
  window.addEventListener('pageshow', on_pageshow, false);
  fix_reset();
}

// add a load
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    "use strict";
    window.removeEventListener("load", load, false);
    on_load();
  }, false)
})();
