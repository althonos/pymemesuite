
var query_motifs = null;
var target_motifs = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "query_motifs") {
    query_motifs = controler;
    // monitor changes to the query motifs
    element.addEventListener('motifs_loaded', on_ch_motifs, false);
  } else if (id == "target_motifs") {
    target_motifs = controler;
  }
}

function check() {
  "use strict";
  var alphabets = null;
  if (query_motifs != null) {
    if (!query_motifs.check()) return false;
    alphabets = query_motifs.get_alphabets();
  }
  // don't allow "alike" alphabets and DNA/RNA are not compatible
  if (target_motifs != null && !target_motifs.check(alphabets, false, true)) return false;
  if (!check_job_details($("instant_run").checked)) return false;
  if ($("thresh_type").value == "1") {
    if (!check_num_value("E-value threshold", "thresh", 0, null, 10)) return false;
  } else { // qvalue
    if (!check_num_value("q-value threshold", "thresh", 0, 0.9999, 0.05)) return false;
  }
  return true;
}


function options_changed() {
  "use strict";
  if ($("comparison_function").value !== "pearson") return true;
  if ($("thresh_type").value !== "1") return true;
  if (!/^\s*10\s*$/.test($("thresh").value)) return true;
  if (!$("complete_scoring").checked) return true;
  if ($("no_rc").checked) return true;
  return false;
}

function options_reset() {
  "use strict";
  $("comparison_function").value = "pearson";
  $("thresh_type").value = "1";
  $("thresh").value = 10;
  $("complete_scoring").checked = true;
  $("no_rc").checked = false;
}

function on_ch_instant_run() {
  "use strict";
  var email_div, instant_run_opt;
  email_div = $('email_section');
  instant_run_opt = $('instant_run');
  if (email_div != null && instant_run_opt != null) {
    email_div.style.display = (instant_run_opt.checked ? 'none' : 'block');
  }
}

function on_ch_motifs(e) {
  "use strict";
  var controler;
  controler = e.detail.controler;
  if (controler.get_motif_count() > 1 && $('instant_run').checked) {
    $('instant_run').checked = false;
    on_ch_instant_run();
  }
}

function fix_reset() {
  // Make sure "Hidden Modifications" gets turned off on form reset.
  var i, more_opts = document.getElementsByClassName("more_opts");
  for (i=0; i<more_opts.length; i++) { toggle_class(more_opts[i], 'modified', false); }
}

function on_form_submit(evt) {
  if (!check()) {
    evt.preventDefault();
  } 
}

function on_form_reset(evt) {
  window.setTimeout(function(evt) {
    fix_reset();
  }, 50);
}

function on_pageshow() {
  query_motifs._source_update();
  target_motifs._source_update();
}

function on_load() {
  "use strict";
  // monitor changes to instant run setting
  on_ch_instant_run();
  $("instant_run").addEventListener('click', on_ch_instant_run, false);
  // add a listener to hide the "processing" message when the user returns to the page via the back button
  window.addEventListener('pageshow', on_pageshow, false);
  // add listener to the form to check the fields before submit
  $("tomtom_form").addEventListener("submit", on_form_submit, false);
  $("tomtom_form").addEventListener("reset", on_form_reset, false);
}

// anon function to avoid polluting global scope
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    "use strict";
    window.removeEventListener("load", load, false);
    on_load();
  }, false);
})();
