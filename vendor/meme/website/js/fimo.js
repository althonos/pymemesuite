
var motifs = null;
var sequences = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "motifs") {
    motifs = controler;
    element.addEventListener("motifs_loaded", function() {
      if (sequences != null) {
        sequences.set_expected_alphabet(motifs.get_alphabet());
      }
    }, false);
  } else if (id == "sequences") {
    sequences = controler;
    if (motifs != null) {
      sequences.set_expected_alphabet(motifs.get_alphabet());
    }
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (motifs != null) {
    alphs = motifs.get_alphabets();
  }
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
  }
  if (!check_job_details()) return false;
  return true;
}

function options_changed() {
  var output_pv = +($("output_pv").value);
  if (typeof output_pv !== "number" || isNaN(output_pv) || output_pv != 1e-4) return true;
  if ($("norc").checked) return true;
  return false;
}

function options_reset(evt) {
  $("output_pv").value = "1e-4";
  $("norc").checked = false;
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
  motifs._source_update();
  sequences._source_update();
}

function on_load() {
  // add listener to the form to check the fields before submit
  $("fimo_form").addEventListener("submit", on_form_submit, false);
  $("fimo_form").addEventListener("reset", on_form_reset, false);
  window.addEventListener('pageshow', on_pageshow, false);
}

// add a load
(function() {
  "use strict";
  window.addEventListener("load", function load(evt) {
    "use strict";
    window.removeEventListener("load", load, false);
    on_load();
  }, false);
})();

