
var sequences = null;
var primary_motif = null;
var secondary_motifs = null;
var background = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "sequences") {
    sequences = controler;
  } else if (id == "primary") {
    primary_motif = controler;
  } else if (id == "secondaries") {
    secondary_motifs = controler;
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if (primary_motif != null) {
    if (!primary_motif.check(alphs)) return false;
    if (alphs == null) alphs = primary_motif.get_alphabets();
  }
  if (secondary_motifs != null) {
    if (!secondary_motifs.check(alphs, true)) return false;
  }
  if (!check_job_details()) return false;
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function options_changed() {
  if (background != null && background.changed()) return true;
  if ($('dumpseqs').checked) return true;
  return false;
}

function options_reset() {
  if (background != null) background.reset();
  $('dumpseqs').checked = false;
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
  sequences._source_update();
  primary_motif._source_update();
  secondary_motifs._source_update();
  background._source_update();
}

function on_load() {
  // add listener to the form to check the fields before submit
  $("spamo_form").addEventListener("submit", on_form_submit, false);
  $("spamo_form").addEventListener("reset", on_form_reset, false);
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
