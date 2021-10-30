
var motifs = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "motifs") {
    motifs = controler;
  }
}

function check() {
  if (motifs != null && !motifs.check()) return false;
  if (!check_job_details()) return false;
  if (!check_num_value("significance threshold", "threshold", 0, 0.5, 0.05)) return false;
  if (!check_int_value("score shuffling rounds", "shuffle_rounds", 1, null, 1000)) return false;
  return true;
}

function options_changed() {
  "use strict";
  if (!/^\s*0\.05\s*$/.test($("threshold").value)) return true;
  if (!/^\s*1000\s*$/.test($("shuffle_rounds").value)) return true;
  if (!$("multi_genome").checked) return true;
  return false;
}

function options_reset() {
  "use strict";
  $("threshold").value = 0.05;
  $("shuffle_rounds").value = 1000;
  $("multi_genome").checked = true;
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
}

function on_load() {
  // add listener to the form to check the fields before submit
  $("gomo_form").addEventListener("submit", on_form_submit, false);
  $("gomo_form").addEventListener("reset", on_form_reset, false);
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
