
var loci = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "loci") {
    loci = controler;
  }
}

function check() {
  if (loci != null && !loci.check()) return false;
  return true;
}

function options_changed() {
  "use strict";
  if (!/^\s*0.05\s*$/.test($("max_pvalue").value)) return true;
  if ($("closest_locus").checked) return true;
  if (! $("closest_tss").checked) return true;
  return false;
}

function options_reset() {
  "use strict";
  $("max_pvalue").value = "0.05";
  $("closest_locus").checked = false;
  $("closest_tss").checked = true;
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
  //loci._source_update();
}

function on_load() {
  // add listener to the form to check the fields before submit
  $("tgene_form").addEventListener("submit", on_form_submit, false);
  $("tgene_form").addEventListener("reset", on_form_reset, false);
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
