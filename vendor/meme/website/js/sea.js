var alphabet = null;
var sequences = null;
var control_sequences = null;
var motifs = null;
var background = null;

function register_component(id, element, controler) {
  "use strict";
  if (id == "alphabet") {
    alphabet = controler;
    element.addEventListener("alphabet_changed", function (e) {
      if (sequences != null) sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
      if (control_sequences != null) control_sequences.set_custom_alphabet(e.detail.has_custom, e.detail.alphabet);
    }, false);
  } else if (id == "sequences") {
    sequences = controler;
    if (alphabet != null) {
      sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "control_sequences") {
    control_sequences = controler;
    if (alphabet != null) {
      control_sequences.set_custom_alphabet(alphabet.has_custom_alphabet(), alphabet.get_custom_alphabet());
    }
  } else if (id == "motifs") {
    motifs = controler;
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  if (alphabet != null) alphs = alphabet.get_alphabets();
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if ($("discr_on").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
  }
  if (motifs != null) {
    if(!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    if (alphs == null) alphs = motifs.get_alphabets();
  }
  if (!check_job_details()) return false;
  if (!check_num_value("E-value threshold", "evt", 0, 10, 1e300)) return false;
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function options_changed() {
  if (!/^\s*10\s*$/.test($("evt").value)) return true;
  if ($("order_enable").checked) return true;
  if (background != null && background.changed()) return true;
  if (!$("align_center").checked) return true;
  return false;
}

function options_reset(evt) {
  $("evt").value = 10;
  $("order_enable").checked = false;
  $("order").value = 0;
  $("order").disabled = true;
  if (background != null) background.reset();
  $("align_center").checked = true;
}

function fix_reset() {
  $('discr_sequences_area').style.display = ($('discr_on').checked ? 'block' : 'none');
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

function on_ch_discr() {
  $('discr_sequences_area').style.display = ($('discr_on').checked ? 'block' : 'none');
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  sequences._source_update();
  control_sequences._source_update();
  motifs._source_update();
  background._source_update();
  on_ch_discr();
}

function on_load() {
  // add listeners for the motif enrichment mode
  $("discr_off").addEventListener("click", on_ch_discr, false);
  $("discr_on").addEventListener("click", on_ch_discr, false);
  // add listener to the form to check the fields before submit
  $("sea_form").addEventListener("submit", on_form_submit, false);
  $("sea_form").addEventListener("reset", on_form_reset, false);
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
  }, false);
})();
