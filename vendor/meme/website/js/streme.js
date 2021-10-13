var alphabet = null;
var sequences = null;
var control_sequences = null;

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
  if (!check_job_details()) return false;
  if (!check_num_value("p-value threshold", "pvt", 0, 0.5, 0.05)) return false;
  if (!check_int_value("motif count", "nmotifs", 1, null, 10)) return false;
  if (!check_int_range("minimum motif width", "minw", 8,
    "maximum motif width", "maxw", 15, 2, 30)) return false;
  return true;
}

function options_changed() {
  if (!/^\s*8\s*$/.test($("minw").value)) return true;
  if (!/^\s*15\s*$/.test($("maxw").value)) return true;
  if (!$("enable_pvt").checked) return true;
  if (!/^\s*0\.05\s*$/.test($("pvt").value)) return true;
  if ($("order_enable").checked) return true;
  if (!$("align_center").checked) return true;
  return false;
}

function options_reset(evt) {
  $("dna2rna").checked = false;
  $("minw").value = 8;
  $("maxw").value = 15;
  $("enable_pvt").checked = true;
  $("pvt").value = 0.05;
  $("pvt").disabled = false;
  $("nmotifs").value = 10;
  $("nmotifs").disabled = true;
  $("order_enable").checked = false;
  $("order").value = 0;
  $("order").disabled = true;
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

function on_ch_sequences() {
  var sequence_alphs = null;
  if (sequences != null) sequence_alphs = sequences.get_alphabets();
  $('dna2rna_area').style.display = (sequence_alphs == "DNA") ? 'block' : 'none';
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  sequences._source_update();
  control_sequences._source_update();
  on_ch_discr();
}

function on_load() {
  // add listeners for the motif discovery mode
  $("discr_off").addEventListener("click", on_ch_discr, false);
  $("discr_on").addEventListener("click", on_ch_discr, false);
  $("sequences").addEventListener("sequences_checked", on_ch_sequences, false);
  // add listener to the form to check the fields before submit
  $("streme_form").addEventListener("submit", on_form_submit, false);
  $("streme_form").addEventListener("reset", on_form_reset, false);
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
