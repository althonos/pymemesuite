
var alphabet = null;
var sequences = null;
var control_sequences = null;
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
  } else if (id == "background") {
    background = controler;
  }
}

function check() {
  "use strict";
  var alphs = null;
  var min_seqs = 2;
  var num_sequences = 0;
  if (alphabet != null) alphs = alphabet.get_alphabets();
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    if (alphs == null) alphs = sequences.get_alphabets();
    if (sequences.source.value == "file") {
      num_sequences = sequences.file_dbh.sequence_count;
    } else if (sequences.source.value == "text") {
      num_sequences = sequences.text_dbh.sequence_count;
    }
  }
  // Make sure there is at least min_seqs sequences unless ANR.
  if (num_sequences < min_seqs && $("dist").value != "anr") {
    alert("Please provide at least " + min_seqs + " primary sequences or\n" +
        "select the 'Any Number of Repetitions' site distribution.\n");
    return false;
  }
  if ( ! $("classic").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
    if (alphs == null) alphs = control_sequences.get_alphabets();
  }
  if (!check_int_value("number of motifs", "nmotifs", 1, null, 3)) return false;
  if (!check_job_details()) return false;
  if (background != null && !background.check(alphs)) return false;
  if (!check_int_range("minimum motif width", "minw", 6,
        "maximum motif width", "maxw", 30, 2, 300)) return false;
  if ($("dist").value !== "oops") {
    if (!check_int_range("minimum sites", "minsites", 2,
          "maximum sites", "maxsites", 600, 2, 600)) return false;
  }
  return true;
}

function options_changed() {
  if (background != null && background.changed()) return true;
  if ($("norc").checked) return true;
  if (!/^\s*6\s*$/.test($("minw").value)) return true;
  if (!/^\s*50\s*$/.test($("maxw").value)) return true;
  if ($("minsites_enable").checked) return true;
  if ($("maxsites_enable").checked) return true;
  if ($("pal").checked) return true;
  if ($("shuffle").checked) return true;
  return false;
}

function options_reset(evt) {
  if (background != null) background.reset();
  $("norc").checked = false;
  $("minw").value = 6;
  $("maxw").value = 50;
  $("sites").style.display = 'block';
  $("minsites_enable").checked = false;
  $("minsites").value = 2;
  $("minsites").disabled = true;
  $("maxsites_enable").checked = false;
  $("maxsites").value = 600;
  $("maxsites").disabled = true;
  $("pal").checked = false;
  $("shuffle").checked = false;
}

function fix_reset() {
  $("sites").style.display = ($('dist').value == 'oops' ? 'none' : 'block');
  $('discr_sequences_area').style.display = ($("classic").checked ? 'none' : 'block');
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

function on_ch_dist() {
  $("sites").style.display = ($('dist').value == 'oops' ? 'none' : 'block');
}

function on_ch_disc_mode() {
  $('discr_sequences_area').style.display = ($("classic").checked ? 'none' : 'block');
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  sequences._source_update();
  control_sequences._source_update();
  background._source_update();
}

function on_load() {
  // add listeners for the motif discovery mode
  $("classic").addEventListener("click", on_ch_disc_mode, false);
  $("psp").addEventListener("click", on_ch_disc_mode, false);
  $("de").addEventListener("click", on_ch_disc_mode, false);
  // add listener for changing the motif distribution
  $("dist").addEventListener("change", on_ch_dist, false);
  // add listener to the form to check the fields before submit
  $("meme_form").addEventListener("submit", on_form_submit, false);
  $("meme_form").addEventListener("reset", on_form_reset, false);
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
