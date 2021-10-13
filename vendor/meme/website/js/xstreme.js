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
  var sequence_alphs = null;
  var motif_alphs = null;
  if (alphabet != null) alphs = alphabet.get_alphabets();
  if (sequences != null) {
    if (!sequences.check(alphs)) return false;
    sequence_alphs = sequences.get_alphabets();
    if (alphs == null) alphs = sequences.get_alphabets();
  }
  if ($("discr_on").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
  }
  if (motifs != null) {
    if (!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    motif_alphs = motifs.get_alphabets();
    if (alphs == null) alphs = motif_alphs;
  }
  if (!check_rna(sequence_alphs, motif_alphs)) return false;
  if (!check_job_details()) return false;
  if (!check_num_value("E-value threshold", "evt", 0, 10, 0.05)) return false;
  if (!check_int_range("minimum motif width", "minw", 8,
    "maximum motif width", "maxw", 15, 2, 30)) return false;
  if (!check_int_value("central trimming region", "ctrim", 2, null, 100)) return false;
  return (general_check(alphs) && streme_check(alphs) && meme_check(alphs));
}

function general_check(alphs) {
  "use strict";
  if (background != null && !background.check(alphs)) return false;
  return true;
}

function general_changed() {
  "use strict";
  if (background != null && background.changed()) return true;
  if (!/^\s*6\s*$/.test($("minw").value)) return true;
  if (!/^\s*15\s*$/.test($("maxw").value)) return true;
  if (!/^\s*0\.05\s*$/.test($("evt").value)) return true;
  if ($("order_enable").checked) return true;
  if ($("ctrim_enable").checked) return true;
  if (!$("align_center").checked) return true;
  return false;
}

function general_reset(evt) {
  "use strict";
  if (background != null) background.reset();
  $("dna2rna").checked = false;
  $("evt").value = 0.05;
  $("minw").value = 6;
  $("maxw").value = 15;
  $("order_enable").checked = false;
  $("order").value = 0;
  $("order").disabled = true;
  $("ctrim_enable").checked = false;
  $("ctrim").value = 100;
  $("ctrim").disabled = true;
  $("align_center").checked = true;
}

function streme_check(alphs) {
  "use strict";
  if (!check_num_value("STREME E-value threshold", "streme_ethresh", 0, null, 0.05)) return false;
  if (!check_int_value("STREME number of motifs", "streme_nmotifs", 0, null, 5)) return false;
  return true;
}

function streme_changed() {
  "use strict";
  if (! $("streme_enable_default_ethresh").checked) return true;
  return false;
}

function streme_reset(evt) {
  "use strict";
  $("streme_enable_default_ethresh").checked = true;
  $("streme_enable_ethresh").checked = false;
  $("streme_ethresh").value = 0.05;
  $("streme_ethresh").disabled = true;
  $("streme_enable_nmotifs").checked = false;
  $("streme_nmotifs").value = 5;
  $("streme_nmotifs").disabled = true;
}

function meme_check(alphs) {
  "use strict";
  if (!check_num_value("MEME E-value threshold", "meme_ethresh", 0, null, 0.05)) return false;
  if (!check_int_value("MEME number of motifs", "meme_nmotifs", 0, null, 5)) return false;
  return true;
}

function meme_changed() {
  "use strict";
  var dist = $("meme_dist");
  if (! $("meme_enable_default_ethresh").checked) return true;
  if (dist.options[dist.selectedIndex].value != "zoops") return true;
  return false;
}

function meme_reset(evt) {
  "use strict";
  $("meme_enable_default_ethresh").checked = true;
  $("meme_enable_ethresh").checked = false;
  $("meme_ethresh").value = 0.05;
  $("meme_ethresh").disabled = true;
  $("meme_enable_nmotifs").checked = false;
  $("meme_nmotifs").value = 5;
  $("meme_nmotifs").disabled = true;
  $("meme_dist").value = "zoops";
}

function sea_changed() {
  "use strict";
  var dist = $("meme_dist");
  if ($("sea_seqs").checked) return true;
  return false;
}

function sea_reset(evt) {
  "use strict";
  $("sea_seqs").checked = false;
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

function check_rna(sequence_alphs, motif_alphs) {
  "use strict";
  var dna2rna = $("dna2rna").checked && (sequence_alphs[0] == "DNA");
  // Uncheck dna2rna if the sequences aren't DNA.
  $("dna2rna").checked = dna2rna;
  // Will find RNA motifs?
  if (sequence_alphs[0] == "RNA" || dna2rna) {
    // Will the discovered motifs be RNA but the database motifs are not?
    if (motif_alphs[0] != "RNA") {
      alert("STREME and MEME will discover RNA motifs but your motif database contains '" + motif_alphs[0] + "' motifs.\n"
         + "Please select (or upload) an RNA motif database.");
      return false;
    }
  } else if (sequence_alphs[0].name == "DNA" && motif_alphs[0].name == "RNA") {
    alert("STREME and MEME will discover 'DNA' motifs but your motif database contains 'RNA' motifs.\n"
       + "Please select (or upload) a 'DNA' motif database or check the 'Convert DNA to RNA' option.");
    return false;
  } else if (sequence_alphs[0].name != motif_alphs[0].name && sequence_alphs[0].like != motif_alphs[0].like) {
    alert("STREME and MEME will discover '" + sequence_alphs[0] + "' motifs but your motif database contains '" + motif_alphs[0] + "' motifs.\n"
       + "Please select (or upload) a '" + motif_alphs[0] + "' motif database.");
    return false;
  }
  return true;
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  sequences._source_update();
  control_sequences._source_update();
  motifs._source_update();
  on_ch_discr();
}

function on_load() {
  // add listeners for the motif discovery mode
  $("discr_off").addEventListener("click", on_ch_discr, false);
  $("discr_on").addEventListener("click", on_ch_discr, false);
  $("sequences").addEventListener("sequences_checked", on_ch_sequences, false);
  // add listener to the form to check the fields before submit
  $("xstreme_form").addEventListener("submit", on_form_submit, false);
  $("xstreme_form").addEventListener("reset", on_form_reset, false);
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
