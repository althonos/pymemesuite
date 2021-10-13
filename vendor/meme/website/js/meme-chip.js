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
    if (alphs == null) alphs = sequence_alphs;
  }
  if (! $("classic").checked && control_sequences != null) {
    if (!control_sequences.check(alphs)) return false;
    if (alphs == null) alphs = control_sequences.get_alphabets();
  }
  if (motifs != null) {
    if (!motifs.check(alphs, alphabet != null && alphabet.get_alphabets() != null)) return false;
    motif_alphs = motifs.get_alphabets();
    if (alphs == null) alphs = motif_alphs;
  }
  if (!check_int_range("minimum motif width", "minw", 6,
        "maximum motif width", "maxw", 15, 3, 30)) return false;
  if (!check_rna(sequence_alphs, motif_alphs)) return false;
  if (!check_job_details()) return false;
  return (general_check(alphs) && meme_check(alphs) && 
      streme_check(alphs) && centrimo_check(alphs));
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
  return false;
}

function general_reset(evt) {
  "use strict";
  if (background != null) background.reset();
  $("dna2rna").checked = false;
  $("minw").value = 6;
  $("maxw").value = 15;
}

function meme_check(alphs) {
  "use strict";
  if (!check_int_value("MEME number of motifs", "meme_nmotifs", 0, null, 3)) return false;
  if ($("meme_dist").value !== "oops") {
    if (!check_int_range("MEME minimum sites", "meme_minsites", 2,
          "MEME maximum sites", "meme_maxsites", 600, 2, 600)) return false;
  }
  return true;
}

function meme_changed() {
  "use strict";
  var dist = $("meme_dist");
  if (dist.options[dist.selectedIndex].value != "zoops") return true;
  if (!/^\s*3\s*$/.test($("meme_nmotifs").value)) return true;
  if ($("meme_minsites_enable").checked) return true;
  if ($("meme_maxsites_enable").checked) return true;
  if ($("meme_pal").checked) return true;
  return false;
}

function meme_reset(evt) {
  "use strict";
  $("meme_dist").value = "zoops";
  $("meme_nmotifs").value = 3;
  $("meme_sites").style.display = 'block';
  $("meme_minsites_enable").checked = false;
  $("meme_minsites").value = 2;
  $("meme_minsites").disabled = true;
  $("meme_maxsites_enable").checked = false;
  $("meme_maxsites").value = 600;
  $("meme_maxsites").disabled = true;
  $("meme_pal").checked = false;
}

function streme_check(alphs) {
  "use strict";
  if (!check_num_value("STREME p-value threshold", "streme_pthresh", 0, 1, 0.05)) return false;
  if (!check_int_value("STREME motif count", "streme_nmotifs", 0, null, 5)) return false;
  return true;
}

function streme_changed() {
  "use strict";
  if (!$("streme_enable_pthresh").checked) return true;
  if (!/^\s*0\.05\s*$/.test($("streme_pthresh").value)) return true;
  return false;
}

function streme_reset(evt) {
  "use strict";
  $("streme_enable_pthresh").checked = true;
  $("streme_pthresh").value = 0.05;
  $("streme_pthresh").disabled = false;
  $("streme_enable_nmotifs").checked = false;
  $("streme_nmotifs").value = 5;
  $("streme_nmotifs").disabled = true;
}

function centrimo_check(alphs) {
  "use strict";
  if (!check_num_value("CentriMo minimum site score", "centrimo_score", null, null, 5)) return false;
  if (!check_int_value("CentriMo maximum region width", "centrimo_maxreg", 1, null, 200)) return false;
  if (!check_num_value("CentriMo E-value threshold", "centrimo_ethresh", 0, null, 10)) return false;
  return true;
}

function centrimo_changed() {
  "use strict";
  if (!/^\s*5\s*$/.test($("centrimo_score").value)) return true;
  if ($("centrimo_maxreg_enable").checked) return true;
  if (!/^\s*10\s*$/.test($("centrimo_ethresh").value)) return true;
  if ($("centrimo_local").checked) return true;
  if (!$("centrimo_store_ids").checked) return true;
  return false;
}

function centrimo_reset(evt) {
  "use strict";
  $("centrimo_score").value = 5;
  $("centrimo_maxreg_enable").checked = false;
  $("centrimo_maxreg").value = 200;
  $("centrimo_maxreg").disabled = true;
  $("centrimo_ethresh").value = 10;
  $("centrimo_local").checked = false;
  $("centrimo_store_ids").checked = true;
}

function fix_reset() {
  "use strict";
  on_ch_dist();
  $('discr_sequences_area').style.display = ($("classic").checked ? 'none' : 'block');
  // Make sure "Hidden Modifications" gets turned off on form reset.
  var i, more_opts = document.getElementsByClassName("more_opts");
  for (i=0; i<more_opts.length; i++) { toggle_class(more_opts[i], 'modified', false); }
}

function on_form_submit(evt) {
  "use strict";
  if (!check()) {
    evt.preventDefault();
  }
}

function on_form_reset(evt) {
  "use strict";
  window.setTimeout(function(evt) {
    "use strict";
    fix_reset();
  }, 50);
}

function on_ch_dist() {
  "use strict";
  $("meme_sites").style.display = ($('meme_dist').value == 'oops' ? 'none' : 'block');
}

function on_ch_disc_mode() {
  $('discr_sequences_area').style.display = ($("classic").checked ? 'none' : 'block');
}

function on_ch_sequences() {
  var sequence_alphs = null;
  if (sequences != null) sequence_alphs = sequences.get_alphabets();
  $('dna2rna_area').style.display = (sequence_alphs == "DNA") ? 'block' : 'none';
}

function check_rna(sequence_alphs, motif_alphs) {
  "use strict";
  var dna2rna = $("dna2rna").checked && (sequence_alphs[0] == "DNA");
  var centrimo_local = $("centrimo_local").checked;
  // Uncheck dna2rna if the sequences aren't DNA.
  $("dna2rna").checked = dna2rna;
  // Will find RNA motifs?
  if (sequence_alphs[0] == "RNA" || dna2rna) {
    // Will the discovered motifs be RNA but the database motifs are not?
    if (motif_alphs[0] != "RNA") {
      alert("MEME and STREME will discover RNA motifs but your motif database contains '" + motif_alphs[0] + "' motifs.\n"
	 + "Please select (or upload) an RNA motif database.");
      return false;
    }
    // Turn on localized search?
    if (!centrimo_local) {
      if (confirm("CentriMo's localized search works well with this option. Enable it now?")) {
	$("centrimo_local").checked = true;
	toggle_class($('centrimo_opts'), 'expanded', true);
      }
    }
  } else if (sequence_alphs[0].name == "DNA" && motif_alphs[0].name == "RNA") {
    alert("MEME and STREME will discover 'DNA' motifs but your motif database contains 'RNA' motifs.\n"
       + "Please select (or upload) a 'DNA' motif database or check the 'Convert DNA to RNA' option.");
    return false;
  } else if (sequence_alphs[0].name != motif_alphs[0].name && sequence_alphs[0].like != motif_alphs[0].like) {
    alert("MEME and STREME will discover '" + sequence_alphs[0] + "' motifs but your motif database contains '" + motif_alphs[0] + "' motifs.\n"
       + "Please select (or upload) a '" + motif_alphs[0] + "' motif database.");
    return false;
  }
  return true;
}

function on_pageshow() {
  alphabet._radio_update(alphabet);
  motifs._source_update();
  sequences._source_update();
  control_sequences._source_update();
}

function on_load() {
  "use strict";
  // add listeners for the motif discovery and enrichment mode
  $("classic").addEventListener("click", on_ch_disc_mode, false);
  $("de").addEventListener("click", on_ch_disc_mode, false);
  $("psp").addEventListener("click", on_ch_disc_mode, false);
  $("sequences").addEventListener("sequences_checked", on_ch_sequences, false);
  // this is needed or refresh closes the control sequences area
  $('discr_sequences_area').style.display = ($("classic").checked ? 'none' : 'block');
  // setup meme dist
  on_ch_dist();
  $("meme_dist").addEventListener("change", on_ch_dist, false);
  // add listener to the form to check the fields before submit
  $("memechip_form").addEventListener("submit", on_form_submit, false);
  $("memechip_form").addEventListener("reset", on_form_reset, false);
  window.addEventListener('pageshow', on_pageshow, false);
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
