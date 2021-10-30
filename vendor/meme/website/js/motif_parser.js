var TagInfo = function (name, skip, input, script) {
  "use strict";
  this.name = name;
  this.skip = skip;
  this.input = input;
  this.script = script;
  this.count = 0;
  this.end_tag = null;
};

TagInfo.prototype.get_end_tag = function() {
  if (this.end_tag == null) {
    this.end_tag = new BMSearchArray("</" + this.name, true);
  }
  return this.end_tag;
};

// Custom error
var MotifError = function(message, reasons) {
  "use strict";
  this.message = message;
  this.stack = Error().stack;
  this.reasons = reasons;
};
MotifError.prototype = Object.create(Error.prototype);
MotifError.prototype.name = "MotifError";
MotifError.prototype.constructor = MotifError;


var MotifUtils = {};

MotifUtils.DEFAULT_SITECOUNT = 20;

MotifUtils.DEFAULT_PSEUDOCOUNT = 0.01;

MotifUtils.uniform_freqs = function (length) {
  "use strict";
  var i, freq, freqs;
  freq = 1.0 / length;
  freqs = [];
  for (i = 0; i < length; i++) {
    freqs[i] = freq;
  }
  return freqs;
};

MotifUtils.scores_to_freqs = function (scores, site_count, background, pseudo_count) {
  "use strict";
  var total_count, col, row, freq, bg_freq, freqs, freqs_row;
  if (site_count == null) site_count = MotifUtils.DEFAULT_SITECOUNT;
  if (background == null) background = MotifUtils.uniform_freqs(scores[0].length);
  if (pseudo_count == null) pseudo_count = MotifUtils.DEFAULT_PSEUDOCOUNT;
  total_count = site_count + pseudo_count;
  freqs = [];
  for (row = 0; row < scores.length; row++) {
    freqs_row = [];
    for (col = 0; col < background.length; col++) {
      bg_freq = background[col];
      // convert to a probability
      freq = Math.pow(2.0, scores[row][col] / 100.0) * bg_freq;
      // remove the pseudo count
      freq = ((freq * total_count) - (bg_freq * pseudo_count)) / site_count;
      // correct for rounding errors
      if (freq < 0) freq = 0;
      else if (freq > 1) freq = 1;
      freqs_row.push(freq);
    }
    freqs.push(freqs_row);
  }
  return freqs;
};

MotifUtils.freqs_to_scores = function (freqs, site_count, background, pseudo_count) {
  "use strict";
  var total_count, col, row, freq, score, scores, scores_row, bg_freq;
  if (freqs.length == 0) return [];
  if (site_count == null) site_count = MotifUtils.DEFAULT_SITECOUNT;
  if (background == null) background = MotifUtils.uniform_freqs(freqs[0].length);
  if (pseudo_count == null) pseudo_count = MotifUtils.DEFAULT_PSEUDOCOUNT;
  total_count = site_count + pseudo_count;
  scores = [];
  for (row = 0; row < freqs.length; row++) {
    scores_row = [];
    for (col = 0; col < background.length; col++) {
      bg_freq = background[col];
      freq = freqs[row][col];
      // apply a pseudo count
      freq = ((pseudo_count * bg_freq) + (freq * site_count)) / total_count;
      // if the background is correct this shouldn't happen
      if (freq <= 0) freq = 0.0000005;
      // convert to a score
      score = (Math.log(freq / bg_freq) / Math.LN2) * 100;
      scores_row.push(score);
    }
    scores.push(scores_row);
  }
  return scores;
};

/*
 * Static function
 * Takes a mapping from code to a list of alphabet indexes as well as the length
 * of the alphabet.
 * Generates a list of pairs containg the a code and the set of frequencies
 * that it implies.
 */
MotifUtils._code_freqs = function(alphabet) {
  var sym, cols, freqs, freq, i, j, indexes;
  cols = [];
  for (i = 0; i < alphabet.get_size_full(); i++) {
    freqs = [];
    for (j = 0; j < alphabet.get_size_core(); j++) {
      freqs[j] = 0;
    }
    sym = alphabet.get_symbol(i);
    indexes = alphabet.get_indexes(sym);
    freq = 1.0 / indexes.length;
    for (j = 0; j < indexes.length; j++) {
      freqs[indexes[j]] += freq;
    }
    cols.push({"code": sym, "freqs": freqs});
  }
  return cols;
};

/*
 * Static Function
 * Determine the best matching code to a set of alphabet frequencies by
 * euclidean distance.
 */
MotifUtils.find_code = function(column, search_space) {
  "use strict";
  var s, i, column2, diff, sum, best, code;
  code = search_space[0].code;
  best = Number.MAX_VALUE;
  // calculate the sum of squares difference
  // to each standard code and pick the best
  for (s = 0; s < search_space.length; s++) {
    column2 = search_space[s].freqs;
    sum = 0;
    for (i = 0; i < column.length; i++) {
      diff = column[i] - column2[i];
      sum += diff*diff;
    }
    if (sum < best) {
      code = search_space[s].code;
      best = sum;
    }
  }
  return code;
};

/*
 * Static Function
 * Convert a pwm into the closest codes.
 */
MotifUtils.find_codes = function(matrix, search_space) {
  var i, codes;
  codes = "";
  for (i = 0; i < matrix.length; i++) {
    codes += MotifUtils.find_code(matrix[i], search_space);
  }
  return codes;
};

MotifUtils.as_pspm = function(motif) {
  "use strict";
  return new Pspm(motif.pwm, motif.id, 0, 0, motif.nsites, motif.evalue, null, motif.alt);
};


var HtmlMotifParser = function (handler) {
  "use strict";
  var tag_names, i;
  this.give_up = false;
  this.stop = false;
  this.handler = handler;
  this.error_list = [];
  this.rating = 1;
  this.process = this._html_ready; // the current state
  this.tag_max = 0; // set below
  this.buf = ""; // temporary storage for bits of text
  this.leading_slash = false; // current tag has a leading slash
  this.trailing_slash = false; // current tag has a trailing slash
  this.partial_match = 0; // current search has a partial match in the previous chunk
  this.partial_match2 = 0; // need a second partial match for finding @JSON_VAR
  this.store = 0; // should a attribute value be stored, 1 = type, 2 = name, 3 = value
  this.gap = false; // is there a gap between the equals and a attribute value
  this.decoder = new Utf8Decoder();
  this.is_end = this._html_is_no_quote; // method to check if a char ends a attribute value when not escaped
  this.input_hidden = false;
  this.input_name = "";
  this.input_value = "";
  this.info = {"version": null, "alphabet": null, "norc": null, "background": null};
  this.motif = {};
  this.tags = {}; // lookup for html tags, set below.
  this.tag = null; // current tag
  tag_names = [ //TODO add HTML5 tag names
    "a", "abbr", "acronym", "address", "applet", "area", "b", "base",
    "basefont", "bdo", "big", "blockquote", "body", "br", "button", "caption",
    "center", "cite", "code", "col", "colgroup", "dd", "del", "dfn", "dir",
    "div", "dl", "dt", "em", "fieldset", "font", "form", "frame", "frameset",
    "h1", "h2", "h3", "h4", "h5", "h6", "head", "hr", "html", "i", "iframe",
    "img", "ins", "isindex", "kbd", "label", "legend", "li", "link", "map",
    "menu", "meta", "noframes", "noscript", "object", "ol", "optgroup",
    "option", "p", "param", "pre", "q", "s", "samp", "select", "small",
    "span", "strike", "strong", "sub", "sup", "table", "tbody", "td",
    "textarea", "tfoot", "th", "thead", "title", "tr", "tt", "u", "ul",
    "var", "xmp"];
  for (i = 0; i < tag_names.length; i++) {
    this.tags[tag_names[i]] = new TagInfo(tag_names[i], false, false, false);
    if (this.tag_max < tag_names[i].length) this.tag_max = tag_names[i].length;
  }
  this.tags["script"] = new TagInfo("script", true, false, true);
  this.tags["style"] = new TagInfo("style", true, false, false);
  this.tags["input"] = new TagInfo("input", false, true, false);
  this.tag_max = Math.max(this.tag_max, 6) + 2; // one extra for the '/' and one after to confirm the end
  this.entities = { // list of html entities and their unicode equalivent
    "nbsp": 160, "iexcl": 161, "cent": 162, "pound": 163, "curren": 164,
    "yen": 165, "brvbar": 166, "sect": 167, "uml": 168, "copy": 169, 
    "ordf": 170, "laquo": 171, "not": 172, "shy": 173, "reg": 174, "macr": 175,
    "deg": 176, "plusmn": 177, "sup2": 178, "sup3": 179, "acute": 180,
    "micro": 181, "para": 182, "middot": 183, "cedil": 184, "sup1": 185,
    "ordm": 186, "raquo": 187, "frac14": 188, "frac12": 189, "frac34": 190,
    "iquest": 191, "Agrave": 192, "Aacute": 193, "Acirc": 194, "Atilde": 195,
    "Auml": 196, "Aring": 197, "AElig": 198, "Ccedil": 199, "Egrave": 200, 
    "Eacute": 201, "Ecirc": 202, "Euml": 203, "Igrave": 204, "Iacute": 205,
    "Icirc": 206, "Iuml": 207, "ETH": 208, "Ntilde": 209, "Ograve": 210,
    "Oacute": 211, "Ocirc": 212, "Otilde": 213, "Ouml": 214, "times": 215,
    "Oslash": 216, "Ugrave": 217, "Uacute": 218, "Ucirc": 219, "Uuml": 220,
    "Yacute": 221, "THORN": 222, "szlig": 223, "agrave": 224, "aacute": 225,
    "acirc": 226, "atilde": 227, "auml": 228, "aring": 229, "aelig": 230,
    "ccedil": 231, "egrave": 232, "eacute": 233, "ecirc": 234, "euml": 235,
    "igrave": 236, "iacute": 237, "icirc": 238, "iuml": 239, "eth": 240,
    "ntilde": 241, "ograve": 242, "oacute": 243, "ocirc": 244, "otilde": 245,
    "ouml": 246, "divide": 247, "oslash": 248, "ugrave": 249, "uacute": 250,
    "ucirc": 251, "uuml": 252, "yacute": 253, "thorn": 254, "yuml": 255,
    "fnof": 402, "Alpha": 913, "Beta": 914, "Gamma": 915, "Delta": 916,
    "Epsilon": 917, "Zeta": 918, "Eta": 919, "Theta": 920, "Iota": 921,
    "Kappa": 922, "Lambda": 923, "Mu": 924, "Nu": 925, "Xi": 926,
    "Omicron": 927, "Pi": 928, "Rho": 929, "Sigma": 931, "Tau": 932,
    "Upsilon": 933, "Phi": 934, "Chi": 935, "Psi": 936, "Omega": 937,
    "alpha": 945, "beta": 946, "gamma": 947, "delta": 948, "epsilon": 949,
    "zeta": 950, "eta": 951, "theta": 952, "iota": 953, "kappa": 954,
    "lambda": 955, "mu": 956, "nu": 957, "xi": 958, "omicron": 959, "pi": 960,
    "rho": 961, "sigmaf": 962, "sigma": 963, "tau": 964, "upsilon": 965,
    "phi": 966, "chi": 967, "psi": 968, "omega": 969, "thetasym": 977,
    "upsih": 978, "piv": 982, "bull": 8226, "hellip": 8230, "prime": 8242,
    "Prime": 8243, "oline": 8254, "frasl": 8260, "weierp": 8472, "image": 8465,
    "real": 8476, "trade": 8482, "alefsym": 8501, "larr": 8592, "uarr": 8593,
    "rarr": 8594, "darr": 8595, "harr": 8596, "crarr": 8629, "lArr": 8656,
    "uArr": 8657, "rArr": 8658, "dArr": 8659, "hArr": 8660, "forall": 8704,
    "part": 8706, "exist": 8707, "empty": 8709, "nabla": 8711, "isin": 8712,
    "notin": 8713, "ni": 8715, "prod": 8719, "sum": 8721, "minus": 8722,
    "lowast": 8727, "radic": 8730, "prop": 8733, "infin": 8734, "ang": 8736,
    "and": 8743, "or": 8744, "cap": 8745, "cup": 8746, "int": 8747,
    "there4": 8756, "sim": 8764, "cong": 8773, "asymp": 8776, "ne": 8800,
    "equiv": 8801, "le": 8804, "ge": 8805, "sub": 8834, "sup": 8835,
    "nsub": 8836, "sube": 8838, "supe": 8839, "oplus": 8853, "otimes": 8855,
    "perp": 8869, "sdot": 8901, "lceil": 8968, "rceil": 8969, "lfloor": 8970,
    "rfloor": 8971, "lang": 9001, "rang": 9002, "loz": 9674, "spades": 9824,
    "clubs": 9827, "hearts": 9829, "diams": 9830, "quot": 34, "amp": 38,
    "lt": 60, "gt": 62, "OElig": 338, "oelig": 339, "Scaron": 352,
    "scaron": 353, "Yuml": 376, "circ": 710, "tilde": 732, "ensp": 8194,
    "emsp": 8195, "thinsp": 8201, "zwnj": 8204, "zwj": 8205, "lrm": 8206,
    "rlm": 8207, "ndash": 8211, "mdash": 8212, "lsquo": 8216, "rsquo": 8217,
    "sbquo": 8218, "ldquo": 8220, "rdquo": 8221, "bdquo": 8222, "dagger": 8224,
    "Dagger": 8225, "permil": 8240, "lsaquo": 8249, "rsaquo": 8250, "euro": 8364
  };
};

HtmlMotifParser.prototype.get_errors = function() {
  return this.error_list;
};

HtmlMotifParser.prototype.parser_rating = function() {
  return this.rating;
};


HtmlMotifParser.prototype._html_ready = function (chunk) {
  "use strict";
  var i;
  for (i = 0; i < chunk.length; i++) {
    if (chunk[i] == '<'.charCodeAt(0)) {
      this.process = this._html_tagname;
      return i+1;
    }
  }
};

HtmlMotifParser.prototype._html_tagname = function (chunk) {
  "use strict";
  var i, existing;
  existing = this.buf.length;
  for (i = 0; i < chunk.length; i++) {
    // whitespace = 9 to 13 or 32, '<' = 60, '>' = 62
    if ((chunk[i] >= 9 && chunk[i] <= 13) || chunk[i] == 32 || 
        chunk[i] == 60 || chunk[i] == 62) break;
    // Note that we also only allow ASCII
    if (chunk[i] > 127) {
      // this can't be a tag so go back to scanning for tags
      this.buf = "";
      this.process = this._html_ready;
      return 0;
    }
  }
  // store up to the max tag name length
  if (existing < this.tag_max) {
    this.buf += String.fromCharCode.apply(null, 
        chunk.subarray(0, Math.min(i, this.tag_max - existing)));
  }
  // when more than 3 chars, check if it's a comment
  if (existing < 3 && this.buf.length >= 3 && this.buf.lastIndexOf("!--", 0) === 0) {
    // found a comment
    this.process = this._html_comment;
    this.buf = "";
    return 3 - existing;
  }
  // check if we found the end of the name
  if (i == chunk.length) return chunk.length;
  // check for leading slash
  if (this.buf.charAt(0) == '/') {
    this.leading_slash = true;
    this.buf = this.buf.substring(1);
  } else {
    this.leading_slash = false;
  }
  // check for a trailing slash
  if (this.buf.charAt(this.buf.length -1) == '/') {
    this.trailing_slash = true;
    this.buf = this.buf.substring(0, this.buf.length - 1);
  } else {
    this.trailing_slash = false;
  }
  // find the tag (if it exists)
  if ((this.tag = this.tags[this.buf.toLowerCase()]) == null) {
    this.process = this._html_ready;
  } else {
    if (this.rating < 2 && this.tag.name === "html") this.rating = 2;
    this.process = this._html_intag;
    this.tag.count++;
  }
  // cleanup
  this.buf = "";
  return i;
};

HtmlMotifParser.prototype._html_intag = function (chunk) {
  "use strict";
  var i;
  // skip whitespace
  for (i = 0; i < chunk.length; i++) {
    if (!(chunk[i] == 32 || (chunk[i] >= 9 && chunk[i] <= 13))) break;
  }
  if (i >= chunk.length) return chunk.length; // all space
  if (chunk[i] == 47) { // check for /
    this.trailing_slash = true;
  } else if (chunk[i] == 60 || chunk[i] == 62) { // 60 = '<', 62 = '>'
    if (this.tag.input && this.input_hidden && this.input_name != "" && this.input_value != "") {
      this._html_hidden_input(this.input_name, this.input_value);
    }
    this.input_hidden = false;
    this.input_name = "";
    this.input_value = "";
    //TODO complete tag
    if (this.tag.skip && !(this.leading_slash || this.trailing_slash)) {
      this.partial_match = 0;
      this.process = this._html_skip;
    } else {
      this.tag = null;
      this.process = this._html_ready;
    }
    if (chunk[i] == 60) return i; // 60 = '<'
  } else {
    this.trailing_slash = false;
    this.process = this._html_attr_name;
    return i;
  }
  return i+1;
};

HtmlMotifParser.prototype._html_attr_name = function (chunk) {
  "use strict";
  var i, existing;
  existing = this.buf.length;
  for (i = 0; i < chunk.length; i++) {
    // '/' = 47, '<' = 60, '=' = 61, '>' = 62
    if ((chunk[i] >= 9 && chunk[i] <= 13) || chunk[i] == 32 || // whitespace
        chunk[i] == 47 || (chunk[i] >= 60 && chunk[i] <= 62)) {
      break;
    }
    if (chunk[i] > 127) {
      // this can't be an attribute (at least not in any html I'd write)
      // so go back to scanning for tags
      this.process = this._html_ready;
      this.buf = "";
      return 0;
    }
  }
  if (this.tag.input) {
    // store just enough of the attribute to tell if it's "type", "name" or "value"
    if (existing < 6) {
      this.buf += String.fromCharCode.apply(null, 
          chunk.subarray(0, Math.min(i, this.tag_max - existing)));
    }
  }
  // read more to find the end of the attribute
  if (i >= chunk.length) return chunk.length;
  // found the end, now test if it is of interest
  if (this.tag.input) {
    this.store = (this.buf === "type" ? 1 : this.buf === "name" ? 2 : this.buf === "value" ? 3 : 0);
  } else {
    this.store = 0;
  }
  // clean up
  this.buf = "";
  // now search for the equals between the attribute and value
  this.process = this._html_attr_equals;
  return i;
};

HtmlMotifParser.prototype._html_attr_equals = function (chunk) {
  "use strict";
  var i;
  // skip whitespace
  for (i = 0; i < chunk.length; i++) {
    if (!(chunk[i] == 32 || (chunk[i] >= 9 && chunk[i] <= 13))) break;
  }
  if (i >= chunk.length) return chunk.length; // all space
  // now check to see if we got an equals
  if (chunk[i] == 61) {
    this.process = this._html_attr_value;
    this.gap = false;
    return i + 1;
  } else {
    //TODO set attribute value to empty
    this.process = this._html_intag;
    return i;
  }
};

HtmlMotifParser.prototype._html_attr_value = function (chunk) {
  "use strict";
  var i;
  // skip whitespace
  for (i = 0; i < chunk.length; i++) {
    if (!(chunk[i] == 32 || (chunk[i] >= 9 && chunk[i] <= 13))) break;
    this.gap = true;
  }
  if (i >= chunk.length) return chunk.length; // all space
  // get ready to read value
  this.escaped = false;
  this.decoder.reset();
  this.buf = "";
  // check type of quote
  switch (chunk[i]) {
    case 62: // '>'
      //TODO set attribute value to empty
      this.process = this._html_intag;
      return i;
    case 34: // '"'
      this.process = this._html_in_value;
      this.is_end = this._html_is_double_quote;
      return i + 1;
    case 39: // '\''
      this.process = this._html_in_value;
      this.is_end = this._html_is_single_quote;
      return i + 1;
    default:
      if (this.gap) {
        //TODO set attribute value to empty
        this.process = this._html_intag;
      } else {
        this.process = this._html_in_value;
        this.is_end = this._html_is_no_quote;
      }
      return i;
  }
};

HtmlMotifParser.prototype._html_in_value = function (chunk) {
  "use strict";
  var code, letter;
  this.decoder.set_source(chunk);
  while ((code = this.decoder.next()) != null) {
    if (this.is_end(code)) break;
    if (this.store != 0) this.buf += String.fromCharCode(code);
  }
  if (code == null) return chunk.length;
  if (this.store != 0) {
    // search for the escaped html entities
    this.buf.replace(/&(?:#([1-9][0-9]*)|([a-zA-Z][a-zA-Z0-9]*));/g, 
      function(m) {
        "use strict";
        var entity_code = null;
        if (typeof m[1] === "string" && m[1].length > 0) {
          entity_code = parseInt(m[1], 10);
        } else {
          entity_code = this.entities[m[2]];
        }
        if (typeof entity_code === "number") {
          return String.fromCharCode(entity_code);
        } else {
          return m[0];
        }
      }
    );
    if (this.store == 1) {
      this.input_hidden = (this.buf === "hidden");
    } else if (this.store == 2) {
      this.input_name = this.buf;
    } else if (this.store == 3) {
      this.input_value = this.buf;
    }
  }
  this.buf = "";
  this.process = this._html_intag;
  return this.decoder.position();
};

HtmlMotifParser.prototype._html_is_single_quote = function (code) {
  "use strict";
  return (code === 39);
};

HtmlMotifParser.prototype._html_is_double_quote = function (code) {
  "use strict";
  return (code === 34);
};

HtmlMotifParser.prototype._html_is_no_quote = function (code) {
  "use strict";
  return ((code >= 9 && code <= 13) || code === 32 || code === 60 || code === 62);
};

HtmlMotifParser.prototype._html_comment = function (chunk) {
  "use strict";
  var search, match;
  if (typeof this._html_comment.search === "undefined") {
    this._html_comment.search = new BMSearchArray("-->");
  }
  search = this._html_comment.search;
  match = search.indexIn(chunk, -this.partial_match);
  this.partial_match = 0;
  if (match != null) {
    if (match.complete) {
      this.process = this._html_ready;
      return match.index + search.length();
    } else {
      this.partial_match = chunk.length - match.index;
    }
  }
  return chunk.length;
};

HtmlMotifParser.prototype._html_jsonvar = function (chunk, end) {
  "use strict";
  var search, match;
  if (typeof this._html_jsonvar.json_search === "undefined") {
    this._html_jsonvar.json_search = new BMSearchArray("@JSON_VAR");
  }
  search = this._html_jsonvar.json_search;
  match =  search.indexIn(chunk, -this.partial_match2);
  this.partial_match2 = 0;
  if (match != null) {
    if (match.complete) {
      this.process = this._html_jsonvar_match;
      return match.index + search.length();
    } else if (!end) {
      this.partial_match2 = chunk.length - match.index;
    }
  }
  return null;
};

HtmlMotifParser.prototype._html_jsonvar_match = function (chunk) {
  "use strict";
  // check for whitespace
  if (chunk[0] == 32 || (chunk[0] >= 9 && chunk[0] <= 13)) {
    this.process = this._html_json_start;
  } else {
    this.process = this._html_skip;
  }
  return 0;
};

HtmlMotifParser.prototype._html_skip = function (chunk) {
  "use strict";
  var search, match, alt_exit;
  search = this.tag.get_end_tag();
  match = search.indexIn(chunk, -this.partial_match);
  this.partial_match = 0; // reset partial match for next search
  if (match != null) {
    // I can only do the search this way because I know any partial match to
    // "</script" can never also be a partial match to "@JSON_VAR" because
    // the first character is different. If I was working with two strings
    // that started the same I would have to do this in a more careful way.
    if (this.tag.script && match.index > 0) {
      alt_exit = this._html_jsonvar(chunk.subarray(0, match.index), true);
      if (alt_exit != null) return alt_exit;
    }
    if (match.complete) {
      this.process = this._html_skip_match;
      return match.index + search.length();
    } else {
      this.partial_match = chunk.length - match.index;
    }
  } else if (this.tag.script) {
    alt_exit = this._html_jsonvar(chunk, false);
    if (alt_exit != null) return alt_exit;
  }
  return chunk.length;
};

HtmlMotifParser.prototype._html_skip_match = function (chunk) {
  "use strict";
  // check for whitespace, < or >
  if (chunk[0] == 32 || (chunk[0] >= 9 && chunk[0] <= 13) || chunk[0] == 60 || chunk[0] == 62) {
    this.leading_slash = true;
    this.trailing_slash = false;
    this.process = this._html_intag;
  } else {
    this.process = this._html_skip;
  }
  return 0;
};

HtmlMotifParser.prototype._html_json_start = function (chunk) {
  "use strict";
  var i;
  // skip until '{'
  for (i = 0; i < chunk.length; i++) {
    if (chunk[i] === 123) {
      this.decoder.reset();
      this.buf = "";
      this.in_json_string = false;
      this.in_json_escape = false;
      this.json_nesting = 0;
      this.process = this._html_json;
      return i; 
    }
  }
  return chunk.length; // no '{' found
};

HtmlMotifParser.prototype._html_json = function (chunk) {
  "use strict";
  var code;
  this.decoder.set_source(chunk);
  while ((code = this.decoder.next()) != null) {
    if (this.in_json_string) {
      if (this.in_json_escape) {
        this.in_json_escape = false;
      } else {
        if (code === 92) { // '\\'
          this.in_json_escape = true;
        } else if (code === 34) { // '"'
          this.in_json_string = false;
        }
      }
    } else {
      if (code === 123) { // '{'
        this.json_nesting++;
      } else if (code === 125) { // '}'
        this.json_nesting--;
      } else if (code === 34) { // '"'
        this.in_json_string = true;
      }
    }
    this.buf += String.fromCharCode(code);
    if (this.json_nesting === 0) {
      try {
        this._html_data(JSON.parse(this.buf));
      } catch (e) {
        // ignore error
      }
      this.partial_match = 0;
      this.process = this._html_skip;
      return this.decoder.position();
    }
  }
  return chunk.length;
};



HtmlMotifParser.prototype._html_hidden_input = function (name, value) {
  "use strict";
  var parts, freqs, letters, i, j, sum, lines, index, is_pssm, keyvals, key, num, matrix;
  // other unused data fields: name, combinedbloack, nmotifs, motifblock_, BLOCKS_
  if (name === "version") {
    if ((parts = /(?:^|\s)MEME\s+version\s+(\d+(?:\.\d+){0,2})(?:\s|$)/.exec(value)) != null) {
      this.info["version"] = parts[1];
    } else {
      throw new MotifError("Bad version value.");
    }
    if (this.rating < 3) this.rating = 3;
  } else if (name === "alphabet") {
    if (this.info["version"] == null) throw new MotifError("Expected version before alphabet");
    if (value == "ACGT") {
      this.info.alphabet = AlphStd.DNA;
    } else if (value === "ACDEFGHIKLMNPQRSTVWY") {
      this.info.alphabet = AlphStd.PROTEIN;
    } else {
      throw new MotifError("Bad alphabet value.");
    }
    if (this.rating < 4) this.rating = 4;
  } else if (name === "strands") {
    if (this.info.version == null) throw new MotifError("Expected version before strands");
    if (this.info.alphabet == null) throw new MotifError("Expected alphabet before strands");
    if (this.info.alphabet.has_complement()) {
      if (value == "+ -" || value == "both") {
        this.info.norc = false;
      } else if (value === "forward") {
        this.info.norc = true;
      } else {
        throw new MotifError("Stranded motifs must have 'both', 'forward' or '+ -' for strands");
      }
    } else {
      if (value == "none") {
        this.info.norc = true;
      } else {
        throw new MotifError("Unstranded motifs must have 'none' value for strands.");
      }
    }
    if (this.rating < 5) this.rating = 5;
  } else if (name === "bgfreq") {
    if (this.info.version == null) throw new MotifError("Expected version before bgfreq");
    if (this.info.alphabet == null) throw new MotifError("Expected alphabet before bgfreq");
    if (this.info.norc == null) throw new MotifError("Expected strands before bgfreq");
    if ((parts = value.match(/\S+/g)) === null || parts.length !== (2 * this.info.alphabet.get_size_core())) {
      throw new MotifError("Incorrect number of parts in bgfreq.");
    }
    freqs = [];
    sum = 0.0;
    for (i = 0; i < this.info.alphabet.get_size_core(); i++) {
      if (i !== this.info.alphabet.get_index(parts[i*2])) {
        throw new MotifError("Expected bgfreq alphabet letter " +
            this.info.alphabet.get_symbol(i) + " but got \"" + parts[i*2] + "\".");
      }
      if (!/^\d(?:\.\d+)?(?:[eE]\d+)?$/.test(parts[i*2 + 1])) {
        throw new MotifError("Expected bgfreq probability for letter " + 
            this.info.alphabet.get_symbol(i) + " but got \"" + parts[i*2 + 1] + "\".");
      }
      freqs[i] = +parts[i*2 + 1];
      if (freqs[i] > 1.0) {
        throw new MotifError("Expected bgfreq probability for letter " + 
            this.info.alphabet.get_symbol(i) + " but got " + freqs[i] + " which is out of range.");
      }
      sum += freqs[i];
    }
    if (Math.abs(sum - 1.0) > 0.01) {
      throw new MotifError("Probabilities in bgfreq don't sum to 1.0");
    }
    this.info.background = freqs;
    if (this.rating < 6) this.rating = 6;
    if (typeof this.handler.meta == "function") {
      this.handler.meta(this.info);
    }
  } else if (name.lastIndexOf("motifname", 0) === 0) {
    if (typeof this.info.version == null) throw new MotifError("Expected version before matrix");
    if (this.info.alphabet == null) throw new MotifError("Expected alphabet before matrix");
    if (this.info.norc == null) throw new MotifError("Expected strands before matrix");
    if (!(this.info.background != null && this.info.background instanceof Array)) {
      throw new MotifError("Expected bgfreq before matrix");
    }
    index = name.substring(9);
    // if this is a new motif then report the previous one
    this._report_motif(index);
    this.motif["id"] = value;
    this.motif["alt"] = "DREME"; // only DREME motifs have this section
  } else if ((is_pssm = (name.lastIndexOf("pssm", 0) === 0)) || name.lastIndexOf("pspm", 0) === 0) {
    if (typeof this.info.version == null) throw new MotifError("Expected version before matrix");
    if (this.info.alphabet == null) throw new MotifError("Expected alphabet before matrix");
    if (this.info.norc == null) throw new MotifError("Expected strands before matrix");
    if (!(this.info.background != null && this.info.background instanceof Array)) {
      throw new MotifError("Expected bgfreq before matrix");
    }
    index = name.substring(4);
    // if this is a new motif then report the previous one
    this._report_motif(index);
    if (typeof this.motif[(is_pssm ? "psm" : "pwm")] !== "undefined") {
      throw new MotifError("Already seen " + (is_pssm ? "PSSM" : "PSPM") + " for this motif.");
    }
    // split into lines
    if ((lines = value.match(/[^\r\n]+/g)) == null) lines = [];
    // filter out lines of only space
    lines = lines.filter(function(item) {return /\S+/.test(item);});
    if (lines.length < 2) {
      throw new MotifError("Insufficient lines to define a motif.");
    }
    // check first line
    if (is_pssm) {
      if (!/^\s*log-odds\s+matrix:/.test(lines[0])) {
        throw new MotifError("First line of PSSM should begin with \"log-odds matrix:\"");
      }
    } else {
      if (!/^\s*letter-probability\s+matrix:/.test(lines[0])) {
        throw new MotifError("First line of PSPM should begin with \"letter-probability matrix:\"");
      }
    }
    // get values out of first line
    keyvals = lines[0].substring(lines[0].indexOf(":")).match(/[a-zA-Z]+\s*=\s*[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?/g);
    if (keyvals != null) {
      for (i = 0; i < keyvals.length; i++) {
        parts = /([a-zA-Z]+)\s*=\s*([\d\.eE+-]+)/.exec(keyvals[i]);
        key = parts[1];
        num = +parts[2];
        if (key == "alength") {
          if (num != this.info.alphabet.get_size_core()) {
            throw new MotifError("Matrix header value 'alength' does not match alphabet.");
          }
        } else if (key == "w") {
          if (num % 1 !== 0) {
            throw new MotifError("Matrix header value 'w' should be a whole number.");
          }
          if (num != (lines.length - 1)) {
            throw new MotifError("Matrix header value 'w' does not match the number of remaining non-empty lines.");
          }
          if (typeof this.motif.len === "number") {
            if (this.motif.len !== num) {
              throw new MotifError("Matrix header value 'w' does not match previously read values.");
            }
          } else {
            this.motif.len = num;
          }
        } else if (key == "nsites") {
          // for some reason lost in the mysts of time we allow non-integers here
          if (typeof this.motif.nsites === "number") {
            if (this.motif.nsites !== num) {
              throw new MotifError("Matrix header value 'nsites' does not match previously read values.");
            }
          } else {
            this.motif.nsites = num;
          }
        } else if (key == "E") {
          if (typeof this.motif.evalue === "number") {
            if (this.motif.evalue !== num) {
              throw new MotifError("Matrix header value 'E' does not match previously read values.");
            }
          } else {
            this.motif.evalue = num;
          }
        }
      }
    }
    // now read the grid of numbers
    matrix = [];
    for (i = 1; i < lines.length; i++) {
      parts = lines[i].match(/\S+/g);
      if (parts == null || parts.length != this.info.alphabet.get_size_core()) {
        throw new MotifError("Matrix row does not have enough entries for the alphabet.");
      }
      sum = 0;
      for (j = 0; j < parts.length; j++) {
        if (!/[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?/.test(parts[j])) {
          throw new MotifError("Matrix row contains an entry \"" + parts[j] + "\" that is not a number.");
        }
        num = +parts[j];
        if (!is_pssm) {
          if (num < 0 || num > 1) {
            throw new MotifError("Matrix row for PSPM contains an entry that is not a probability.");
          }
          sum += num;
        }
        parts[j] = num;
      }
      if (!is_pssm && Math.abs(sum - 1.0) > 0.00001) { // the output uses 6 digits of precision
        throw new MotifError("Matrix row for PSPM does not sum to 1.");
      }
      matrix[i-1] = parts;
    }
    if (typeof this.motif.len !== "number") {
      this.motif.len = matrix.length;
    }
    // store
    if (is_pssm) {
      this.motif.psm = matrix;
    } else {
      this.motif.pwm = matrix;
    }
    // now if we have both then report the motif
    if (typeof this.motif.psm !== "undefined" &&
        typeof this.motif.pwm !== "undefined") {
      this._report_motif();
    }
  } else if (name === "nmotifs") {
    console.log("# motifs: " + value);
  }
};

HtmlMotifParser.prototype._report_motif = function(new_index) {
  if (typeof new_index !== "string") new_index = null;
  if (typeof this.motif["index"] !== "string") {
    // no motif to report
    if (new_index != null) this.motif = {"index": new_index};
  } else if (this.motif["index"] !== new_index) {
    //TODO check things, generate psm when missing?
    if (typeof this.motif["id"] !== "string") {
      this.motif["id"] = this.motif["index"];
      this.motif["alt"] = "MEME";
    }
    // report motif
    if (typeof this.handler.motif == "function") {
      this.handler.motif(this.motif);
    }
    if (this.rating < 10) this.rating = 10;
    if (new_index != null) {
      this.motif = {"index": new_index};
    } else {
      this.motif = {};
    }
  }
};

HtmlMotifParser.prototype._html_data = function (data) {
  "use strict";
  var motifs, i, motif;
  var version = data.version;
  if (version == null || !/^\d+(?:\.\d+){0,2}$/.test(version)) {
    throw new MotifError("Invalid version string \"" + version + "\".");
  }
  var program = data.program.toUpperCase();
  if (program == null || (program != "MEME" && program != "DREME" && program != "STREME")) {
    throw new MotifError("Unrecognized program, expected MEME, STREME or DREME.");
  }
  this.rating = 3;
  var alphabetData = data.alphabet;
  if (alphabetData == null) {
    throw new MotifError("No alphabet defined");
  }
  var alphabet;
  var background;
  var norc;
  if (typeof alphabetData.symbols === "string") {
    if (alphabetData.symbols == "ACGT") {
      alphabet = AlphStd.DNA;
      this.rating = 4;
      if (alphabetData.strands == "both") {
        norc = false;
      } else if (alphabetData.strands == "forward") {
        norc = true;
      } else {
        throw new MotifError("Invalid alphabet strands \"" + alphabetData.strands + "\".");
      }
      this.rating = 5;
    } else if (alphabetData.symbols == "ACDEFGHIKLMNPQRSTVWY") {
      alphabet = AlphStd.PROTEIN;
      this.rating = 4;
      if (alphabetData.strands != "none") {
        throw new MotifError("Invalid alphabet strands \"" + alphabetData.strands + "\".");
      }
      norc = true;
      this.rating = 5;
    } else {
      throw new MotifError("Unrecognized alphabet string \"" + alphabetData.symbols + "\".");
    }
    if (typeof alphabetData.freqs != "object" || !(alphabetData.freqs instanceof Array)) {
      throw new MotifError("Missing background");
    }
    if (alphabetData.freqs.length != alphabet.get_size_core()) {
      throw new MotifError("background freqs is not the right length. " +
          "Expected length " + alphabet.get_size_core() + " but got length " +
          alphabetData.freqs.length + ".");
    }
    for (i = 0; i < alphabetData.freqs.length; i++) {
      if (typeof alphabetData.freqs[i] != "number") {
        throw new MotifError("Background freqs array index " + i + " is not a number");
      }
    }
    background = alphabetData.freqs;
    this.rating = 6;
  } else {
    if (typeof alphabetData.ncore != "number") {
      throw new MotifError("Missing alphabet core count");
    }
    if (typeof alphabetData.symbols != "object" || !(alphabetData.symbols instanceof Array)) {
      throw new MotifError("Missing alphabet symbols");
    }
    var symbolData;
    for (i = 0; i < alphabetData.symbols.length; i++) {
      symbolData = alphabetData.symbols[i];
      if (typeof symbolData.symbol != "string" || !/[A-Za-z0-9.*-?]/.test(symbolData.symbol)) {
        throw new MotifError("Alphabet symbol is outside expected posibilities");
      }
    }
    alphabet = new Alphabet(alphabetData);
    this.rating = 4;
    if (program == "STREME") {
      norc = true;
    } else {
      if (typeof data.options != "object" || typeof data.options.revcomp != "boolean") {
        throw new MotifError("Missing revcomp");
      }
      norc = !data.options.revcomp;
    }
    this.rating = 5;
    if (typeof data.background != "object" || 
        typeof data.background.freqs != "object" ||
        !(data.background.freqs instanceof Array)) {
      throw new MotifError("Missing background");
    }
    if (data.background.freqs.length != alphabet.get_size_core()) {
      throw new MotifError("Background freqs array is not the right length. " +
          "Expected length " + alphabet.get_size_core() + " but got length " +
          data.background.freqs.length + ".");
    }
    for (i = 0; i < data.background.freqs.length; i++) {
      if (typeof data.background.freqs[i] != "number") {
        throw new MotifError("Background freqs array index " + i + " is not a number");
      }
    }
    background = data.background.freqs;
    this.rating = 6;
  }
  if (typeof this.handler.meta == "function") {
    this.handler.meta({
      "version": version, 
      "alphabet": alphabet, 
      "norc": norc,
      "background": background
    });
  }
  motifs = data["motifs"];
  for (i = 0; i < motifs.length; i++) {
    motif = motifs[i];
    // TODO test motif
    if (typeof this.handler.motif == "function") this.handler.motif(motif);
    this.rating = 10;
  }
};

HtmlMotifParser.prototype.process_chunk = function (chunk) {
  "use strict";
  var offset, consumed;
  if (this.stop) return;
  offset = 0;
  while (offset < chunk.length) {
    try {
      consumed = this.process(chunk.subarray(offset));
    } catch (e) {
      if (e instanceof MotifError) {
        this.error_list.push({"is_error": true, "message": e.message, "reasons": e.reasons});
        this.stop = true;
        return;
      }
      throw e;
    }
    if (consumed < 0) throw new Error("Can not consume a negative amount!");
    offset += consumed;
  }
};

HtmlMotifParser.prototype.process_end = function() {
  "use strict";
  if (typeof this.motif["index"] !== "undefined") {
    this._report_motif();
  }
};

var TextMotifParser = function(handler) {
  "use strict";
  this.handler = handler;
  this.decoder = new Utf8Decoder();
  this.alph_parser = null;
  this.line = "";
  this.process = this._find_version;
  this.stop = false;
  this.sent_meta = false;
  this.is_html = false;
  this.started_with_bl_line = null;
  this.counter = 0;
  this.freq_array_ref;
  this.sum = 0;
  this.version = null;
  this.data_file = null;
  this.alphabet = null;
  this.norc = null;
  this.background_source = null;
  this.letter_freqs = null;
  this.background = null;
  this.motif_id = null;
  this.motif_alt = null;
  this.motif_len = null;
  this.motif_sites = null;
  this.motif_evalue = null;
  this.motif_url = null;
  this.motif_pspm = null;
  this.motif_pssm = null;
  this.error_list = [];
  this.rating = 1;
};

TextMotifParser.prototype.get_errors = function() {
  return this.error_list;
};

TextMotifParser.prototype.parser_rating = function() {
  return this.rating;
};

TextMotifParser.prototype._find_version = function (line) {
  "use strict";
  var match, ver;
  if (/<html>/i.test(line)) this.is_html = true;
  if ((match = /^\s*MEME\s+version\s+(\d+(?:\.\d+){0,2}).*$/.exec(line)) != null) {
    if (this.is_html) {
      ver = match[1].split(".").map(function(val){return +val;});
      while (ver.length < 3) ver.push(0);
      if (ver[0] > 4 || (ver[0] == 4 && ver[1] > 3) || 
          (ver[0] == 4 && ver[1] == 3 && ver[2] > 2)) {
        this.stop = true;
        return;
      }
    }
    this.version = match[1];
    this.process = this._pre_motif;
    this.rating = 3;
  }
  return true;
};

TextMotifParser.prototype._pre_motif = function (line) {
  "use strict";
  var match, line_part;
  if ((match = /^\s*ALPHABET\s*=\s*(\S*)\s*$/i.exec(line)) != null) {
    var alphabet;
    if (match[1] == "ACGT") {
      alphabet = AlphStd.DNA;
    } else if (match[1] == "ACGU") {
      alphabet = AlphStd.RNA;
    } else if (match[1] == "ACDEFGHIKLMNPQRSTVWY") {
      alphabet = AlphStd.PROTEIN;
    } else {
      throw new MotifError("Invalid alphabet " + match[1]);
    }
    if (this.alphabet == null) {
      this.alphabet = alphabet;
    } else if (!this.alphabet.equals(alphabet)) {
      throw new MotifError("Alphabet definition does not match previous!");
    }
    if (this.rating < 4) this.rating = 4;
  } else if (/^\s*ALPHABET(?:\s+\"|\s*$)/.test(line)) {
    this.counter = 0;
    this.alph_parser = new AlphabetParser();
    this.process = this._in_alphabet;
    return false;
  } else if ((match = /^\s*strands\s*:([\s\+\-]*)$/i.exec(line)) != null) {
    if (match[1].indexOf("+") != -1) {
      if (match[1].indexOf("-") != -1) {
        this.norc = false;
      } else {
        this.norc = true;
      }
      if (this.rating < 5) this.rating = 5;
    } else {
      throw new MotifError("Invalid strands value \"" + match[1] + "\"."); 
    }
  } else if (/^Letter frequencies in dataset:$/i.test(line)) {
    this.counter = 0;
    this.sum = 0;
    this.letter_freqs = [];
    this.freq_array_ref = this.letter_freqs;
    this.process = this._in_freqs;
  } else if ((match = /^\s*Background\s+letter\s+frequencies(\s.*)?$/i.exec(line)) != null) {
    if ((match = /^\s+\(from\s+(.*)\):.*$/i.exec(match[1])) != null) {
      this.background_source = match[1];
    }
    this.counter = 0;
    this.sum = 0;
    this.background = [];
    this.freq_array_ref = this.background;
    this.process = this._in_freqs;
  } else if ((match = /^\s*(BL\s+)?MOTIF\s*(\S+)(\s.*)?$/.exec(line)) != null) {
    this.process = this._in_motif;
    return false;
  }
  
  return true;
};

TextMotifParser.prototype._in_alphabet = function (line) {
  var done;
  try {
    if (/^\s*END ALPHABET\s*$/.test(line)) {
      done = true;
    } else if (/^\*{80}$/.test(line)) {
      if (this.counter >= 2) done = true;
    } else {
      this.alph_parser.parse_line(line);
    }
    this.counter++;
    if (done) {
      var alphData = this.alph_parser.parse_done();
      var alphabet = new Alphabet(alphData);
      if (this.alphabet == null) {
        this.alphabet = alphabet;
      } else if (!this.alphabet.equals(alphabet)) {
        throw new MotifError("Alphabet definition does not match previous!");
      }
      this.process = this._pre_motif;
      if (this.rating < 4) this.rating = 4;
    }
  } catch (e) {
    if (e instanceof AlphabetError) {
      throw new MotifError(e.message, e.reasons);
    } else {
      throw e;
    }
  }
  return true;
};

TextMotifParser.prototype._alphabet_test = function (index, letter) {
  if (letter.length != 1) throw new MotifError("Letter value should be only one character");
  if (this.alphabet == null) {
    switch (index) {
      case 0:
        return (letter == "A");
      case 1:
        return (letter == "C");
      case 2:
        if (letter == "D") {
          this.alphabet = AlphStd.PROTEIN;
          return true;
        } else if (letter == "G") {
          this.alphabet = AlphStd.DNA;
          return true;
        }
        return false;
      default:
        throw new Error("Should not still be attempting to guess by the 3th letter (index = " + index + ").");
    }
  } else {
    return (index < this.alphabet.get_size_core() && this.alphabet.get_symbol(index) == letter);
  }
};

TextMotifParser.prototype._in_freqs = function (line) {
  "use strict";
  var parts, i, letter, freq;
  // split into words separated by spaces
  parts = line.match(/\S+/g);
  if (parts == null) parts = [];
  if (parts.length % 2 != 0) {
    throw new MotifError("Expected letter and frequency pairs. Got \"" + line + "\".");
  }
  for (i = 0; i < parts.length; i += 2) {
    letter = parts[i];
    if (!this._alphabet_test(this.counter, letter)) {
      throw new MotifError("Unrecognized alphabet letter \"" + letter + "\" at position " + this.counter + ".");
    }
    if (!/^[\+]?\d*\.?\d+(?:[eE][\-\+]?\d+)?$/.test(parts[i + 1])) {
      throw new MotifError("Letter frequency is not the accepted number format \"" + parts[i + 1] + "\" at position " + this.counter + ".");
    }
    freq = +parts[i + 1];
    this.freq_array_ref[this.counter] = freq;
    this.sum += freq;
    this.counter++;
  }
  if (this.alphabet != null && this.counter >= this.alphabet.get_size_core()) {
    if (Math.abs(this.sum - 1.0) > 0.1) {
      throw new MotifError("The letter frequencies do not sum to 1.");
    }
    this.process = this._pre_motif;
  }
  if (this.rating < 6) this.rating = 6;
  return true;
};

TextMotifParser.prototype._enqueue_motif = function() {
  "use strict";
  var i, uniform;
  if (this.motif_id == null) return;
  if (!this.sent_meta) {
    // generate missing fields
    if (this.alphabet == null) throw new Error("Expected to have worked out the alphabet by now!");
    if (this.norc == null) this.norc = !this.alphabet.has_complement();
    if (this.background == null)
      this.background = MotifUtils.uniform_freqs(this.alphabet.get_size_core());
    if (typeof this.handler.meta == "function") {
      this.handler.meta({
        "version": this.version,
        "alphabet": this.alphabet,
        "norc": this.norc,
        "background": this.background
      });
    }
    this.sent_meta = true;
  }
  // check all motif fields are available
  if (this.motif_pspm == null && this.motif_pssm == null) {
    throw new MotifError("Missing motif score and probability matrix");
  }
  // set defaults when not available
  if (this.motif_sites == null) this.motif_sites = 20;
  if (this.motif_evalue == null) this.motif_evalue = 0;
  // convert 
  if (this.motif_pspm == null) {
    this.motif_pspm = MotifUtils.scores_to_freqs(this.motif_pssm, 
        this.motif_sites, this.background);
  } else if (this.motif_pssm == null) {
    this.motif_pssm = MotifUtils.freqs_to_scores(this.motif_pspm,
        this.motif_sites, this.background);
  }
  // send motif
  if (typeof this.handler.motif == "function") {
    this.handler.motif({
      "id": this.motif_id,
      "alt": this.motif_alt,
      "len": this.motif_len,
      "nsites": this.motif_sites,
      "evalue": this.motif_evalue,
      "psm": this.motif_pssm,
      "pwm": this.motif_pspm 
    });
  }
  if (this.rating < 10) this.rating = 10;

  // reset everything
  this.motif_id = null;
  this.motif_alt = null;
  this.motif_len = null;
  this.motif_sites = null;
  this.motif_evalue = null;
  this.motif_url = null;
  this.motif_pspm = null;
  this.motif_pssm = null;
};

TextMotifParser.prototype._parse_keyvals = function(text, alpha_key, width_key, sites_key, evalue_key) {
  "use strict";
  var parts, keyval_re, i, keyval, key, value;
  keyval_re = /([a-zA-Z]+)\s*=\s*([\+]?\d*\.?\d+(?:[eE][\-\+]?\d+)?)/g;
  parts = text.match(keyval_re);
  for (i = 0; i < parts.length; i++) {
    keyval_re.lastIndex = 0; // as RE is global, need to reset it
    keyval = keyval_re.exec(parts[i]);
    key = keyval[1];
    value = +keyval[2];
    if (key === alpha_key) {
      if (this.alphabet == null) {
        if (value == 4) {
          this.alphabet = AlphStd.DNA;
        } else if (value == 20) {
          this.alphabet = AlphStd.PROTEIN;
        } else {
          throw new MotifError("Expected the motif to have either 4 or 20 columns.");
        }
      } else if (this.alphabet.get_size_core() != value) {
        throw new MotifError("Alphabet length does not match previously specified alphabet.");
      }
    } else if (key === width_key) {
      if (value % 1 != 0 || value < 1) throw new MotifError("Motif width must be a positive whole number.");
      this.motif_len = value;
    } else if (key === sites_key) {
      if (value <= 0) throw new MotifError("Motif sites must be positive.");
      this.motif_sites = value;
    } else if (key === evalue_key) {
      if (value < 0) throw new MotifError("Motif evalue must be positive.");
      this.motif_evalue = value;
    }
  }
};

TextMotifParser.prototype._in_motif = function (line) {
  "use strict";
  var match, bl_line;
  if ((match = /^\s*(BL\s+)?MOTIF\s*(\S+)(\s.*)?$/.exec(line)) != null) {
    bl_line = (match[1] != null && match[1].length > 2);
    if (this.motif_id != null && match[2] == this.motif_id) {
      // check if we've seen this motif intro before as the old html format 
      // has a nasty habit of repeating the intro
      if (this.is_html) return true;
      // Additionally it's possible to have one or both of:
      // MOTIF 1 width=18 seqs=18
      // BL    MOTIF 1 width=18 seqs=18
      // So only conclude this is a new motif if we see the same type twice
      if (this.started_with_bl_line != bl_line) return true;
    }
    // enqueue existing motif
    this._enqueue_motif();
    // create a new motif
    this.started_with_bl_line = bl_line;
    this.motif_id = match[2];
    // try to find a alternative name
    if (!this.is_html) {
      if (match[3] != null && match[3].length > 1 && 
          (match = /^\s+([^\s=]+)(?:\s+(?:[^\s=].*)?)?$/.exec(match[3])) != null) {
        this.motif_alt = match[1];
      }
    } else { // only really old HTML files
      this.motif_alt = "MEME";
    }
  } else if ((match = /^\s*letter\s*-\s*probability\s+matrix\s*:(.*)$/.exec(line)) != null) {
    // assert that we've seen a motif id
    if (this.motif_id == null) throw new Error("Impossible state");
    // check this is not a repeat
    if (this.motif_pspm != null) {
      throw new MotifError("Repeated \"letter-probability matrix\" section in motif " + this.motif_id + ".");
    }
    this._parse_keyvals(match[1], "alength", "w", "nsites", "E");
    this.process = this._in_pspm;
  } else if ((match = /^\s*log\s*-\s*odds\s+matrix\s*:(.*)$/.exec(line)) != null) {
    // assert that we've seen a motif id
    if (this.motif_id == null) throw new Error("Impossible state");
    // check this is not a repeat
    if (this.motif_pssm != null) {
      throw new MotifError("Repeated \"log-odds matrix\" section in motif " + this.motif_id + ".");
    }
    this._parse_keyvals(match[1], "alength", "w", null, "E");
    this.process = this._in_pssm;
  } else if ((match = /^\s*URL\s*(\S*)\s*/.exec(line)) != null) {
    this.motif_url = match[1];
  } else if (/^\tCombined block diagrams: non-overlapping sites with p-value < [\+]?\d*\.?\d+(?:[eE][\-\+]?\d+)?$/.test(line)) {
    // end of motifs
    this._enqueue_motif();
    this.process = this._after_motifs;
  }
  return true;
};

TextMotifParser.prototype._parse_nums = function (line, is_prob) {
  "use strict";
  var parts, i, sum;
  parts = line.match(/\S+/g);
  if (parts == null || parts.length == 0) return null;
  sum = 0;
  for (i = 0; i < parts.length; i++) {
    if (!/^[\-\+]?\d*\.?\d+(?:[eE][\-\+]?\d+)?$/.test(parts[i])) {
      throw new MotifError("Value is not a number.");
    }
    parts[i] = +parts[i];
    if (is_prob) {
      if (parts[i] < 0 || parts[i] > 1) throw new MotifError("Value is not a probability.");
      sum += parts[i];
    }
  }
  if (is_prob && Math.abs(sum - 1.0) > 0.1) {
    throw new MotifError("Array of probabilitys do not sum to 1.");
  }
  // check the column count, determine alphabet if we haven't already
  if (this.alphabet != null) {
    if (this.alphabet.get_size_core() != parts.length) {
      throw new MotifError("Expected row to have " + this.alphabet.get_size_core() + " columns.");
    }
  } else {
    if (parts.length != 4 && parts.length != 20) 
      throw new MotifError("Expected the row to have either 4 or 20 columns.");
    this.alphabet = (parts.length == 4 ? AlphStd.DNA : AlphStd.PROTEIN);
  }
  return parts;
};

TextMotifParser.prototype._in_pspm = function (line) {
  "use strict";
  var row;
  if (this.motif_pspm == null) this.motif_pspm = [];
  // parse the row
  if ((row = this._parse_nums(line, true)) != null) {
    this.motif_pspm.push(row);
  } else if (this.motif_len == null && this.motif_pspm.length > 0) {
    this.motif_len = this.motif_pspm.length;
  }
  if (this.motif_len != null && this.motif_pspm.length >= this.motif_len) {
    this.process = this._in_motif;
  }
  return true;
};

TextMotifParser.prototype._in_pssm = function (line) {
  "use strict";
  var row;
  if (this.motif_pssm == null) this.motif_pssm = [];
  // parse the row
  if ((row = this._parse_nums(line, false)) != null) {
    this.motif_pssm.push(row);
  } else if (this.motif_len == null && this.motif_pssm.length > 0) {
    this.motif_len = this.motif_pssm.length;
  }
  if (this.motif_len != null && this.motif_pssm.length >= this.motif_len) {
    this.process = this._in_motif;
  }
  return true;
};

TextMotifParser.prototype._after_motifs = function (line) {
  "use strict";
  return true;
};

TextMotifParser.prototype.process_chunk = function(chunk) {
  "use strict";
  var code, letter, done, loop_check;
  if (this.stop) return;
  this.decoder.set_source(chunk);
  while ((code = this.decoder.next()) != null) {
    letter = String.fromCharCode(code);
    if (letter == '\n') {
      // return the line - first remove any leading or trailing carrage return chars
      // as I'm assuming that they're part of the newline
      if (this.line.charAt(0) == '\r') {
        this.line = this.line.slice(1);
      }
      if (this.line.charAt(this.line.length - 1) == '\r') {
        this.line = this.line.slice(0, -1);
      }
      loop_check = 0;
      do {
        if (loop_check++ > 6) throw new Error("Infinite Loop");
        try {
          done = this.process(this.line);
        } catch (e) {
          if (e instanceof MotifError) {
            this.error_list.push({"is_error": true, "message": e.message, "reasons": e.reasons});
            this.stop = true;
            return;
          }
          throw e;
        }
      } while (!done);
      this.line = "";
      continue;
    }
    // add a letter to the line
    this.line += letter;
  }

};

TextMotifParser.prototype.process_end = function() {
  "use strict";
  this._enqueue_motif();
};

var XmlMotifParser = function(handler) {
  "use strict";
  this.handler = handler;
  this.error_list = [];
  this.rating = 1;
};

XmlMotifParser.prototype.get_errors = function() {
  return this.error_list;
};

XmlMotifParser.prototype.parser_rating = function() {
  return this.rating;
};

XmlMotifParser.prototype._meme_alphabet_array = function(array_node, id_order, isPWA) {
  "use strict";
  var i, id, value;
  var id_to_value_map = {};
  for (i = 0; i < array_node.childNodes.length; i++) {
    if (array_node.childNodes[i].nodeType != 1) continue;
    id = array_node.childNodes[i].getAttribute("letter_id");
    if (id == null) throw new MotifError("MEME XML missing value alphabet_array/value/@letter_id .");
    value = array_node.childNodes[i].textContent;
    if (!/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(value)) {
      throw new MotifError("MEME XML alphabet_array contains a value that is not a number \"" + value + "\".");
    }
    value = +value;
    if (isPWA) {
      if (value < 0 || value > 1) {
        throw new MotifError("MEME XML alphabet_array contains a value that is not a probability.");
      }
    }
    id_to_value_map[id] = value;
  }
  var sum = 0;
  var array = [];
  for (i = 0; i < id_order.length; i++) {
    id = id_order[i];
    value = id_to_value_map[id];
    if (value == null) throw new MotifError("MEME XML alphabet_array missing value for id \"" + id + "\".");
    array.push(value);
    sum += value;
  }
  if (isPWA && Math.abs(sum - 1.0) > 0.01) {
    throw new MotifError("MEME alphabet_array contains probabilities that do not sum to 1.0.");
  }
  return array;
};

XmlMotifParser.prototype._meme_alphabet_matrix = function (am_node, id_order, isPWA) {
  "use strict";
  var i;
  var matrix = [];
  var aa_nodes = am_node.querySelectorAll("alphabet_array");
  for (i = 0; i < aa_nodes.length; i++) {
    matrix.push(this._meme_alphabet_array(aa_nodes[i], id_order, isPWA));
  }
  return matrix;
};

XmlMotifParser.prototype._meme_motif = function (motif_node, symbol_ids) {
  "use strict";
  // read the name
  var id = motif_node.getAttribute("name");
  if (id == null) throw new MotifError("MEME XML motif does not have a name attribute.");
  var alt = "MEME";
  // read the width
  var len = motif_node.getAttribute("width");
  if (len == null || !/^\d+$/.test(len)) throw new MotifError("MEME XML motif does not have a numeric width attribute.");
  len = +len;
  if (len < 1) throw new MotifError("MEME XML motif specifies a width less than 1.");
  // read the nsites
  var nsites = motif_node.getAttribute("sites");
  if (nsites == null || !/^\d+$/.test(nsites)) throw new MotifError("MEME XML motif does not have a numeric sites attribute.");
  nsites = +nsites;
  if (nsites < 2) throw new MotifError("MEME XML motif specifies a site count less than 2.");
  // read the evalue
  var evalue = motif_node.getAttribute("e_value");
  if (evalue == null || !/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(evalue)) throw new MotifError("MEME XML motif does not have valid e_value attribute.");
  evalue = +evalue;
  // read the score matrix
  var scores_am_node = motif_node.querySelector("scores > alphabet_matrix");
  if (scores_am_node == null) throw new MotifError("MEME XML motif missing scores matrix");
  var psm = this._meme_alphabet_matrix(scores_am_node, symbol_ids, false);
  if (psm.length != len) throw new MotifError("MEME XML motif PSM is not the expected length");
  // read the probability matrix
  var probs_am_node = motif_node.querySelector('probabilities > alphabet_matrix');
  if (probs_am_node == null) throw new MotifError("MEME XML motif missing probabilities matrix");
  var pwm = this._meme_alphabet_matrix(probs_am_node, symbol_ids, true);
  if (pwm.length != len) throw new MotifError("MEME XML motif PWM is not the expected length");
  // send the motif data
  if (typeof this.handler.motif == "function") {
    this.handler.motif({
      "id": id,
      "alt": alt,
      "len": len,
      "nsites": nsites,
      "evalue": evalue,
      "psm": psm,
      "pwm": pwm
    });
  }
  if (this.rating < 10) this.rating = 10;
};

// special version for internet explorer which doesn't use evaluate
XmlMotifParser.prototype.process_meme_doc = function (doc) {
  "use strict";
  var i, version, alphabet, strands, background;
  var meme = doc.querySelector("MEME");
  if (meme == null) throw new MotifError("MEME XML missing /MEME .");
  // get the version
  version = meme.getAttribute("version");
  if (version == null) throw new MotifError("MEME XML missing /MEME/@version .");
  if (!/^\d+(\.\d+){0,2}$/.test(version)) {
    throw new MotifError("MEME XML value /MEME/@version is not formatted correctly \"" + version + "\".");
  }
  this.rating = 2;
  // get the alphabet
  var ele_alphabet = doc.querySelector("MEME > training_set > alphabet");
  if (ele_alphabet == null) throw new MotifError("MEME XML missing /MEME/training_set/alphabet .");
  var letter_elems = doc.querySelectorAll("MEME > training_set > alphabet > letter");
  var letter_id, letter_symbol, letter_aliases, letter_name, letter_colour, letter_complement, letter_equals;
  var alphabet_id = ele_alphabet.getAttribute("id");
  if (alphabet_id != null) {
    if (alphabet_id === "nucleotide") {
      alphabet = AlphStd.DNA;
    } else if (alphabet_id === "amino-acid") {
      alphabet = AlphStd.PROTEIN;
    } else {
      throw new MotifError("MEME XML value /MEME/training_set/alphabet/@id is unrecognized \"" + alphabet_id + "\".");
    }
  } else {
    try {
      var alph_parser = new AlphabetParser();
      alph_parser.parse_header(ele_alphabet.getAttribute("name"));
      for (i = 0; i < letter_elems.length; i++) {
        letter_symbol = letter_elems[i].getAttribute("symbol");
        if (letter_symbol == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[" + i + "]/@symbol .");
        letter_aliases = letter_elems[i].getAttribute("aliases");
        letter_name = letter_elems[i].getAttribute("name");
        letter_colour = letter_elems[i].getAttribute("colour");
        letter_complement = letter_elems[i].getAttribute("complement");
        letter_equals = letter_elems[i].getAttribute("equals");
        alph_parser.parse_symbol(letter_symbol, {
          "aliases": letter_aliases, "name": letter_name, "colour": letter_colour,
          "complement": letter_complement, "comprise": letter_equals
        });
      }
      var alph_data = alph_parser.parse_done();
      alphabet = new Alphabet(alph_data);
    } catch (e) {
      if (e instanceof AlphabetError) {
        throw new MotifError(e.message, e.reasons);
      } else {
        throw e;
      }
    }
  }
  this.rating = 3;
  // make a list of the letter_ids in order of the alphabet
  var symbol_map = {};
  for (i = 0; i < letter_elems.length; i++) {
    letter_symbol = letter_elems[i].getAttribute("symbol");
    if (letter_symbol == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[" + i + "]/@symbol .");
    letter_id = letter_elems[i].getAttribute("id");
    if (letter_id == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[" + i + "]/@id .");
    symbol_map[letter_symbol] = letter_id;
  }
  var symbol_ids = [];
  for (i = 0; i < alphabet.get_size_core(); i++) {
    var id = symbol_map[alphabet.get_symbol(i)];
    if (id == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[@symbol='" + alphabet.get_symbol(i) + "']/@id .");
    symbol_ids.push(id);
  }
  // get the strands
  var strands_ele = doc.querySelector("MEME > model > strands");
  if (strands_ele == null) throw new MotifError("MEME XML missing /MEME/model/strands .");
  strands = strands_ele.textContent;
  if (strands == "both" || strands == "forward" || strands == "none") {
    if (alphabet.has_complement()) {
      if (strands == "none") throw new MotifError(
          "Strand type \"none\" is not allowed to be used for \"" + alphabet.get_alphabet_name() + "\".");
    } else {
      if (strands != "none") throw new MotifError(
          "Strand type \"" + strands + "\" is not allowed to be used for \"" +
          alphabet.get_alphabet_name() + "\".");
    }
  } else {
    throw new MotifError("MEME XML contains an unrecognized strand type \"" + strands + "\".");
  }
  this.rating = 4;
  // get the background
  var background_ele = doc.querySelector("MEME > model > background_frequencies > alphabet_array");
  if (background_ele == null) throw new MotifError("MEME XML missing /MEME/model/background_frequencies/alphabet_array .");
  background = this._meme_alphabet_array(background_ele, symbol_ids, true);
  this.rating = 5;
  // send the meta data
  if (typeof this.handler.meta == "function") {
    this.handler.meta({
      "version": version,
      "alphabet": alphabet,
      "norc": (strands == "forward" || strands == "none"),
      "background": background
    });
  }
  // now iterate over motifs
  var motifs = doc.querySelectorAll("MEME > motifs > motif");
  for (i = 0; i < motifs.length; i++) {
    this._meme_motif(motifs[i], symbol_ids);
  }
};

XmlMotifParser.prototype._dreme_probability_array = function (probs_node, alph_ids) {
  "use strict";
  var i, value;
  var array = [];
  var sum = 0;
  for (i = 0; i < alph_ids.length; i++) {
    value = probs_node.getAttribute(alph_ids[i]);
    if (value == null || !/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(value)) {
      throw new MotifError("DREME probability attribute " + alph_ids[i] + 
          " contains a value that is not a number \"" + value + "\".");
    }
    value = +value;
    if (value < 0 || value > 1) {
      throw new MotifError("DREME probability attribute " + alph_ids[i] + 
          " contains a value that is not a probabilty \"" + value + "\".");
    }
    sum += value;
    array.push(value);
  }
  if (Math.abs(sum - 1.0) > 0.01) {
    throw new MotifError("DREME alphabet_array contains probabilities that do not sum to 1.0.");
  }
  return array;
};

XmlMotifParser.prototype._streme_probability_array = function (probs_node, alph_ids) {
  "use strict";
  var i, value;
  var array = [];
  var sum = 0;
  for (i = 0; i < alph_ids.length; i++) {
    value = probs_node.getAttribute(alph_ids[i]);
    if (value == null || !/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(value)) {
      throw new MotifError("STREME probability attribute " + alph_ids[i] + 
          " contains a value that is not a number \"" + value + "\".");
    }
    value = +value;
    if (value < 0 || value > 1) {
      throw new MotifError("STREME probability attribute " + alph_ids[i] + 
          " contains a value that is not a probabilty \"" + value + "\".");
    }
    sum += value;
    array.push(value);
  }
  if (Math.abs(sum - 1.0) > 0.01) {
    throw new MotifError("STREME alphabet_array contains probabilities that do not sum to 1.0.");
  }
  return array;
};

XmlMotifParser.prototype._dreme_motif = function (motif_node, alph_ids) {
  "use strict";
  var i;
  // read the id
  var id = motif_node.getAttribute("seq");
  if (id == null) throw new MotifError("DREME XML motif seq attribute is missing.");
  var alt = "DREME";
  // read the length
  var len = motif_node.getAttribute("length");
  if (len == null || !/^\d+$/.test(len)) throw new  MotifError("DREME XML motif length attribute is not a number.");
  len = +len;
  if (len < 1) throw new MotifError("DREME XML motif length attribute should be at least 1.");
  // read the nsites
  var nsites = motif_node.getAttribute("nsites");
  if (nsites == null || !/^\d+$/.test(nsites)) throw new MotifError("DREME XML motif nsites attribute is not a number.");
  nsites = +nsites;
  if (nsites < 2) throw new MotifError("DREME XML motif nsites attribute should be at least 2.");
  // read the evalue
  var evalue = motif_node.getAttribute("evalue");
  if (evalue == null || !/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(evalue)) {
    throw new MotifError("DREME XML motif does not have valid evalue attribute.");
  }
  evalue = +evalue;
  // read the pwm
  var pwm = [];
  var probs = motif_node.querySelectorAll("pos");
  for (i = 0; i < probs.length; i++) {
    pwm.push(this._dreme_probability_array(probs[i], alph_ids));
  }
  if (pwm.length != len) {
    throw new MotifError("Expected length of pwm to match the stated length of the motif.");
  }
  // send the motif data
  if (typeof this.handler.motif == "function") {
    this.handler.motif({
      "id": id,
      "alt": alt,
      "len": len,
      "nsites": nsites,
      "evalue": evalue,
      "pwm": pwm
    });
  }
  if (this.rating < 10) this.rating = 10;
};

XmlMotifParser.prototype._streme_motif = function (motif_node, alph_ids) {
  "use strict";
  var i;
  // read the id
  var id = motif_node.getAttribute("id");
  if (id == null) throw new MotifError("STREME XML motif id attribute is missing.");
  var alt = "STREME";
  // read the length
  var len = motif_node.getAttribute("width");
  if (len == null || !/^\d+$/.test(len)) throw new  MotifError("STREME XML motif width attribute is not a number.");
  len = +len;
  if (len < 1) throw new MotifError("STREME XML motif length attribute should be at least 1.");
  // read the nsites
  var nsites = motif_node.getAttribute("train_pos_count");
  if (nsites == null || !/^\d+$/.test(nsites)) throw new MotifError("STREME XML motif train_pos_count attribute is not a number.");
  nsites = +nsites;
  if (nsites < 1) throw new MotifError("STREME XML motif train_pos_count attribute should be at least 1.");
  // read the evalue
  var pvalue = motif_node.getAttribute("test_pvalue");
  if (pvalue == null || !/^[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(pvalue)) {
    throw new MotifError("STREME XML motif does not have valid test_pvalue attribute.");
  }
  pvalue = +pvalue;
  // read the pwm
  var pwm = [];
  var probs = motif_node.querySelectorAll("pos");
  for (i = 0; i < probs.length; i++) {
    pwm.push(this._dreme_probability_array(probs[i], alph_ids));
  }
  if (pwm.length != len) {
    throw new MotifError("Expected length of pwm to match the stated length of the motif.");
  }
  // send the motif data
  if (typeof this.handler.motif == "function") {
    this.handler.motif({
      "id": id,
      "alt": alt,
      "len": len,
      "nsites": nsites,
      "evalue": pvalue,
      "pwm": pwm
    });
  }
  if (this.rating < 10) this.rating = 10;
};

XmlMotifParser.prototype.process_dreme_doc = function (doc) {
  "use strict";
  var i, version, alphabet, norc, background;
  // get the version
  var dreme_ele = doc.querySelector("dreme");
  if (dreme_ele == null) throw new MotifError("DREME XML missing /dreme .");
  version = dreme_ele.getAttribute("version");
  if (version == null) throw new MotifError("DREME XML missing /dreme/@version .");
  if (!/^\d+(?:\.\d+){0,2}$/.test(version)) {
    throw new MotifError("DREME XML value /dreme/@version does not match expected pattern; got \"" + version + "\".");
  }
  this.rating = 2;
  // get the alphabet
  var alphabet_ele = doc.querySelector("dreme > model > alphabet");
  var alphabet_ids;
  if (alphabet_ele != null) {
    var letter_elems = doc.querySelectorAll("dreme > model > alphabet > letter");
    var letter_id, letter_symbol, letter_aliases, letter_name, letter_colour, letter_complement, letter_equals;
    var sym2id = {};
    try {
      var alph_parser = new AlphabetParser();
      alph_parser.parse_header(alphabet_ele.getAttribute("name"));
      for (i = 0; i < letter_elems.length; i++) {
        letter_id = letter_elems[i].getAttribute("id");
        if (letter_id == null) throw new MotifError("DREME XML missing value /dreme/model/alphabet/letter[" + i + "]/@id .");
        letter_symbol = letter_elems[i].getAttribute("symbol");
        if (letter_symbol == null) throw new MotifError("DREME XML missing value /dreme/model/alphabet/letter[" + i + "]/@symbol .");
        sym2id[letter_symbol] = letter_id;
        letter_aliases = letter_elems[i].getAttribute("aliases");
        letter_name = letter_elems[i].getAttribute("name");
        letter_colour = letter_elems[i].getAttribute("colour");
        letter_complement = letter_elems[i].getAttribute("complement");
        letter_equals = letter_elems[i].getAttribute("equals");
        alph_parser.parse_symbol(letter_symbol, {
          "aliases": letter_aliases, "name": letter_name, "colour": letter_colour,
          "complement": letter_complement, "comprise": letter_equals
        });
      }
      var alph_data = alph_parser.parse_done();
      alphabet = new Alphabet(alph_data);
    } catch (e) {
      if (e instanceof AlphabetError) {
        throw new MotifError(e.message, e.reasons);
      } else {
        throw e;
      }
    }
    alphabet_ids = [];
    for (i = 0; i < alphabet.get_size_core(); i++) {
      alphabet_ids.push(sym2id[alphabet.get_symbol(i)]);
    }
  } else {
    alphabet = AlphStd.DNA;
    alphabet_ids = ["A", "C", "G", "T"];
  }
  this.rating = 3;
  // early versions of DREME did not have norc so check it exists first
  var norc_ele = doc.querySelector("dreme > model > norc");
  if (norc_ele == null) {
    norc = true;
  } else if (/^(TRUE|FALSE)$/.test(norc_ele.textContent)) {
    norc = (norc_ele.textContent == "TRUE");
  } else {
    throw new MotifError("DREME XML contains an unrecognized value for the norc \"" + norc + "\".");
  }
  this.rating = 4;
  // get the background
  var background_ele = doc.querySelector('dreme > model > background');
  background = this._dreme_probability_array(background_ele, alphabet_ids);
  this.rating = 5;
  // send the meta data
  if (typeof this.handler.meta == "function") {
    this.handler.meta({
      "version": version,
      "alphabet": alphabet,
      "norc": norc,
      "background": background
    });
  }
  // now iterate over motifs
  var motifs = doc.querySelectorAll("dreme  > motifs > motif");
  for (i = 0; i < motifs.length; i++) {
    this._dreme_motif(motifs[i], alphabet_ids);
  }
};

XmlMotifParser.prototype.process_streme_doc = function (doc) {
  "use strict";
  var i, version, alphabet, norc, background;
  // get the version
  var streme_ele = doc.querySelector("STREME");
  if (streme_ele == null) throw new MotifError("STREME XML missing /STREME.");
  version = streme_ele.getAttribute("version");
  if (version == null) throw new MotifError("STREME XML missing /STREME/@version .");
  if (!/^\d+(?:\.\d+){0,2}$/.test(version)) {
    throw new MotifError("STREME XML value /STREME/@version does not match expected pattern; got \"" + version + "\".");
  }
  this.rating = 2;
  // get the alphabet
  var alphabet_ele = doc.querySelector("STREME > model > alphabet");
  var alphabet_ids;
  if (alphabet_ele != null) {
    var letter_elems = doc.querySelectorAll("STREME > model > alphabet > letter");
    var letter_id, letter_symbol, letter_aliases, letter_name, letter_colour, letter_complement, letter_equals;
    var sym2id = {};
    try {
      var alph_parser = new AlphabetParser();
      alph_parser.parse_header(alphabet_ele.getAttribute("name"));
      for (i = 0; i < letter_elems.length; i++) {
        letter_id = letter_elems[i].getAttribute("id");
        if (letter_id == null) throw new MotifError("STREME XML missing value /STREME/model/alphabet/letter[" + i + "]/@id .");
        letter_symbol = letter_elems[i].getAttribute("symbol");
        if (letter_symbol == null) throw new MotifError("STREME XML missing value /STREME/model/alphabet/letter[" + i + "]/@symbol .");
        sym2id[letter_symbol] = letter_id;
        letter_aliases = letter_elems[i].getAttribute("aliases");
        letter_name = letter_elems[i].getAttribute("name");
        letter_colour = letter_elems[i].getAttribute("colour");
        letter_complement = letter_elems[i].getAttribute("complement");
        letter_equals = letter_elems[i].getAttribute("equals");
        alph_parser.parse_symbol(letter_symbol, {
          "aliases": letter_aliases, "name": letter_name, "colour": letter_colour,
          "complement": letter_complement, "comprise": letter_equals
        });
      }
      var alph_data = alph_parser.parse_done();
      alphabet = new Alphabet(alph_data);
    } catch (e) {
      if (e instanceof AlphabetError) {
        throw new MotifError(e.message, e.reasons);
      } else {
        throw e;
      }
    }
    alphabet_ids = [];
    for (i = 0; i < alphabet.get_size_core(); i++) {
      alphabet_ids.push(sym2id[alphabet.get_symbol(i)]);
    }
  } else {
    alphabet = AlphStd.DNA;
    alphabet_ids = ["A", "C", "G", "T"];
  }
  this.rating = 3;
  // make a list of the letter_ids in order of the alphabet
  var symbol_map = {};
  for (i = 0; i < letter_elems.length; i++) {
    letter_symbol = letter_elems[i].getAttribute("symbol");
    if (letter_symbol == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[" + i + "]/@symbol .");
    letter_id = letter_elems[i].getAttribute("id");
    if (letter_id == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[" + i + "]/@id .");
    symbol_map[letter_symbol] = letter_id;
  }
  var symbol_ids = [];
  for (i = 0; i < alphabet.get_size_core(); i++) {
    var id = symbol_map[alphabet.get_symbol(i)];
    if (id == null) throw new MotifError("MEME XML missing value /MEME/training_set/alphabet/letter[@symbol='" + alphabet.get_symbol(i) + "']/@id .");
    symbol_ids.push(id);
  }
  this.rating = 4;
  // try to get the background_frequencies (post 5.3.0)
  var background_ele = doc.querySelector("STREME > model > background_frequencies > alphabet_array");
  if (background_ele != null) {
    background = this._meme_alphabet_array(background_ele, symbol_ids, true);
  } else {
    // try to get the background (removed after 5.3.0)
    var background_ele = doc.querySelector('STREME > model > background');
    if (background_ele != null) {
      background = this._streme_probability_array(background_ele, alphabet_ids);
    } else {
      throw new MotifError("STREME XML missing /STREME/model/background_frequencies/alphabet_array .");
    }
  }
  this.rating = 5;
  // send the meta data
  if (typeof this.handler.meta == "function") {
    this.handler.meta({
      "version": version,
      "alphabet": alphabet,
      "norc": true,
      "background": background
    });
  }
  // now iterate over motifs
  var motifs = doc.querySelectorAll("STREME > motifs > motif");
  for (i = 0; i < motifs.length; i++) {
    this._streme_motif(motifs[i], alphabet_ids);
  }
};

XmlMotifParser.prototype.process_doc = function (doc) {
  "use strict";
    try {
      // IE does not support evaluate so we have to do things the messy way
      if (doc.querySelector("MEME")) {
        this.process_meme_doc(doc);
      } else if (doc.querySelector("dreme")) {
        this.process_dreme_doc(doc);
      } else if (doc.querySelector("STREME")) {
        this.process_streme_doc(doc);
      } else {
        throw new MotifError("Not a recognized XML motif format.");
      }
    } catch (e) {
      if (e instanceof MotifError) {
        this.error_list.push({"is_error": true, "message": e.message, "reasons": e.reasons});
      } else {
        throw e;
      }
    }
};

var MotifParser = function(handler) {
  "use strict";
  var me, h2;
  this.handler = handler;
  this.stop_xml = false;
  this.motif_count = 0;
  me = this;
  h2 = {
    "meta": function (info) {
      me.handler.meta(info);
    },
    "motif": function (motif_data) {
      me.stop_xml = true;
      me.motif_count++;
      me.handler.motif(motif_data);
    }
  };
  this.html_parser = new HtmlMotifParser(h2);
  this.text_parser = new TextMotifParser(h2);
  this.xml_parser = new XmlMotifParser(h2);
  this.streamers = [this.html_parser, this.text_parser];
};

MotifParser.prototype.get_errors = function () {
  "use strict";
  var i, parsers, best_rating, messages;
  best_rating = 0;
  messages = [{
    "is_error": true, 
    "message": "The file did not match any supported motif formats.", 
    "reasons": []
  }];
  parsers = [this.html_parser, this.text_parser, this.xml_parser];
  for (i = 0; i < parsers.length; i++) {
    if (parsers[i].parser_rating() <= 1) continue;
    if (parsers[i].parser_rating() > best_rating) {
      messages = parsers[i].get_errors();
      best_rating = parsers[i].parser_rating();
    }
  }
  if (this.motif_count == 0) {
    messages.push({
      "is_error": true,
      "message": "No motifs found.",
      "reasons": []
    });
  }
  return messages;
};

MotifParser.prototype._check_nulls = function (chunk) {
  "use strict";
  var i;
  for (i = 0; i < chunk.length; i++) {
    if (chunk[i] == 0x0) return true;
  }
  return false;
};

MotifParser.prototype._report_errors = function (error_list) {
  "use strict";
  var i, error;
  if (typeof this.handler.error !== "function") return;
  for (i = 0; i < error_list.length; i++) {
    error = error_list[i];
    this.handler.error(error["is_error"], error["message"], error["reasons"]);
  }
};

/*
 * If we fail to find motif by the streaming method
 * and we don't encounter anything that indicates 
 * it is binary then try to load the full file in as XML.
 */
MotifParser.prototype._process_blob2 = function (blob) {
  "use strict";
  var me, reader;
  me = this;
  reader = new FileReader();
  reader.onload = function(e) {
    "use strict";
    var parser, parseErrorNS, doc;
    console.log("File is loaded, now attempting to parse as XML.\n" +
      "This will throw an XML parsing error which can be safely ignored.\n");
    parser = new DOMParser();
    // this produces a warning about invalid XML, it can be safely ignored.
    // IE throws an error when the syntax is wrong so we have to put a try block
    try { parseErrorNS = parser.parseFromString('INVALID', 'text/xml').childNodes[0].namespaceURI; } catch (e) { }
    try {
      doc = null;
      doc = parser.parseFromString(e.target.result, "application/xml");
    } catch (e) {
      // Note that IE is the only one that throws a error on syntax problems
      // other browsers return a special structure
    }
    if (doc != null && doc.getElementsByTagNameNS(parseErrorNS, 'parsererror').length == 0) {
      me.xml_parser.process_doc(doc);
    }
    me._report_errors(me.get_errors());
    if (typeof me.handler.end == "function") me.handler.end();
  };
  reader.readAsText(blob, "UTF-8");
};

MotifParser.prototype.process_blob = function (blob, offset, chunk_size) {
  "use strict";
  var reader, me;
  if (this.give_up) return;
  // set default parameter values
  if (typeof offset === "undefined") offset = 0;
  if (typeof chunk_size === "undefined") chunk_size = 1 << 12;
  // update the progress
  if (offset === 0 && typeof this.handler.begin == "function")
    this.handler.begin(blob.size);
  if (typeof this.handler.progress == "function")
    this.handler.progress(offset / blob.size);
  // setup the reader
  me = this; // so we can access 'this' inside the closure
  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
    var chunk, i, file_name, format, message;
    if (me.give_up) return;
    // process the loaded chunk
    chunk = new Uint8Array(reader.result);
    if (offset == 0) {
      // check magic numbers
      if ((format = unusable_format(chunk, 40, blob.name)) != null) {
        switch (format.type) {
          case FileType.ENCODING:
            message = "This is encoded as " + format.name + " which is not parsable as a motif.";
            break;
          case FileType.BINARY:
          case FileType.COMPRESSED:
          default:
            message = "This is a " + format.name + " file which is not parsable as a motif.";
            break;
        }
        if (typeof me.handler.progress == "function") me.handler.progress(1.0);
        if (typeof me.handler.error == "function") me.handler.error(true, message);
        if (typeof me.handler.end == "function") me.handler.end();
        return;
      }
    } else if (me._check_nulls(chunk)) {
        message = "This is a unknown binary file which is not parsable as a motif.";
        if (typeof me.handler.progress == "function") me.handler.progress(1.0);
        if (typeof me.handler.error == "function") me.handler.error(true, message);
        if (typeof me.handler.end == "function") me.handler.end();
        return;
    }
    for (i = 0; i < me.streamers.length; i++) {
      try { me.streamers[i].process_chunk(chunk); } catch (e) {
        // remove the parser which errored
        me.streamers.splice(i, 1); i--;
      }
    }
    if ((offset + chunk_size) >= blob.size) {
      // at the end, notify all remaining parsers
      for (i = 0; i < me.streamers.length; i++) {
        try { me.streamers[i].process_end(); } catch (e) {
          // remove the parser which errored
          me.streamers.splice(i, 1); i--;
        }
      }
      // update the streaming progress
      if (typeof me.handler.progress == "function") me.handler.progress(1.0);
      // if we've ruled out xml then finish, otherwise call the xml parser
      if (me.stop_xml) {
        if (typeof me.handler.end == "function") me.handler.end(me.get_errors());
      } else {
        // try parsing as XML
        me._process_blob2(blob);
      }
    } else {
      // check if any of the parser still work
      if (me.streamers.length > 0) {
        // start loading the next chunk
        me.process_blob(blob, offset + chunk_size, chunk_size);
      } else {
        // all the streaming parsers ran into errors 
        if (typeof me.handler.progress == "function") me.handler.progress(1.0);
        // try parsing as XML
        if (!me.stop_xml) me._process_blob2(blob);
      }
    }
  };
  // read the next chunk
  reader.readAsArrayBuffer(blob.slice(offset, offset + chunk_size));
};

MotifParser.prototype.cancel = function () {
  "use strict";
  this.give_up = true;
  this.handler = {};
};

/*
 * An enumeration of attributes a token in a simple motif can have.
 *
 * SPACE - the token contains only whitespace
 * LINE - the token contains newline characters (often combined with SPACE)
 * MATRIX - the token is a number and probably part of a matrix motif
 * IUPAC - the token is a IUPAC code and probably part of a IUPAC motif
 * WARN - the token should be highlighted to signify a warning
 * ERROR - the token should be highlighted to signify an error
 * EAST - the token's indicator text should be shown on the east side of the token
 * SOUTH - the token's indicator text should be shown on the south side of the token
 */
var SimpleMotifTokEn = {"SPACE":1, "LINE":2, "MATRIX":4, "IUPAC":8,
  "WARN":16, "ERROR":32, "EAST": 64, "SOUTH" : 128};
if (Object.freeze) Object.freeze(SimpleMotifTokEn);

/*
 * Stores a chunk of text that makes up a simple motif along with the type of
 * content and a optional indicator to give hints about problems.
 */
var SimpleMotifToken = function(text, type) {
  this.text = text;
  this.type = type;
  this.indicator = "";
};

var SimpleMotif = function(text, alphabet, motif_num) {
  "use strict";
  // initilize object vars
  this.text = text;
  this.alphabet = alphabet;
  this.tokens = [];
  // generic motif values
  this.id = null;
  this.alt = null;
  this.len = null;
  this.nsites = null;
  this.evalue = 0; // we can't determine an E-value for a simple motif
  this.pwm = null;
  this.psm = null;

  // separate off the first line and see if the user defined a name for the motif
  var motif_text = text;
  var name_tokens = [];
  var first_line_end = text.indexOf("\n") + 1;
  var first_line = text.substring(0, first_line_end);
  if ((match = /^\s*>(\S+)(?:\s+(\S+))?/.exec(first_line)) != null) {
    this.id = match[1];
    this.alt = match[2];
    name_tokens.push(new SimpleMotifToken(first_line, 0));
    motif_text = text.substring(first_line_end);
  }

  // RULES
  // What if an alphabet contains numbers?...
  // How do we tell if a motif is a matrix of numbers or an iupac motif?
  // 1) Does the line only contain valid characters for an IUPAC expression?
  // 2) Does the line only contain possible numbers?
  // 3) What if it fits both criteria? -> use spaces in middle as differenciator?
  
  var i, match, last, has_nl, token;
  var word_re = /(?:(\s+)|(\S+))/g;
  var matrix_tokens = [];
  var matrix_errors = 0;
  var matrix_numbers = 0;
  var iupac_tokens = [];
  var iupac_errors = 0;
  var midline_gaps = 0;
  // parse the text attempting both possible ways, meanwhile collecting information
  while ((match = word_re.exec(motif_text)) != null) {
    if (match[1] != null) { // space
      has_nl = (match[0].indexOf('\n') != -1);
      token = new SimpleMotifToken(match[0],
          SimpleMotifTokEn.SPACE | (has_nl ? SimpleMotifTokEn.LINE : 0));
      if (matrix_tokens.length >= 1 && !has_nl) {
        midline_gaps++;
      }
      matrix_tokens.push(token);
      iupac_tokens.push(token);
    } else if (match[2] != null) {
      // parse as if IUPAC only
      for (i = 0; i < match[0].length; i++) {
        if (alphabet.get_index(match[0].charAt(i)) != -1 || match[0].charAt(i) == '[' || match[0].charAt(i) == ']') {
          iupac_tokens.push(new SimpleMotifToken(match[0].charAt(i), SimpleMotifTokEn.IUPAC));
        } else {
          iupac_tokens.push(new SimpleMotifToken(match[0].charAt(i), SimpleMotifTokEn.ERROR));
          iupac_errors++;
        }
      }
      // parse as if matrix only
      if (/^[+]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?$/.test(match[0])) {
        matrix_tokens.push(new SimpleMotifToken(match[0], SimpleMotifTokEn.MATRIX));
        matrix_numbers++;
      } else {
        matrix_tokens.push(new SimpleMotifToken(match[0], SimpleMotifTokEn.ERROR));
        matrix_errors++;
      }
    }
  }
  var pair;
  var is_iupac = false;
  if (iupac_errors == 0 && matrix_errors == 0) {
    is_iupac = (matrix_numbers == 0 || midline_gaps == 0);
  } else if (iupac_errors == 0) {
    is_iupac = true;
  } else if (matrix_errors == 0) {
    is_iupac = false;
  } else {
    is_iupac = (matrix_numbers == 0 || midline_gaps == 0);
  }
  this.tokens.push.apply(this.tokens, name_tokens);
  if (is_iupac) {
    this.tokens.push.apply(this.tokens, iupac_tokens);
    pair = SimpleMotif._parse_iupac(alphabet, this.tokens);
  } else {
    this.tokens.push.apply(this.tokens, matrix_tokens);
    pair = SimpleMotif._parse_matrix(alphabet, this.tokens);
  }

  this.pwm = pair["pwm"];
  this.nsites = pair["nsites"];
  this.len = this.pwm.length;
  if (this.id == null) {
    this.id = "" + motif_num;
    this.alt = MotifUtils.find_codes(this.pwm, MotifUtils._code_freqs(alphabet));
  }
  this.psm = MotifUtils.freqs_to_scores(this.pwm, this.nsites);
};


/*
 * Static Function
 * Convert a set of tokens into a pwm and a site count.
 */
SimpleMotif._parse_iupac = function(alphabet, tokens) {
  "use strict";
  var i, j, k, token, any, lines, line, group_start, group_value, value;
  var expected, last;
  // make a value for any
  any = [];
  for (i = 0; i < alphabet.get_size_core(); i++) any.push(i);
  // convert into a matrix of positions
  lines = [];
  line = [];
  group_start = -1;
  expected = -1;
  last = -1;
  for (i = 0; i < tokens.length; i++) {
    token = tokens[i];
    if ((token.type & SimpleMotifTokEn.LINE) != 0) {
      if (group_start != -1) {
        if (group_value.length > 0) line.push(group_value);
        // a group was started but not finished
        tokens[group_start].type |= SimpleMotifTokEn.ERROR;
        // reset
        group_start = -1;
      }
      if (lines.length == 0) {
        expected = line.length;
      } else if (line.length < expected) {
        tokens[last].indicator = "\u2605"; // filled 5 point star
        tokens[last].type |= SimpleMotifTokEn.EAST;
      }
      lines.push(line);
      line = [];
      continue;
    } else if ((token.type & SimpleMotifTokEn.MATRIX) != 0) {
      token.type |= SimpleMotifTokEn.ERROR;
    }
    if ((token.type & SimpleMotifTokEn.IUPAC) == 0) continue;
    // an IUPAC token or a bracket past this point
    last = i;
    // mark anything past the expected size as an error
    if (expected != -1 && line.length >= expected) {
      token.type |= SimpleMotifTokEn.ERROR;
    }
    if (token.text == '[') {
      if (group_start == -1) {
        group_start = i;
        group_value = [];
      } else {
        token.type |= SimpleMotifTokEn.ERROR;
      }
    } else if (token.text == ']') {
      if (group_start != -1) {
        group_start = -1;
        if (group_value.length > 0) line.push(group_value);
      } else {
        token.type |= SimpleMotifTokEn.ERROR;
      }
    } else {
      value = alphabet.get_indexes(token.text);
      if (typeof value === "undefined") { // treat as any
        token.type |= SimpleMotifTokEn.ERROR;
        value = any;
      }
      if (group_start == -1) {
        line.push(value);
      } else {
        // merge with group value
        j = 0; k = 0;
        while (j < group_value.length && k < value.length) {
          if (group_value[j] < value[k]) {
            j++;
          } else if (group_value[j] > value[k]) {
            group_value.splice(j, 0, value[k]);
            j++; k++;
          } else {
            j++; k++;
          }
        }
        for (;k < value.length; k++) {
          group_value.push(value[k]);
        }
      }
    }
  }
  if (group_start != -1) {
    if (group_value.length > 0) line.push(group_value);
    // a group was started but not finished
    tokens[group_start].type |= SimpleMotifTokEn.ERROR;
  }
  if (lines.length > 0 && line.length < expected && last != -1) {
    tokens[last].indicator = "\u2605"; // filled 5 point star
    tokens[last].type |= SimpleMotifTokEn.EAST;
  }
  lines.push(line);
  // convert into a probability matrix
  var size = lines[0].length;
  var pwm = [];
  for (i = 0; i < size; i++) {
    pwm[i] = [];
    // initialise pwm with pseudocount
    for (j = 0; j < alphabet.get_size_core(); j++) {
      pwm[i][j] = 0;
    }
    // sum up counts from each line
    var counts = 0;
    for (j = 0; j < lines.length; j++) {
      line = lines[j];
      if (line.length <= i) continue;
      value = line[i];
      var count = 1.0 / value.length;
      for (k = 0; k < value.length; k++) {
        pwm[i][value[k]] += count;
      }
      counts++;
    }
    // normalise
    for (j = 0; j < alphabet.get_size_core(); j++) {
      pwm[i][j] /= counts;
    }
  }
  return {"pwm": pwm, "nsites": lines.length};
};

SimpleMotif._parse_matrix = function(alphabet, tokens) {
  "use strict";
  var i, j, token, lines, line, row_matrix, get_alen, get_entry;
  var size, entry_alen, full, pwm, counts, expected_counts, diff_counts, delta;
  var site_counts;
  var alen = alphabet.get_size_core();
  // sort matrix tokens into lines
  lines = [];
  line = [];
  for (i = 0; i < tokens.length; i++) {
    token = tokens[i];
    if ((token.type & SimpleMotifTokEn.LINE) != 0) {
      lines.push(line);
      line = [];
    } else if ((token.type & SimpleMotifTokEn.IUPAC) != 0) {
      token.type |= SimpleMotifTokEn.ERROR; // mark IUPAC tokens as illegal
    }
    if ((token.type & SimpleMotifTokEn.MATRIX) == 0) continue;
    // matrix tokens only past here
    line.push(token);
  }
  lines.push(line);

  // determine the orientation
  if (lines[0].length == alen) {
    // row matricies are preferred
    row_matrix = true;
  } else if (lines.length == alen) {
    row_matrix = false;
  } else {
    // uhoh, it doesn't seem to match the alphabet...
    // pick the closest orientation that is still smaller, 
    // if they're equal then assume it's a column matrix
    // because they're more likely to be typing more rows
    // then more columns.
    if (lines[0].length < alen && lines.length < alen) {
      if (lines[0].length > lines.length) {
        row_matrix = true;
      } else {
        row_matrix = false;
      }
    } else if (lines[0].length < alen) {
      row_matrix = true;
    } else if (lines.length < alen) {
      row_matrix = false;
    } else {
      // well there goes the theory that they're still typing it...
      // maybe they have the wrong alphabet selected?
      row_matrix = true;
    }
  }

  // generate pwm
  get_alen = (row_matrix ? SimpleMotif._row_alen : SimpleMotif._col_alen);
  get_entry = (row_matrix ? SimpleMotif._row_accessor : SimpleMotif._col_accessor);
  size = (row_matrix ? lines.length : lines[0].length);
  pwm = [];
  expected_counts = -1;
  delta = 0.1;
  for (i = 0; i < size; i++) {
    pwm[i] = [];
    counts = 0;
    // copy counts
    full = true;
    for (j = 0; j < alen; j++) {
      token = get_entry(lines, i, j);
      if (token != null) {
        pwm[i][j] = +(token.text);
      } else {
        pwm[i][j] = 0;
        full = false;
      }
      counts += pwm[i][j];
    }
    entry_alen = get_alen(lines, i);
    if (full && entry_alen == alen) {
      // check that the counts are the same
      if (expected_counts == -1) {
        expected_counts = counts;
      } else if (Math.abs(expected_counts - counts) > delta) {
        for (j = 0; j < alen; j++) {
          get_entry(lines, i, j).type |= SimpleMotifTokEn.WARN;
        }
        token = get_entry(lines, i, alen-1);
        token.indicator = SimpleMotif._count_format(counts - expected_counts);
        token.type |= (row_matrix ? SimpleMotifTokEn.EAST : SimpleMotifTokEn.SOUTH);
      }
    } else if (entry_alen > alen) {
      // check for extra numbers
      for (j = alen; j < entry_alen; j++) {
        token = get_entry(lines, i, j);
        if (token != null) token.type |= SimpleMotifTokEn.ERROR;
      }
    } else {
      token = get_entry(lines, i, entry_alen-1);
      token.indicator = "\u2605"; // filled 5 point star
      token.type |= (row_matrix ? SimpleMotifTokEn.EAST : SimpleMotifTokEn.SOUTH);
    }
    // check for zero counts
    if (counts == 0) {
      // evenly distribute a pseudocount
      // to ensure the pspm is not zero
      for (j = 0; j < alen; j++) {
        pwm[i][j] = 1.0 / alen;
      }
    } else {
      // normalise
      for (j = 0; j < alen; j++) {
        pwm[i][j] /= counts;
      }
    }
  }
  if (!row_matrix) {
    for (i = 1; i < lines.length; i++) {
      for (j = size; j < lines[i].length; j++) {
        lines[i][j].type |= SimpleMotifTokEn.ERROR;
      }
    }
  }
  site_counts = Math.max(1, Math.round(expected_counts));
  return {"pwm": pwm, "nsites": site_counts};
};

SimpleMotif._row_alen = function(lines, entry_pos) {
  return lines[entry_pos].length;
};

SimpleMotif._col_alen = function(lines, entry_pos) {
  var i;
  for (i = lines.length -1; i >= 0; i--) {
    if (lines[i].length > entry_pos) break;
  }
  return i+1;
};

SimpleMotif._row_accessor = function(lines, entry_pos, alpha_pos) {
  if (entry_pos >= lines.length) return null;
  if (alpha_pos >= lines[entry_pos].length) return null;
  return lines[entry_pos][alpha_pos];
};

SimpleMotif._col_accessor = function(lines, entry_pos, alpha_pos) {
  if (alpha_pos >= lines.length) return null;
  if (entry_pos >= lines[alpha_pos].length) return null;
  return lines[alpha_pos][entry_pos];
};

SimpleMotif._count_format = function(difference) {
   return (difference > 0 ? "+" : "") + difference.toFixed(1);
};
