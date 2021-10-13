var BedCheckerUtil = {};

//******************************************************************************
// BED format Checker
//******************************************************************************
var BedChecker = function (handler) {
  "use strict";
  var i;
  // store a reference to the handler
  this.handler = handler;
  // current parsing function
  this.process = this._process_start;
  // abort flag
  this.give_up = false;
};

BedChecker.prototype._process_start = function (code, type) {
  "use strict";
  return true;
};

// When we're done, call the approprate functions on the handler
BedChecker.prototype._signal_stop = function() {
  if (typeof this.handler.end == "function") this.handler.end();
};

//******************************************************************************
// Public functions
//******************************************************************************

BedChecker.prototype.process_file = function (file_input) {
  "use strict";
  var err, message, num_comments, reader;
  err = false;
  num_comments = 0;
  if (this.give_up) return;
  reader = new FileReader();
  reader.onload = function(evt) {
    "use strict";
     var lines = (this.result).trim().split('\n');
     for(var line = 0; line < lines.length; line++) {
        if (lines[line].charAt(0) != '#') {
          var tabs = lines[line].split('\t');
          if (tabs.length < 3 || tabs.length > 12) {
            err = true;
            message = "It has " + tabs.length + " tab delimited fields.\n"
            message += "It should have between 3 and 9."
            break;
          }
          var start = Number(tabs[1].trim());
          var stop = Number(tabs[2].trim());
            // The 2nd and 3rd fields thould be positive integers
          if ( isNaN(start) || (Math.floor(start) != start) || start <= 0 || 
               isNaN(stop) || (Math.floor(stop)  != stop) || stop <= 0) {
            err = true;
            message = "Fields 2 and 3 should be positive integers."
            break;
          }
          if (start > stop) {
            err = true;
            message = "Fields 2 and 3 should be positive integers and\n"
            message += "field 3 should be greater then field 2."
            break;
          }
      } else {
        num_comments += 1;
      }
    }
    if (lines.length === 0 || num_comments === lines.length) {
            err = true;
            message = "The file doesn't contain any BED data.\n"
    }
    if (err) {
      alert("The file " + file_input.files[0].name +" doesn't seem to be a valid BED file.\n" +
       "Invalid format found on line " + (line + 1) + ".\n" + message);
      file_input.value = "";
      return;
    }
  };
  if (file_input.files[0]) {
    reader.readAsText(file_input.files[0]);
  }
};

BedChecker.prototype.cancel = function () {
  "use strict";
  this.give_up = true;
  this.handler = {};
};

//******************************************************************************
// Bed Handler
//******************************************************************************
var BedHandler = function () {
  this.reset();
};

BedHandler.prototype.reset = function () {
  // have the file details changed?
  this.updated = false;
  // the part of the file processed
  this.fraction = 0;
  // fasta details
  this.file_size = 0;
  this.file_symbols = "";
  // keep track of problems found
  //this.missing_name = new FileFaults();
};

BedHandler.prototype.summary = function () {
  "use strict";
  var error, warning, messages, reason, reasons, letters, add;
  var help;
  // setup
  error = false;
  warning = false;
  messages = [];
  // create closure to add messages
  add = function(is_error, message, reasons) {
    "use strict";
    messages.push({"is_error": is_error, "message": message, "reasons": reasons});
    if (is_error) error = true;
    else warning = true;
  };
  // file size warning
  if (this.max_file_size != null && this.file_size > this.max_file_size) {
    add(false, "Large file. ", ["File is " + (Math.round(this.file_size / (1<<20) )) + "MB"] + ". ")
  }
  // clear updated state
  this.updated = false;
  // return state
  return {"error": error, "warning": warning, "messages": messages};
};

// Reading of the file has begun
BedHandler.prototype.begin = function (file_size) {
  "use strict";
  this.reset();
  this.file_size = file_size;
  this.updated = true;
};

// Reading of the file has finished (perhaps early due to an error)
BedHandler.prototype.end = function () {
  "use strict";
  this.updated = true;
};

// Parsing has stopped due to an unreadable file format
BedHandler.prototype.error_format = function (type, name) {
  "use strict";
  this.unusable_format_type = type;
  this.unusable_format_name = name;
  this.updated = true;
};

