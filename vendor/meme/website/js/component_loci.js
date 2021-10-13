/******************************************************************************
 * Construct the manager of the loci input.
 ******************************************************************************/
var LociInput = function(container) {
  "use strict";
  var me, box;
  // make 'this' accessable in inner scopes
  me = this;
  // store the parameters
  this.container = container;
  // get the file related parts
  this.file_surround = this.container.querySelector("span.loci_file");
  this.file_input = this.file_surround.querySelector("input");
  this.file_popup = this.file_surround.querySelector("div.popup");
  // create a list of submittable fields so we can disable the ones we are not using
  this.submittables = [this.file_input];
  // other things
  this.parser = null;
  this.timer = null;
  // initialise
  this._file_update();
  // add listeners
  // detect file changes
  this.file_input.addEventListener('change', function() {
    me._file_update();
  }, false);
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype.check = function() {
  "use strict";
  var id, exists, handler, summary;
  exists = this.file_input.value.length > 0;
  return true;
};

/******************************************************************************
 *
 ******************************************************************************/
// reset the fields to a consistant state
LociInput.prototype.reset = function() {
  "use strict";
  if (this.parser) this.parser.cancel();
  this.parser = null;
  if (this.timer) window.clearTimeout(this.timer);
  this.timer = null;
  this.file_input.value = "";
  this.file_dbh.reset();
  this.text_area.value = "";
  this.text_dbh.reset();
  if (this.db_listing != null) this._clear_select(this.db_listing);
  if (this.db_version != null) this._clear_select(this.db_version);
};


LociInput.prototype._fire_loci_checked_event = function() {
  "use strict";
  var me;
  me = this;
  try {
    // IE sometimes has problems with this line.
    // I think they are related to the page not being fully loaded.
    this.container.dispatchEvent(new CustomEvent("fire_checked", {detail: {controler: me}}));
  } catch (e) {
    if (e.message && e.name && window.console) {
      console.log("Suppressed exception " + e.name + ": " + e.message);
    }
  }
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype._file_update = function() {
  "use strict";
   this._file_validate();
};

/******************************************************************************
 *
 ******************************************************************************/
LociInput.prototype._file_validate = function() {
  var file;
  if (this.parser) this.parser.cancel();
  // file = this.file_input.files[0];
  this.parser = new BedChecker(this.file_dbh);
  this.parser.process_file(this.file_input);
};
