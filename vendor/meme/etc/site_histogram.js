var SiteHistogramGraphMetrics = function(top_edge, bottom_edge,
    left_mark, left_val, right_mark, right_val) {
  this.top_edge = top_edge;
  this.bottom_edge = bottom_edge;
  this.left_mark = left_mark;
  this.left_val = left_val;
  this.right_mark = right_mark;
  this.right_val = right_val;
};

/*
 * SiteHistogramGraph
 */
var SiteHistogramGraph = function(counts) {
  "use strict";
  this.counts = [];
  this.start = 0;
  this.end = counts.length;
  this.ymax = 100;
  this.yinc = 20;
  this.ydigits = 0

  // Normalize the counts.
  var i, total=0;
  for (i=0; i<counts.length; i++) total += counts[i];
  for (i=0; i<counts.length; i++) this.counts[i] =  counts[i] *= 100/total;
} // SiteHistogramGraph

/*
 * draw_graph
 *
 * Draws a site count histogram.
 */
SiteHistogramGraph.prototype.draw_graph = function (ctx, w, h, x_label, y_label, color) {
  "use strict";
  var gap = 10;
  var l_margin = gap + 15 + 10 * (2 + this.ydigits);
  var t_margin = 10;
  var b_margin = 40;
  var r_margin = 30;

  // draw graph
  ctx.save();
  // draw border
  ctx.beginPath();
  ctx.moveTo(l_margin - 0.5, t_margin +0.5);
  ctx.lineTo(l_margin - 0.5, h - (b_margin - 0.5));
  ctx.lineTo(w - r_margin - 0.5, h - (b_margin - 0.5));
  ctx.lineTo(w - r_margin - 0.5, t_margin +0.5);
  ctx.closePath();
  ctx.stroke();
  // draw y axis
  ctx.save();
  ctx.translate(l_margin, t_margin);
  this.draw_y_axis(ctx, 0, h - (t_margin + b_margin));
  ctx.restore();
  // draw y labels
  ctx.save();
  ctx.translate(gap, t_margin);
  this.draw_y_axis_label(ctx, l_margin, h - (t_margin + b_margin), y_label);
  ctx.restore();
  // draw x axis
  ctx.save();
  ctx.translate(l_margin, h - b_margin);
  this.draw_x_axis(ctx, w - (l_margin + r_margin), b_margin);
  ctx.restore();
  // draw top axis
  ctx.save();
  ctx.translate(l_margin, t_margin);
  this.draw_top_axis(ctx, w - (l_margin + r_margin), 0);
  ctx.restore();
  // draw x axis labels
  ctx.save();
  ctx.translate(l_margin, h - b_margin);
  this.draw_x_axis_label(ctx, w - (l_margin + r_margin), b_margin - gap, x_label);
  ctx.restore();
  // draw bars
  ctx.save();
  ctx.translate(l_margin, t_margin);
  this.draw_bars(ctx, w - (l_margin + r_margin), h - (t_margin + b_margin), color);
  ctx.restore();
  return new SiteHistogramGraphMetrics(t_margin, h - b_margin,
    l_margin, this.start, w - r_margin, this.end);
} // SiteHistogramGraph.prototype.draw_graph

/*
 * draw_x_axis
 *
 */
SiteHistogramGraph.prototype.draw_x_axis = function(ctx, w, h) {
  "use strict";
  var scale_x, length, tic_inc, tic_min, tic_max, i, x;
  // pixels per unit
  scale_x = w / (this.end - this.start);
  // calculate a good tic increment
  length = this.end - this.start;
  if (length > 50) {
    // multiple of 10
    tic_inc = Math.max(1, Math.round(length / 100)) * 10;
  } else if (length > 25) {
    tic_inc = 5;
  } else if (length > 10) {
    tic_inc = 2;
  } else {
    tic_inc = 1;
  }
  // work out the min and max values within the start and end
  tic_min = Math.round(this.start / tic_inc) * tic_inc;
  if (tic_min < this.start) {
    tic_min += tic_inc;
  }
  tic_max = Math.round(this.end / tic_inc) * tic_inc;
  if (tic_max > this.end) {
    tic_max -= tic_inc;
  }
  // draw the tics
  ctx.save();
  ctx.font = (length > 10) ? "7pt Helvetica" : "9pt Helvetica";
  ctx.textBaseline = "top";
  ctx.textAlign = "center";
  for (i = tic_min; i <= tic_max; i+= tic_inc) {
    x = Math.round((i - this.start) * scale_x) + 0.5;
    ctx.fillText(""+i, x, 5);
    if (i == this.start || i == this.end) continue;
    ctx.beginPath();
    ctx.moveTo(x, -5);
    ctx.lineTo(x, 3);
    ctx.stroke();
  }
  ctx.restore();
};

/*
 * draw_top_axis
 *
 */
SiteHistogramGraph.prototype.draw_top_axis = function(ctx, w, h) {
  "use strict";
  var scale_x, length, tic_inc, tic_min, tic_max, i, x;
  // pixels per unit
  scale_x = w / (this.end - this.start);
  // calculate a good tic increment
  length = this.end - this.start;
  if (length > 50) {
    // multiple of 10
    tic_inc = Math.max(1, Math.round(length / 100)) * 10;
  } else if (length > 25) {
    tic_inc = 5;
  } else {
    tic_inc = 1;
  }
  // work out the min and max values within the start and end
  tic_min = Math.round(this.start / tic_inc) * tic_inc;
  if (tic_min < this.start) {
    tic_min += tic_inc;
  }
  tic_max = Math.round(this.end / tic_inc) * tic_inc;
  if (tic_max > this.end) {
    tic_max -= tic_inc;
  }
  // draw the tics
  ctx.save();
  ctx.font = "9pt Helvetica";
  ctx.textBaseline = "top";
  ctx.textAlign = "center";
  for (i = tic_min; i <= tic_max; i+= tic_inc) {
    if (i == this.start || i == this.end) continue;
    x = Math.round((i - this.start) * scale_x) + 0.5;
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, 8);
    ctx.stroke();
  }
  ctx.restore();
};

/*
 * draw_x_axis_label
 */
SiteHistogramGraph.prototype.draw_x_axis_label = function(ctx, w, h, label) {
  "use strict";
  ctx.save();
  ctx.font = "9pt Helvetica";
  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.fillText(label, w/2, h);
  ctx.restore();
};

/*
 * draw_y_axis
 *
 */
SiteHistogramGraph.prototype.draw_y_axis = function(ctx, w, h) {
  "use strict";
  var y_scale, p, y;
  ctx.save();
  ctx.font = "9pt Helvetica";
  ctx.textBaseline = "middle";
  ctx.textAlign = "right";

  y_scale = h / this.ymax;
  for (p = 0; p < this.ymax; p += this.yinc) {
    y = Math.round(h - p * y_scale) + 0.5;
    this.draw_y_tic(ctx, this.ydigits, y, p);
  }
  this.draw_y_tic(ctx, this.ydigits, 0.5, this.ymax);

  ctx.restore();
};

/*
 * draw_y_tic
 */
SiteHistogramGraph.prototype.draw_y_tic = function(ctx, digits, y, p) {
  "use strict";
  ctx.beginPath();
  ctx.moveTo(5, y);
  ctx.lineTo(-3, y);
  ctx.stroke();
  ctx.fillText(p.toFixed(digits), -5, y);
};

/*
 * draw_y_axis_label
 */
SiteHistogramGraph.prototype.draw_y_axis_label = function(ctx, w, h, label) {
  "use strict";
  ctx.save();
  ctx.translate(0, h/2);
  ctx.rotate(-Math.PI / 2);
  ctx.font = "9pt Helvetica";
  ctx.textBaseline = "top";
  ctx.textAlign = "center";
  ctx.fillText(label, 0, 0);
  ctx.restore();
};

/*
 * draw_bar
 */
function draw_bar(ctx, upperLeftCornerX, upperLeftCornerY, width, height, color){
    ctx.save();
    ctx.fillStyle=color;
    ctx.fillRect(upperLeftCornerX,upperLeftCornerY,width,height);
    ctx.restore();
}

/*
 * draw_bars
 */
SiteHistogramGraph.prototype.draw_bars = function(ctx, w, h, color) {
  "use strict";
  var i;
  var scale_y = h / this.ymax;
  var scale_x = w / (this.end - this.start);
  var numberOfBars = this.counts.length;
  var spacing = 5;
  var barSize = w/numberOfBars;
  var maxValue = 100;
  for (i=0; i<numberOfBars; i++) {
    var val = this.counts[i];
    var barHeight = Math.round( h * val/maxValue) ;
    draw_bar(
      ctx,
      (i - 0.5) * barSize + spacing/2,
      h - barHeight,
      barSize - spacing,
      barHeight,
      color
    );
  }
  ctx.restore();
};
