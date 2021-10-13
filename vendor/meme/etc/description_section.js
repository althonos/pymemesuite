//
// description_section.js
//

//
// Create a description element taking into account newlines (line breaks) 
// and multiple newlines (paragraph breaks) in the source text.
//
function make_description(container, description) {
  "use strict";
  var header, box, paragraphs, p, lines, i, j;
  container.innerHTML = "";
  if (description) {
    header = document.createElement("h2");
    header.className = "mainh pad2";
    header.id = "description";
    header.appendChild(document.createTextNode("Description"));
    container.appendChild(header);
    box = document.createElement("div");
    box.className = "box";
    container.appendChild(box);
    var text = description.replace(/\\n/g, "\n");
    paragraphs = text.split(/\n\n+/);
    for (i = 0; i < paragraphs.length; i++) {
      p = document.createElement("p");
      lines = paragraphs[i].split("\n");
      for (j = 0; j < lines.length; j++) {
        if (j != 0) p.appendChild(document.createElement('br'));
        p.appendChild(document.createTextNode(lines[j]));
      }
      box.appendChild(p);
    }
  }
}
