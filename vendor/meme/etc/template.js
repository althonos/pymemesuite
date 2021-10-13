//
// template.js
//

/*
 * Fill in a template variable
 */
function set_tvar(template, tvar, value) {
  var node;
  node = find_child(template, tvar);
  if (node === null) {
    throw new Error("Template does not contain variable " + tvar);
  }
  node.innerHTML = "";
  if (typeof value !== "object") {
    node.appendChild(document.createTextNode(value));
  } else {
    node.appendChild(value);
  }
} // set_tvar

/*
 * Get the text contained within the element.
 */
function elem_text(elem, separator) {
  if (separator === undefined) separator = '';
  return text_nodes(elem).map(node_text).join(separator);
}

/*
 * Get the text out of a specific text node.
 */
function node_text(node) {
  if (node === undefined) {
    return '';
  } else if (node.textContent) {
    return node.textContent;
  } else if (node.innerText) {
    return node.innerText;
  } else {
    return '';
  }
}

/*
 * Find all text nodes in the given container.
 */
function text_nodes(container) {
  var textNodes = [];
  var stack = [container];
  // depth first search to maintain ordering when flattened
  while (stack.length > 0) {
    var node = stack.pop();
    if (node.nodeType == Node.TEXT_NODE) {
      textNodes.push(node);
    } else {
      for (var i = node.childNodes.length-1; i >= 0; i--) {
        stack.push(node.childNodes[i]);
      }
    }
  }
  return textNodes;
}

/*
 * Create a button designed to contain a single symbol
 */
function make_sym_btn(symbol, title, action) {
  var box, sbox;
  box = document.createElement("div");
  box.tabIndex = 0;
  box.className = "sym_btn";
  sbox = document.createElement("span");
  if (typeof symbol === "string") {
    sbox.appendChild(document.createTextNode(symbol));
  } else {
    sbox.appendChild(symbol);
  }
  box.appendChild(sbox);
  box.title = title;
  box.addEventListener('click', action, false);
  box.addEventListener('keydown', action, false);
  return box;
}

/*
 * Create a pair of text spans with different classes.
 * This is useful when using CSS to only display one of them.
 */
function text_pair(txt1, cls1, txt2, cls2) {
  var container, part1, part2;
  container = document.createElement("span");
  part1 = document.createElement("span");
  part1.appendChild(document.createTextNode(txt1));
  part1.className = cls1;
  container.appendChild(part1);
  part2 = document.createElement("span");
  part2.appendChild(document.createTextNode(txt2));
  part2.className = cls2;
  container.appendChild(part2);
  return container;
}
