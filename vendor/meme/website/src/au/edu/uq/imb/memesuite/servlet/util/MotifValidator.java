package au.edu.uq.imb.memesuite.servlet.util;

import java.io.*;
import java.math.BigDecimal;
import java.util.*;
import java.util.regex.*;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;

import au.edu.uq.imb.memesuite.data.Alph;
import au.edu.uq.imb.memesuite.data.Alph.Symbol;
import au.edu.uq.imb.memesuite.data.AlphStd;
import au.edu.uq.imb.memesuite.data.MotifStats;
import au.edu.uq.imb.memesuite.io.alph.AlphParser;
import au.edu.uq.imb.memesuite.io.html.HDataEventType;
import au.edu.uq.imb.memesuite.io.html.HDataParseException;
import au.edu.uq.imb.memesuite.io.html.HDataPullParser;

import static au.edu.uq.imb.memesuite.io.html.HDataPullParser.getEntry;
import static au.edu.uq.imb.memesuite.io.html.HDataPullParser.getOptEntry;
import static au.edu.uq.imb.memesuite.io.html.HDataPullParser.hasEntry;
import static javax.xml.stream.XMLStreamConstants.END_ELEMENT;
import static javax.xml.stream.XMLStreamConstants.START_ELEMENT;

import static au.edu.uq.imb.memesuite.io.html.HDataEventType.*;

public class MotifValidator {

  public static class Attr {
    public final String name;
    protected boolean required;
    protected String value;
    public Attr(String attributeName) {
      name = attributeName;
      value = null;
      required = true;
    }
    public Attr(String attributeName, boolean optional) {
      name = attributeName;
      value = null;
      required = !optional;
    }
    public boolean isRequired() {
      return required;
    }
    public boolean store(String value) {
      this.value = value;
      return true;
    }
    public String getValue() {
      return value;
    }
  }
  public static class AttrRegex extends Attr {
    protected  Pattern regex;
    public AttrRegex(String attributeName, Pattern attributeRegex, boolean optional) {
      super(attributeName, optional);
      regex = attributeRegex;
    }
    public AttrRegex(String attributeName, Pattern attributeRegex) {
      this(attributeName, attributeRegex, false);
    }
    public AttrRegex(String attributeName, String attributeRegex, boolean optional) {
      this(attributeName, Pattern.compile(attributeRegex), optional);
    }
    public AttrRegex(String attributeName, String attributeRegex) {
      this(attributeName, Pattern.compile(attributeRegex), false);
    }
    @Override
    public boolean store(String value) {
      Matcher m;
      if ((m = regex.matcher(value)).matches()) {
        this.value = (m.groupCount() >= 1 ? m.group(1) : value);
        return true;
      }
      return false;
    }
  }
  public static class MotifParseException extends Exception {
    public MotifParseException(String message) {
      super(message);
    }

    public MotifParseException(String message, Throwable cause) {
      super(message, cause);
    }

    public MotifParseException(Throwable cause) {
      super(cause);
    }
  }

  protected MotifValidator() {

  }

  private static void nextTag(XMLStreamReader cursor, FeedbackHandler feedback) throws XMLStreamException{
    //boolean trace = true;		// Set true for debugging parsing problems
    boolean trace = false;		// Set false for normal
    cursor.nextTag();
    if (trace && feedback != null) feedback.whine(cursor.getLocalName());
  }

  private static void tagAtrs(XMLStreamReader cursor, boolean disallowUnknown, FeedbackHandler feedback, Attr... attributes) throws XMLStreamException{
    final int len = cursor.getAttributeCount();
    Map<String,Attr> attributesLookup = new TreeMap<String, Attr>();
    for (Attr attribute : attributes) attributesLookup.put(attribute.name, attribute);
    for (int i = 0; i < len; i++) {
      String name = cursor.getAttributeLocalName(i);
      String value = cursor.getAttributeValue(i);
      Attr attr = attributesLookup.get(name);
      if (attr == null) {
        //if (disallowUnknown) throw new XMLStreamException("Unknown attribute " + name, cursor.getLocation());
        if (disallowUnknown && feedback != null) feedback.whine("Unknown attribute " + name + " " + cursor.getLocation());
        continue;
      }
      if (!attr.store(value)) {
        //throw new XMLStreamException("Attribute " + name + " has invalid value", cursor.getLocation());
        if (feedback != null) feedback.whine("Attribute " + name + " has invalid value " + cursor.getLocation());
      }
    }
    for (Attr attribute : attributes) {
      if (attribute.isRequired() && attribute.getValue() == null) {
        //throw new XMLStreamException("Missing attribute " + attribute.name, cursor.getLocation());
        if (feedback != null) feedback.whine("Missing attribute " + attribute.name + " " + cursor.getLocation());
      }
    }
  }

  private static BigDecimal[] parseAlphabetArrayMemeXML(XMLStreamReader cursor, FeedbackHandler feedback, Map<String,String> id2sym, Alph alph) throws XMLStreamException {
    BigDecimal[] values = new BigDecimal[alph.sizeCore()];
    Arrays.fill(values, null);
    cursor.require(START_ELEMENT, null, "alphabet_array");
    nextTag(cursor, feedback);
    for (int i = 0; i < alph.sizeCore(); i++) {
      cursor.require(START_ELEMENT, null, "value");
      Attr letterId = new Attr("letter_id");
      tagAtrs(cursor, true, null, letterId);
      String sym = id2sym.get(letterId.getValue());
      if (sym == null) throw new XMLStreamException("Unknown letter_id!", cursor.getLocation());
      Symbol symbol = alph.find(sym);
      if (!symbol.isCore()) throw new XMLStreamException("letter_id does not refer to a core symbol!", cursor.getLocation());
      if (values[symbol.index] != null) throw new XMLStreamException("The symbol " + symbol.sym + " already has a value assigned!", cursor.getLocation());
      String letterValue = cursor.getElementText();
      try {
        values[symbol.index] = new BigDecimal(letterValue.trim());
      } catch (NumberFormatException e) {
        throw new XMLStreamException("Value is not a number!", e);
      }
      cursor.require(END_ELEMENT, null, "value");
      nextTag(cursor, feedback);
    }
    cursor.require(END_ELEMENT, null, "alphabet_array");
    return values;
  }

  private static BigDecimal[][] parseAlphabetMatrixMemeXML(XMLStreamReader cursor, FeedbackHandler feedback, Map<String,String> id2sym, Alph alph, int length) throws XMLStreamException {
    List<BigDecimal[]> values = new ArrayList<BigDecimal[]>();
    cursor.require(START_ELEMENT, null, "alphabet_matrix");
    nextTag(cursor, feedback);
    for (int i = 0; i < length; i++) {
      cursor.require(START_ELEMENT, null, "alphabet_array");
      values.add(parseAlphabetArrayMemeXML(cursor, feedback, id2sym, alph));
      cursor.require(END_ELEMENT, null, "alphabet_array");
      nextTag(cursor, feedback);
    }
    cursor.require(END_ELEMENT, null, "alphabet_matrix");
    return values.toArray(new BigDecimal[values.size()][]);
  }

  private static void parseMemeXML(MotifStats stats, XMLStreamReader cursor, FeedbackHandler feedback) throws XMLStreamException, AlphParser.AlphParseException {
    Map<String,String> id2sym = new TreeMap<String, String>();
    cursor.require(START_ELEMENT, null, "MEME");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "training_set");
    tagAtrs(cursor, true, feedback, 
      // post 5.0 format
      new Attr("primary_sequences", true), new Attr("primary_count", true), new Attr("primary_positions", true),
      new Attr("control_sequences", true), new Attr("control_count", true), new Attr("control_positions", true),
      // pre 5.0 format
      new Attr("datafile", true), new Attr("length", true));
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "alphabet");
    Attr alphabetId = new AttrRegex("id", "^(amino-acid|nucleotide)$", true),
        alphabetLength = new AttrRegex("length", "^(\\d+)$", true),
        alphabetName = new Attr("name", true),
        alphabetLike = new AttrRegex("like", "^(?:rna|dna|protein)$", true);
    tagAtrs(cursor, true, null, alphabetId, alphabetLength, alphabetName, alphabetLike);
    Alph alph;
    AlphParser parser = null;
    if (alphabetId.getValue() != null) {
      if (alphabetId.getValue().equals("amino-acid")) {
        alph = AlphStd.PROTEIN.getAlph();
      } else {
        alph = AlphStd.DNA.getAlph();
      }
    } else {
      alph = null;
      parser = new AlphParser();
      parser.parseHeader(alphabetName.getValue(), (alphabetLike.getValue() != null ? alphabetLike.getValue().toUpperCase() : null));
    }
    nextTag(cursor, feedback);
    while (!(cursor.isEndElement() && cursor.getLocalName().equals("alphabet"))) {
      cursor.require(START_ELEMENT, null, "letter");
      Attr letterAliases = new Attr("aliases", true), letterColour = new Attr("colour", true),
          letterComplement = new Attr("complement", true), letterEquals = new Attr("equals", true),
          letterId = new Attr("id"), letterName = new Attr("name", true), letterSymbol = new AttrRegex("symbol", "^(.)$");
      tagAtrs(cursor, true, null, letterAliases, letterColour, letterComplement, letterEquals, letterId, letterName, letterSymbol);
      nextTag(cursor, feedback);
      id2sym.put(letterId.getValue(), letterSymbol.getValue());
      if (parser != null) {
        parser.parseSymbol(letterSymbol.getValue(), letterName.getValue(), letterColour.getValue(),
            letterComplement.getValue(), letterEquals.getValue(), letterAliases.getValue());
      }
      cursor.require(END_ELEMENT, null, "letter");
      nextTag(cursor, feedback);
    }
    if (parser != null) {
      alph = parser.parseDone();
    }
    stats.setAlphabet(alph);
    cursor.require(END_ELEMENT, null, "alphabet");
    nextTag(cursor, feedback);
    if (alphabetId.getValue() != null) {
      cursor.require(START_ELEMENT, null, "ambigs");
      nextTag(cursor, feedback);
      while (!(cursor.isEndElement() && cursor.getLocalName().equals("ambigs"))) {
        cursor.require(START_ELEMENT, null, "letter");
        tagAtrs(cursor, true, null, new Attr("id"), new Attr("symbol"));
        nextTag(cursor, feedback);
        cursor.require(END_ELEMENT, null, "letter");
        nextTag(cursor, feedback);
      }
      cursor.require(END_ELEMENT, null, "ambigs");
      nextTag(cursor, feedback);
    }
    while (cursor.isStartElement() && cursor.getLocalName().equals("sequence")) {
    //do {
      cursor.require(START_ELEMENT, null, "sequence");
      tagAtrs(cursor, true, feedback, new Attr("id"), new Attr("length"), new Attr("name"), new Attr("weight"));
      nextTag(cursor, feedback);
      cursor.require(END_ELEMENT, null, "sequence");
      nextTag(cursor, feedback);
    //} while (cursor.isStartElement() && cursor.getLocalName().equals("sequence"));
    }
    cursor.require(START_ELEMENT, null, "letter_frequencies");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "alphabet_array");
    parseAlphabetArrayMemeXML(cursor, feedback, id2sym, alph);
    cursor.require(END_ELEMENT, null, "alphabet_array");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "letter_frequencies");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "training_set");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "command_line");
    for (cursor.next(); !(cursor.isEndElement() && cursor.getLocalName().equals("command_line")); cursor.next()) {
      if (!(cursor.isStartElement() && cursor.getLocalName().equals("arg") ||
          cursor.isEndElement() && cursor.getLocalName().equals("arg") || cursor.isCharacters())) {
        throw new XMLStreamException("Expected content of command_line", cursor.getLocation());
      }
    }
    cursor.require(END_ELEMENT, null, "command_line");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "host");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "host");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "type");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "type");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "nmotifs");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "nmotifs");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "evalue_threshold");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "evalue_threshold");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "object_function");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "object_function");
    nextTag(cursor, feedback);
    try {
      // post 5.0
      cursor.require(START_ELEMENT, null, "spfun");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "spfun");
    } catch (XMLStreamException e) {
      // pre 5.0
      cursor.require(START_ELEMENT, null, "use_llr");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "use_llr");
    }
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "min_width");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "min_width");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "max_width");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "max_width");
    nextTag(cursor, feedback);
    if (cursor.isStartElement() && cursor.getLocalName().equals("minic")) {
      cursor.require(START_ELEMENT, null, "minic");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "minic");
      nextTag(cursor, feedback);
    }

    if (cursor.isStartElement() && cursor.getLocalName().equals("wg")) {
      cursor.require(START_ELEMENT, null, "wg");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "wg");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "ws");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "ws");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "endgaps");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "endgaps");
      nextTag(cursor, feedback);
    }
    if (cursor.isStartElement() && cursor.getLocalName().equals("substring")) {
      cursor.require(START_ELEMENT, null, "substring");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "substring");
      nextTag(cursor, feedback);
    }
    cursor.require(START_ELEMENT, null, "minsites");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "minsites");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "maxsites");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "maxsites");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "wnsites");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "wnsites");
    nextTag(cursor, feedback);
    if (cursor.isStartElement() && cursor.getLocalName().equals("prob")) {
      cursor.require(START_ELEMENT, null, "prob");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "prob");
      nextTag(cursor, feedback);
    }
    cursor.require(START_ELEMENT, null, "spmap");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "spmap");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "spfuzz");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "spfuzz");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "prior");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "prior");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "beta");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "beta");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "maxiter");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "maxiter");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "distance");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "distance");
    nextTag(cursor, feedback);
    try {
      // pre 5.0
      cursor.require(START_ELEMENT, null, "num_sequences");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "num_sequences");
      nextTag(cursor, feedback);
    } catch (XMLStreamException e) {
      // post 5.0
    }
    cursor.require(START_ELEMENT, null, "num_positions");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "num_positions");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "seed");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "seed");
    nextTag(cursor, feedback);
    try {
      // post 5.0
      cursor.require(START_ELEMENT, null, "hsfrac");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "hsfrac");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "searchsize");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "searchsize");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "maxsize");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "maxsize");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "norand");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "norand");
      nextTag(cursor, feedback);
      if (cursor.isStartElement() && cursor.getLocalName().equals("csites")) {
	cursor.require(START_ELEMENT, null, "csites");
	cursor.getElementText();
	cursor.require(END_ELEMENT, null, "csites");
        nextTag(cursor, feedback);
      }
    } catch (XMLStreamException e) {
      // pre 5.0
      cursor.require(START_ELEMENT, null, "ctfrac");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "ctfrac");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "maxwords");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "maxwords");
      nextTag(cursor, feedback);
    }
    cursor.require(START_ELEMENT, null, "strands");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "strands");
    nextTag(cursor, feedback);
    try {
      // post 5.0
      cursor.require(START_ELEMENT, null, "brief");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "brief");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "psp_file");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "psp_file");
      nextTag(cursor, feedback);
    } catch (XMLStreamException e) {
      // pre 5.0
    }
    cursor.require(START_ELEMENT, null, "priors_file");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "priors_file");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "reason_for_stopping");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "reason_for_stopping");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "background_frequencies");
    tagAtrs(cursor, true, null, new Attr("source"), new Attr("order"));
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "alphabet_array");
    parseAlphabetArrayMemeXML(cursor, feedback, id2sym, alph);
    cursor.require(END_ELEMENT, null, "alphabet_array");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "background_frequencies");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "motifs");
    nextTag(cursor, feedback);
    while (!(cursor.isEndElement() && cursor.getLocalName().equals("motifs"))) {
      cursor.require(START_ELEMENT, null, "motif");
      Attr aMotifWidth = new AttrRegex("width", "^(\\d+)$");
      tagAtrs(cursor, true, feedback, new Attr("id"), new Attr("name"), new Attr("alt", true), aMotifWidth, new Attr("sites"), new Attr("ic"), new Attr("re"), new Attr("llr"), new Attr("p_value", true), new Attr("e_value"), new Attr("bayes_threshold"), new Attr("elapsed_time"), new Attr("url", true), new Attr("total_sites", true), new Attr("site_distr", true));
      int motifWidth;
      try {
        motifWidth = Integer.parseInt(aMotifWidth.getValue());
      } catch (NumberFormatException e) {
        throw new XMLStreamException("Failed to parse motif width as a number.", cursor.getLocation(), e);
      }
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "scores");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "alphabet_matrix");
      parseAlphabetMatrixMemeXML(cursor, feedback, id2sym, alph, motifWidth);
      cursor.require(END_ELEMENT, null, "alphabet_matrix");
      nextTag(cursor, feedback);
      cursor.require(END_ELEMENT, null, "scores");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "probabilities");
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "alphabet_matrix");
      parseAlphabetMatrixMemeXML(cursor, feedback, id2sym, alph, motifWidth);
      cursor.require(END_ELEMENT, null, "alphabet_matrix");
      nextTag(cursor, feedback);
      cursor.require(END_ELEMENT, null, "probabilities");
      nextTag(cursor, feedback);
      if (cursor.isStartElement() && cursor.getLocalName().equals("regular_expression")) {
        cursor.require(START_ELEMENT, null, "regular_expression");
        cursor.getElementText();
        cursor.require(END_ELEMENT, null, "regular_expression");
        nextTag(cursor, feedback);
      }
      cursor.require(START_ELEMENT, null, "contributing_sites");
      nextTag(cursor, feedback);
      //do {
      while (cursor.isStartElement() && cursor.getLocalName().equals("contributing_site")) {
        cursor.require(START_ELEMENT, null, "contributing_site");
        tagAtrs(cursor, true, null, new Attr("position"), new Attr("pvalue"), new Attr("sequence_id"), new Attr("strand", true));
        nextTag(cursor, feedback);
        cursor.require(START_ELEMENT, null, "left_flank");
        cursor.getElementText();
        cursor.require(END_ELEMENT, null, "left_flank");
        nextTag(cursor, feedback);
        cursor.require(START_ELEMENT, null, "site");
        nextTag(cursor, feedback);
        do {
          cursor.require(START_ELEMENT, null, "letter_ref");
          tagAtrs(cursor, true, null, new Attr("letter_id"));
          nextTag(cursor, feedback);
          cursor.require(END_ELEMENT, null, "letter_ref");
          nextTag(cursor, feedback);
        } while (!(cursor.isEndElement() && cursor.getLocalName().equals("site")));
        cursor.require(END_ELEMENT, null, "site");
        nextTag(cursor, feedback);
        cursor.require(START_ELEMENT, null, "right_flank");
        cursor.getElementText();
        cursor.require(END_ELEMENT, null, "right_flank");
        nextTag(cursor, feedback);
        cursor.require(END_ELEMENT, null, "contributing_site");
        nextTag(cursor, feedback);
      //} while (!(cursor.isEndElement() && cursor.getLocalName().equals("contributing_sites")));
      }
      cursor.require(END_ELEMENT, null, "contributing_sites");
      nextTag(cursor, feedback);
      cursor.require(END_ELEMENT, null, "motif");
      nextTag(cursor, feedback);
      stats.addMotif(motifWidth);
    }
    cursor.require(END_ELEMENT, null, "motifs");
    nextTag(cursor, feedback);
    if (cursor.isStartElement() && cursor.getLocalName().equals("scanned_sites_summary")) {
      cursor.require(START_ELEMENT, null, "scanned_sites_summary");
      tagAtrs(cursor, true, null, new Attr("p_thresh"));
      nextTag(cursor, feedback);
      do {
        cursor.require(START_ELEMENT, null, "scanned_sites");
        tagAtrs(cursor, true, null, new Attr("num_sites"), new Attr("pvalue"), new Attr("sequence_id"));
        nextTag(cursor, feedback);
        while (!(cursor.isEndElement() && cursor.getLocalName().equals("scanned_sites"))) {
          cursor.require(START_ELEMENT, null, "scanned_site");
          tagAtrs(cursor, true, null, new Attr("motif_id"), new Attr("position"), new Attr("pvalue"),
              new AttrRegex("strand", "^(minus|none|plus)$"));
          nextTag(cursor, feedback);
          cursor.require(END_ELEMENT, null, "scanned_site");
          nextTag(cursor, feedback);
        }
        cursor.require(END_ELEMENT, null, "scanned_sites");
        nextTag(cursor, feedback);
      } while (!(cursor.isEndElement() && cursor.getLocalName().equals("scanned_sites_summary")));
      cursor.require(END_ELEMENT, null, "scanned_sites_summary");
      nextTag(cursor, feedback);
    }
    //if (feedback != null) feedback.whine("at end of MEME")
    cursor.require(END_ELEMENT, null, "MEME");
  }

  private static void parseDremeXML(MotifStats stats, XMLStreamReader cursor, FeedbackHandler feedback) throws XMLStreamException, AlphParser.AlphParseException {
    cursor.require(START_ELEMENT, null, "dreme");
    tagAtrs(cursor, true, null, new AttrRegex("version", "(\\d+\\.\\d+\\.\\d+)"), new Attr("release"));
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "command_line");
    for (cursor.next(); !(cursor.isEndElement() && cursor.getLocalName().equals("command_line")); cursor.next()) {
      if (!(cursor.isStartElement() && cursor.getLocalName().equals("arg") || cursor.isEndElement() && cursor.getLocalName().equals("arg") || cursor.isCharacters())) {
        throw new XMLStreamException("Expected content of command_line", cursor.getLocation());
      }
    }
    cursor.require(END_ELEMENT, null, "command_line");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "positives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^(\\d+)$"), new Attr("file"),
        new Attr("last_mod_date"), new Attr("name"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "positives");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "negatives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^\\d+$"), new Attr("file", true),
        new AttrRegex("from", "^(?:file|shuffled)$"), new Attr("last_mod_date", true), new Attr("name"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "negatives");
    Alph alph;
    nextTag(cursor, feedback);
    if (cursor.getLocalName().equals("alphabet")) {
      AlphParser alphParser = new AlphParser();
      cursor.require(START_ELEMENT, null, "alphabet");
      Attr name = new Attr("name", true);
      Attr like = new AttrRegex("like", "^(?:rna|dna|protein)$", true);
      tagAtrs(cursor, true, null, name, like);
      alphParser.parseHeader(name.getValue(), (like.getValue() != null ? like.getValue().toUpperCase() : null));
      nextTag(cursor, feedback);
      while (cursor.getLocalName().equals("letter")) {
        cursor.require(START_ELEMENT, null, "letter");
        Attr aAliases = new Attr("aliases", true), aColour = new Attr("colour", true),
            aComplement = new Attr("complement", true), aEquals = new Attr("equals", true), aId = new Attr("id"),
            aName = new Attr("name", true), aSymbol = new AttrRegex("symbol", "^(.)$");
        tagAtrs(cursor, true, null, aAliases, aColour, aComplement, aEquals, aId, aName, aSymbol);
        alphParser.parseSymbol(aSymbol.getValue(), aName.getValue(), aColour.getValue(), aComplement.getValue(),
            aEquals.getValue(), aAliases.getValue());
        nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "letter");
        nextTag(cursor, feedback);
      }
      cursor.require(END_ELEMENT, null, "alphabet");
      nextTag(cursor, feedback);
      alph = alphParser.parseDone();
      cursor.require(START_ELEMENT, null, "strands");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "strands");
      nextTag(cursor, feedback);
    } else {
      alph = AlphStd.DNA.getAlph(); // if the alphabet is not specified then assume DNA
    }
    stats.setAlphabet(alph);
    cursor.require(START_ELEMENT, null, "background");
    nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "background");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "stop");
    tagAtrs(cursor, true, null, new Attr("count", true), new Attr("evalue", true), new Attr("time", true));
    nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "stop");
    nextTag(cursor, feedback);
    if (cursor.getLocalName().equals("norc")) {
      cursor.require(START_ELEMENT, null, "norc");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "norc");
      nextTag(cursor, feedback);
    }
    cursor.require(START_ELEMENT, null, "ngen");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "ngen");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "add_pv_thresh");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "add_pv_thresh");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "seed");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "seed");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "host");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "host");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "when");
    cursor.getElementText();
    cursor.require(END_ELEMENT, null, "when");
    nextTag(cursor, feedback);
    if (cursor.isStartElement() && cursor.getLocalName().equals("description")) {
      cursor.require(START_ELEMENT, null, "description");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "description");
      nextTag(cursor, feedback);
    }
    cursor.require(END_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "motifs");
    nextTag(cursor, feedback);
    while (!(cursor.isEndElement() && cursor.getLocalName().equals("motifs"))) {
      cursor.require(START_ELEMENT, null, "motif");
      Attr aLength = new AttrRegex("length", "^(\\d+)$");
      tagAtrs(cursor, true, null, new Attr("evalue"), new Attr("id"), aLength, new Attr("n"),
          new Attr("nsites"), new Attr("p"), new Attr("pvalue"), new Attr("seq"), new Attr("unerased_evalue"));
      int motifLength;
      try {
        motifLength = Integer.parseInt(aLength.getValue(), 10);
      } catch (NumberFormatException e) {
        throw new XMLStreamException("Failed to parse length attribute as an integer", cursor.getLocation(), e);
      }
      nextTag(cursor, feedback);
      for (int i = 0; i < motifLength; i++) {
        cursor.require(START_ELEMENT, null, "pos");
        nextTag(cursor, feedback);
        cursor.require(END_ELEMENT, null, "pos");
        nextTag(cursor, feedback);
      }
      while (!(cursor.isEndElement() && cursor.getLocalName().equals("motif"))) {
        cursor.require(START_ELEMENT, null, "match");
        tagAtrs(cursor, true, null, new Attr("evalue"), new Attr("n"), new Attr("p"), new Attr("pvalue"), new Attr("seq"));
        nextTag(cursor, feedback);
        cursor.require(END_ELEMENT, null, "match");
        nextTag(cursor, feedback);
      }
      cursor.require(END_ELEMENT, null, "motif");
      nextTag(cursor, feedback);
      stats.addMotif(motifLength);
    }
    cursor.require(END_ELEMENT, null, "motifs");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "run_time");
    tagAtrs(cursor, true, null, new Attr("cpu"), new Attr("real"), new AttrRegex("stop", "^(count|evalue|time)$"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "run_time");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "dreme");
  }

  private static void parseStremeXML(MotifStats stats, XMLStreamReader cursor, FeedbackHandler feedback) throws XMLStreamException, AlphParser.AlphParseException {
    Map<String,String> id2sym = new TreeMap<String, String>();
    cursor.require(START_ELEMENT, null, "STREME");
    tagAtrs(cursor, true, null, new AttrRegex("version", "(\\d+\\.\\d+\\.\\d+)"), new Attr("release"));
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "command_line");
    for (cursor.next(); !(cursor.isEndElement() && cursor.getLocalName().equals("command_line")); cursor.next()) {
      if (!(cursor.isStartElement() && cursor.getLocalName().equals("arg") || cursor.isEndElement() && cursor.getLocalName().equals("arg") || cursor.isCharacters())) {
        throw new XMLStreamException("Expected content of command_line", cursor.getLocation());
      }
    }
    cursor.require(END_ELEMENT, null, "command_line");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "train_positives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^(\\d+)$"), new AttrRegex("positions", "^(\\d+)$"), new Attr("file"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "train_positives");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "train_negatives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^(\\d+)$"), new AttrRegex("positions", "^(\\d+)$"),
      new AttrRegex("from", "^(?:file|shuffled)$"), new Attr("file", true));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "train_negatives");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "test_positives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^(\\d+)$"), new AttrRegex("positions", "^(\\d+)$"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "test_positives");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "test_negatives");
    tagAtrs(cursor, true, null, new AttrRegex("count", "^(\\d+)$"), new AttrRegex("positions", "^(\\d+)$"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "test_negatives");
    nextTag(cursor, feedback);
    if (cursor.getLocalName().equals("sequence_db")) {		// Moved down after 5.3.0
      cursor.require(START_ELEMENT, null, "sequence_db"); nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "sequence_db");
      nextTag(cursor, feedback);
    }
    Alph alph;
    if (cursor.getLocalName().equals("alphabet")) {
      AlphParser alphParser = new AlphParser();
      cursor.require(START_ELEMENT, null, "alphabet");
      Attr name = new Attr("name", true);
      Attr like = new AttrRegex("like", "^(?:rna|dna|protein)$", true);
      tagAtrs(cursor, true, null, name, like);
      alphParser.parseHeader(name.getValue(), (like.getValue() != null ? like.getValue().toUpperCase() : null));
      nextTag(cursor, feedback);
      while (cursor.getLocalName().equals("letter")) {
        cursor.require(START_ELEMENT, null, "letter");
        Attr aAliases = new Attr("aliases", true), aColour = new Attr("colour", true),
	  aComplement = new Attr("complement", true), aEquals = new Attr("equals", true), aId = new Attr("id"),
	  aName = new Attr("name", true), aSymbol = new AttrRegex("symbol", "^(.)$");
        tagAtrs(cursor, true, null, aAliases, aColour, aComplement, aEquals, aId, aName, aSymbol);
        id2sym.put(aId.getValue(), aSymbol.getValue());
        alphParser.parseSymbol(aSymbol.getValue(), aName.getValue(), aColour.getValue(), aComplement.getValue(),
            aEquals.getValue(), aAliases.getValue());
        nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "letter");
        nextTag(cursor, feedback);
      }
      cursor.require(END_ELEMENT, null, "alphabet");
      nextTag(cursor, feedback);
      alph = alphParser.parseDone();
      cursor.require(START_ELEMENT, null, "strands");
      cursor.getElementText();
      cursor.require(END_ELEMENT, null, "strands");
      nextTag(cursor, feedback);
    } else {
      alph = AlphStd.DNA.getAlph(); // if the alphabet is not specified then assume DNA
    }
    stats.setAlphabet(alph);
    if (cursor.getLocalName().equals("sequence_db")) {		// Moved here after 5.3.0
      cursor.require(START_ELEMENT, null, "sequence_db");
      nextTag(cursor, feedback);
      nextTag(cursor, feedback);
    }
    if (cursor.getLocalName().equals("background")) {		// 5.3.0 
      cursor.require(START_ELEMENT, null, "background");
      nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "background");
    } else {
      cursor.require(START_ELEMENT, null, "background_frequencies");
      tagAtrs(cursor, true, null, new Attr("source"), new Attr("order"));
      nextTag(cursor, feedback);
      cursor.require(START_ELEMENT, null, "alphabet_array");
      parseAlphabetArrayMemeXML(cursor, feedback, id2sym, alph);
      cursor.require(END_ELEMENT, null, "alphabet_array");
      nextTag(cursor, feedback);
      cursor.require(END_ELEMENT, null, "background_frequencies");
    }
    nextTag(cursor, feedback); 
    cursor.require(START_ELEMENT, null, "stop");
    tagAtrs(cursor, true, null, new Attr("count", true), new Attr("evalue", true), new Attr("time", true));
    nextTag(cursor, feedback); cursor.require(END_ELEMENT, null, "stop");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "objfun"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "objfun");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "test"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "test");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "minw"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "minw");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "maxw"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "maxw");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "kmer"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "kmer");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "hofract"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "hofract");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "neval"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "neval");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "nref"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "nref");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "niter"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "niter");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "patience"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "patience");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "seed"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "seed");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "useer"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "useer");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "minscore"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "minscore");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "ignore_depth"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "ignore_depth");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "nsubsets"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "nsubsets");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "min_pal_ratio"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "min_pal_ratio");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "max_pal_ed"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "max_pal_ed");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "cand"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "cand");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "experimental"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "experimental");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "totallength"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "totallength");
    nextTag(cursor, feedback); 
    if (cursor.getLocalName().equals("align")) {		// added in 5.4.0 
      cursor.require(START_ELEMENT, null, "align"); 
      cursor.getElementText(); 
      cursor.require(END_ELEMENT, null, "align");
    }
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "host"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "host");
    nextTag(cursor, feedback); 
    if (cursor.getLocalName().equals("description")) {		// optional
      cursor.require(START_ELEMENT, null, "description"); 
      cursor.getElementText(); 
      cursor.require(END_ELEMENT, null, "description");
      nextTag(cursor, feedback);
    }
    cursor.require(END_ELEMENT, null, "model");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "motifs");
    nextTag(cursor, feedback);
    while (!(cursor.isEndElement() && cursor.getLocalName().equals("motifs"))) {
      cursor.require(START_ELEMENT, null, "motif");
      Attr aLength = new AttrRegex("width", "^(\\d+)$");
      tagAtrs(cursor, true, null, new Attr("id"), aLength, new Attr("train_pos_count"),
        new Attr("train_neg_count"), new Attr("test_pvalue"));
      int motifLength;
      try {
        motifLength = Integer.parseInt(aLength.getValue(), 10);
      } catch (NumberFormatException e) {
        throw new XMLStreamException("Failed to parse length attribute as an integer", cursor.getLocation(), e);
      }
      nextTag(cursor, feedback);
      for (int i = 0; i < motifLength; i++) {
        cursor.require(START_ELEMENT, null, "pos");
        nextTag(cursor, feedback);
        cursor.require(END_ELEMENT, null, "pos");
        nextTag(cursor, feedback);
      }
      cursor.require(END_ELEMENT, null, "motif");
      nextTag(cursor, feedback);
      stats.addMotif(motifLength);
    }
    cursor.require(END_ELEMENT, null, "motifs");
    nextTag(cursor, feedback); cursor.require(START_ELEMENT, null, "reason_for_stopping"); cursor.getElementText(); cursor.require(END_ELEMENT, null, "reason_for_stopping");
    nextTag(cursor, feedback);
    cursor.require(START_ELEMENT, null, "run_time");
    tagAtrs(cursor, true, null, new Attr("cpu"));
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "run_time");
    nextTag(cursor, feedback);
    cursor.require(END_ELEMENT, null, "STREME");
  }

  public static MotifStats validateAsXML(InputStream in, FeedbackHandler feedback) throws IOException, MotifParseException {
    MotifStats stats = new MotifStats();
    try {
      XMLInputFactory xmlif = XMLInputFactory.newFactory();
      XMLStreamReader cursor = xmlif.createXMLStreamReader(in,"UTF-8");
      while(!cursor.isStartElement()) cursor.next();
      if (cursor.getLocalName().equals("MEME")) {
        parseMemeXML(stats, cursor, feedback);
      } else if (cursor.getLocalName().equals("STREME")) {
        parseStremeXML(stats, cursor, feedback);
      } else if (cursor.getLocalName().equals("dreme")) {
        parseDremeXML(stats, cursor, feedback);
      } else {
        throw new MotifParseException("Expected a MEME, STREME or DREME XML file.");
      }
    } catch (XMLStreamException e) {
      Throwable nestErr = e.getNestedException();
      if (nestErr != null && nestErr instanceof IOException) {
        throw (IOException)nestErr;
      }
      throw new MotifParseException(e);
    } catch (AlphParser.AlphParseException e) {
      throw new MotifParseException(e);
    }
    stats.freeze();
    return stats;
  }

  public static MotifStats parseJson(HDataPullParser cursor, FeedbackHandler feedback) throws IOException, AlphParser.AlphParseException, HDataParseException {
    MotifStats stats = new MotifStats();
    cursor.next(JSON_START_DATA);
    cursor.next(JSON_START_OBJECT);
    String program = cursor.nextStringProperty("program");
    boolean isMeme = program.equalsIgnoreCase("meme");
    boolean isDreme = program.equalsIgnoreCase("dreme");
    boolean isStreme = program.equalsIgnoreCase("streme");
    cursor.nextStringProperty("version");
    cursor.nextStringProperty("release");
    if (cursor.isProperty("stop_reason")) {
      cursor.nextStringProperty("stop_reason");
    }
    if (cursor.isProperty("description")) {
      cursor.nextStringProperty("description");
    }
    cursor.nextListProperty("cmd", String.class);
    cursor.skipProperty("options"); // skip over options
    Object alphabetObj = cursor.nextProperty("alphabet");
    if (hasEntry(alphabetObj, "symbols", String.class)) {
      stats.setAlphabet(AlphStd.fromCore((String) getEntry(alphabetObj, "symbols", String.class)));
    } else {
      AlphParser alphParser = new AlphParser();
      String alphName = (String) getEntry(alphabetObj, "name", String.class);
      String alphLike = (String)getOptEntry(alphabetObj, null, "like", String.class);
      alphParser.parseHeader(alphName, (alphLike != null ? alphLike.toUpperCase() : null));
      List symbolsObj = (List) getEntry(alphabetObj, "symbols", List.class);
      for (Object symbolObj : symbolsObj) {
        String sym = (String) getEntry(symbolObj, "symbol", String.class);
        String name = (String) getOptEntry(symbolObj, null, "name", String.class);
        String colour = (String) getOptEntry(symbolObj, null, "colour", String.class);
        String complement = (String) getOptEntry(symbolObj, null, "complement", String.class);
        String equals = (String) getOptEntry(symbolObj, null, "equals", String.class);
        String aliases = (String) getOptEntry(symbolObj, null, "aliases", String.class);
        alphParser.parseSymbol(sym, name, colour, complement, equals, aliases);
      }
      stats.setAlphabet(alphParser.parseDone());
    }
    if (cursor.isProperty("background")) { 		// start background
      cursor.next();
      cursor.next(JSON_START_OBJECT);
      if (cursor.isProperty("source")) cursor.nextStringProperty("source");
      if (cursor.isProperty("order")) cursor.nextNumberProperty("order");
      cursor.nextListProperty("freqs", Number.class);
      cursor.next(JSON_END_OBJECT);
      cursor.next(JSON_END_PROPERTY); 			// end background
    }
    if (isMeme) {
      cursor.next("sequence_db"); 			// start sequence_db
      cursor.next(JSON_START_OBJECT);
      if (cursor.isProperty("source")) {		// pre 5.0
        cursor.nextStringProperty("source");
      } else {						// post 5.0
        cursor.nextStringProperty("primary_source");
        cursor.nextNumberProperty("primary_count");
        cursor.nextNumberProperty("primary_positions");
        if (cursor.isProperty("control_source")) {	// optional
        cursor.nextStringProperty("control_source");
        cursor.nextNumberProperty("control_count");
        cursor.nextNumberProperty("control_positions");
        }
      }
      if (cursor.isProperty("psp_source")) cursor.nextStringProperty("psp_source");
      cursor.nextListProperty("freqs", Number.class);
      if (cursor.isProperty("sequences")) { // start sequences
        cursor.next("sequences"); // start sequences
	cursor.next(JSON_START_LIST);
	while (cursor.getEventType() != JSON_END_LIST) {
	  cursor.next(JSON_START_OBJECT);
	  cursor.nextStringProperty("name");
	  cursor.nextIntegerProperty("length");
	  cursor.nextNumberProperty("weight");
	  cursor.next(JSON_END_OBJECT);
	}
	cursor.next(JSON_END_LIST);
	cursor.next(JSON_END_PROPERTY); // end sequences
      }
      cursor.next(JSON_END_OBJECT);
      cursor.next(JSON_END_PROPERTY); // end sequence_db
    } else if (isDreme) {
      cursor.skipProperty("sequence_db");
      cursor.skipProperty("control_db");
    } else if (isStreme) {
      cursor.skipProperty("train_positives");
      cursor.skipProperty("train_negatives");
      cursor.skipProperty("test_positives");
      cursor.skipProperty("test_negatives");
      cursor.skipProperty("sequence_db");
    }

    cursor.next("motifs"); // start motifs
    cursor.next(JSON_START_LIST);
    while (cursor.getEventType() != JSON_END_LIST) {
      cursor.next(JSON_START_OBJECT);
      cursor.nextIntegerProperty("db");
      cursor.nextStringProperty("id");
      cursor.nextStringProperty("alt");
      int motifLen;
      if (isStreme) {
        motifLen = (int)cursor.nextIntegerProperty("width");
      } else {
        motifLen = (int)cursor.nextIntegerProperty("len");
	cursor.nextIntegerProperty("nsites");
	cursor.nextStringProperty("evalue");
      }
      if (isMeme) {
        cursor.nextNumberProperty("ic");
        cursor.nextNumberProperty("re");
        cursor.nextNumberProperty("llr");
        cursor.nextNumberProperty("bt");
        cursor.nextNumberProperty("time");
        cursor.nextMatrixProperty("psm", Number.class);
      } else if (isDreme) {
        cursor.nextNumberProperty("p");
        cursor.nextNumberProperty("n");
        cursor.nextStringProperty("pvalue");
        cursor.nextStringProperty("unerased_evalue");
      } else if (isStreme) {
        cursor.nextNumberProperty("initial_width");
        cursor.nextStringProperty("seed");
        cursor.nextNumberProperty("score_threshold");
        cursor.nextNumberProperty("train_pos_count");
        cursor.nextNumberProperty("train_neg_count");
        cursor.nextNumberProperty("train_log_pvalue");
        cursor.nextStringProperty("train_pvalue");
        cursor.nextNumberProperty("train_dtc");
        cursor.nextNumberProperty("train_bernoulli");
        cursor.nextNumberProperty("test_pos_count");
        cursor.nextNumberProperty("test_neg_count");
        cursor.nextNumberProperty("test_log_pvalue");
        cursor.nextStringProperty("test_pvalue");
        if (cursor.isProperty("test_log_evalue")) {	// post 5.4
          cursor.skipProperty("test_log_evalue");
        }
        if (cursor.isProperty("test_evalue")) {		// post 5.4
          cursor.skipProperty("test_evalue");
        }
        cursor.nextNumberProperty("test_dtc");
        cursor.nextNumberProperty("test_bernoulli");
        cursor.nextNumberProperty("elapsed_time");
	if (cursor.isProperty("total_sites")) {		// post 5.4
	  cursor.skipProperty("total_sites");
	}
	if (cursor.isProperty("site_distr")) {		// post 5.4
	  cursor.skipProperty("site_distr");
	}
	if (cursor.isProperty("max_sites")) {		// post 5.4
	  cursor.skipProperty("max_sites");
	}
	if (cursor.isProperty("site_hist")) {		// post 5.4
	  cursor.skipProperty("site_hist");
	}
        cursor.nextNumberProperty("len");
        cursor.nextNumberProperty("nsites");
        cursor.nextStringProperty("evalue");
      }
      cursor.nextMatrixProperty("pwm", Number.class);
      if (isMeme) {
        cursor.next("sites"); // start sites
        cursor.next(JSON_START_LIST);
        while (cursor.getEventType() != JSON_END_LIST) {
          cursor.next(JSON_START_OBJECT);
          cursor.nextIntegerProperty("seq");
          cursor.nextIntegerProperty("pos");
          cursor.nextBooleanProperty("rc");
          cursor.nextNumberProperty("pvalue");
          cursor.nextStringProperty("lflank");
          cursor.nextStringProperty("match");
          cursor.nextStringProperty("rflank");
          cursor.next(JSON_END_OBJECT);
        }
        cursor.next(JSON_END_LIST);
        cursor.next(JSON_END_PROPERTY); // end sites
      } else if (isDreme) {
        cursor.skipProperty("matches");
      }
      cursor.next(JSON_END_OBJECT);
      stats.addMotif(motifLen);
    }
    cursor.next(JSON_END_LIST);
    cursor.next(JSON_END_PROPERTY); // end motifs

    if (isMeme) {
      if (cursor.isProperty("scan")) { // start scan
        cursor.next("scan"); // start scan
	cursor.next(JSON_START_LIST);
	while (cursor.getEventType() != JSON_END_LIST) {
	  cursor.next(JSON_START_OBJECT);
	  cursor.nextNumberProperty("pvalue");
	  cursor.next("sites"); // start scan_sites
	  cursor.next(JSON_START_LIST);
	  while (cursor.getEventType() != JSON_END_LIST) {
	    cursor.next(JSON_START_OBJECT);
	    cursor.nextIntegerProperty("motif");
	    cursor.nextIntegerProperty("pos");
	    cursor.nextBooleanProperty("rc");
	    cursor.nextNumberProperty("pvalue");
	    cursor.next(JSON_END_OBJECT);
	  }
	  cursor.next(JSON_END_LIST);
	  cursor.next(JSON_END_PROPERTY); // end scan_sites
          cursor.next(JSON_END_OBJECT);
	}
	cursor.next(JSON_END_LIST);
	cursor.next(JSON_END_PROPERTY); // end scan
      }
    } else if (isDreme) {
      cursor.next("runtime"); // start runtime
      cursor.next(JSON_START_OBJECT);
      cursor.nextStringProperty("host");
      cursor.nextStringProperty("when");
      cursor.nextNumberProperty("cpu");
      cursor.nextNumberProperty("real");
      cursor.nextStringProperty("stop");
      cursor.next(JSON_END_OBJECT);
      cursor.next(JSON_END_PROPERTY); // end runtime
    } else if (isStreme) {
      cursor.nextStringProperty("stop_reason");
      cursor.next("runtime"); // start runtime
      cursor.next(JSON_START_OBJECT);
      cursor.nextStringProperty("host");
      cursor.nextNumberProperty("cpu");
      cursor.next(JSON_END_OBJECT);
      cursor.next(JSON_END_PROPERTY); // end runtime
    }
    cursor.next(JSON_END_OBJECT);
    cursor.expect(JSON_END_DATA);
    return stats;
  }

  public static MotifStats validateAsHTML(InputStream in, FeedbackHandler feedback) throws IOException, MotifParseException {
    MotifStats stats = null;
    Pattern versionPattern = Pattern.compile("^\\s*MEME\\s+version\\s+(\\d+)(?:\\.(\\d+)(?:\\.(\\d+))?)?\\s*$");
    Pattern pspmNamePattern = Pattern.compile("^pspm(\\d+)$", Pattern.CASE_INSENSITIVE);
    Pattern motifWidthPattern = Pattern.compile("\\sw\\s*=\\s*(\\d+(?:\\.\\d+)?)\\s");
    HDataPullParser cursor = null;
    try {
      cursor = new HDataPullParser(new InputStreamReader(in, "UTF-8"));
      while (cursor.getEventType() != HDataEventType.FILE_END) {
        switch (cursor.getEventType()) {
          case HIDDEN_FIELD: {
            if (stats == null) stats = new MotifStats();
            String name = cursor.getName();
            String value = (String)cursor.getValue();
            // other fields that are not checked for are: strands, bgfreq, name, combinedblock, nmotifs, motifname#, motifblock#, pssm#, BLOCKS#
            if ("version".equalsIgnoreCase(name)) {
              Matcher lineMatch = versionPattern.matcher(value);
              if (lineMatch.matches()) {
                int major = Integer.parseInt(lineMatch.group(1));
                int minor = (lineMatch.group(2) != null ? Integer.parseInt(lineMatch.group(2)) : 0);
                int patch = (lineMatch.group(3) != null ? Integer.parseInt(lineMatch.group(3)) : 0);
                if (!(major > 4 || (major == 4 && minor > 3) || (major == 4 && minor == 3 && patch >= 2))) return null;
              }
            } else if ("alphabet".equalsIgnoreCase(name)) {
              stats.setAlphabet(AlphStd.fromCore(value.trim()));
            } else if (pspmNamePattern.matcher(name).matches()) {
              Matcher m = motifWidthPattern.matcher(value);
              if (m.find()) stats.addMotif(Integer.parseInt(m.group(1)));
            }
          }
          break;
          case JSON_START_DATA: {
            try {
              return parseJson(cursor, feedback);
            } catch (AlphParser.AlphParseException e) {
	      feedback.whine("AlphParseException"); //debug
              throw new MotifParseException(e);
            } catch (HDataParseException e) {
	      feedback.whine("HDParseException"); //debug
              throw new MotifParseException(e);
            }
          }
        }
        cursor.next();
      }
      cursor.close();
    } finally {
      if (cursor != null) {
        try {
          cursor.close();
        } catch (IOException e) { /* ignore */ }
      }
    }
    return stats;
  }

  public static MotifStats validateAsText(InputStream inStream) throws IOException, MotifParseException {
    TextMotifParser parser = new TextMotifParser(); 
    BufferedReader in = null;
    try {
      in = new BufferedReader(new InputStreamReader(inStream, "UTF-8"));
      String line;
      while ((line = in.readLine()) != null) {
        parser.update(line);
      }
      in.close();
      in = null;
    } catch (AlphParser.AlphParseException e) {
      throw new MotifParseException(e);
    }finally {
      if (in != null) {
        try {
          in.close();
        } catch (IOException e) { /* ignore */ }
      }
    }
    return parser.getStats();
  }

  private static class TextMotifParser {
    private MotifStats stats;
    private int[] version;
    private int stage;
    private int alphLineCount;
    private AlphParser parser;

    private static final Pattern verLineRe = Pattern.compile(
        "^\\s*MEME\\s+version\\s+(\\d+)(?:\\.(\\d+)(?:\\.(\\d+))?)?");
    private static final Pattern oldAlphabetRe = Pattern.compile(
        "^\\s*ALPHABET\\s*=\\s*(\\S+)");
    private static final Pattern alphabetStartRe = Pattern.compile(
        "^\\s*ALPHABET\\s*(?:\"|$)", Pattern.CASE_INSENSITIVE);
    private static final Pattern alphabetEndRe = Pattern.compile(
        "^\\s*END\\s+ALPHABET\\s*$", Pattern.CASE_INSENSITIVE);
    private static final Pattern dividerRe = Pattern.compile("^\\*{80}$");

    private static final Pattern motifWidthRe = Pattern.compile(
        "^\\s*letter-probability\\s+matrix:.*\\s+w\\s*=\\s+(\\d+)");

    public TextMotifParser() {
      this.stats = new MotifStats();
      this.version = new int[]{0, 0, 0};
      this.stage = 0;
      this.parser = null;
      this.alphLineCount = 0;
    }

    public void update(String line) throws AlphParser.AlphParseException {
      Matcher m;
      switch (this.stage) {
        case 0:
          m = verLineRe.matcher(line);
          if (m.find()) {
            this.version[0] = Integer.parseInt(m.group(1));
            if (m.group(2) != null && m.group(2).length() > 0) {
              this.version[1] = Integer.parseInt(m.group(2));
            } else {
              this.version[1] = 0;
            }
            if (m.group(3) != null && m.group(3).length() > 0) {
              this.version[2] = Integer.parseInt(m.group(3));
            } else {
              this.version[2] = 0;
            }
            this.stage = 1;
          }
          break;
        case 1:
          if ((m = oldAlphabetRe.matcher(line)).find()) {
            String alphStr = m.group(1);
            if (alphStr.equals("ACGU")) {
              this.stats.setAlphabet(AlphStd.RNA.getAlph());
            } else if (alphStr.equals("ACGT")) {
              this.stats.setAlphabet(AlphStd.DNA.getAlph());
            } else if (alphStr.equals("ACDEFGHIKLMNPQRSTVWY")) {
              this.stats.setAlphabet(AlphStd.PROTEIN.getAlph());
            }
            this.stage = 3;
          } else if (alphabetStartRe.matcher(line).find()) {
            stage = 2;
            parser = new AlphParser();
            parser.parseLine(line);
            alphLineCount = 0;
          }
          break;
        case 2:
          if (alphabetEndRe.matcher(line).find() || (alphLineCount > 0 && dividerRe.matcher(line).find())) {
            this.stats.setAlphabet(parser.parseDone());
            stage = 3;
          } else {
            parser.parseLine(line);
            alphLineCount++;
          }
          break;
        case 3:
          m = motifWidthRe.matcher(line);
          if (m.find()) {
            int motifCols = Integer.parseInt(m.group(1));
            this.stats.addMotif(motifCols);
          }
          break;
      }
    }
    
    public MotifStats getStats() {
      if (this.stage != 3) return null;
      this.stats.freeze();
      return this.stats;
    }
  }

  public static MotifStats validate(File file, FeedbackHandler feedback) throws IOException {
  //public static MotifStats validate(File file) throws IOException {
    MotifStats stat = null;
    try {
      stat = validateAsHTML(new FileInputStream(file), feedback);
      if (stat != null) return stat;
    } catch (MotifParseException e) { /* ignore */ }
    try {
      stat = validateAsXML(new FileInputStream(file), feedback);
      if (stat != null) return stat;
    } catch (MotifParseException e) { /* ignore */ }
    try {
      stat = validateAsText(new FileInputStream(file));
    } catch (MotifParseException e) { /* ignore */ }
    return stat;
  }

  public static void main(String[] args) throws IOException, MotifParseException {
    if (args.length == 1) {
      MotifStats stats = validate(new File(args[0]), null);
      //MotifStats stats = validate(new File(args[0]));
      System.out.print(stats);
    }
  }
}
