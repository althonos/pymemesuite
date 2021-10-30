/**************************************************************************
 * FILE: xml-io.c
 * CREATE DATE: 8/3/2007
 * AUTHOR: William Stafford Noble, Tim Bailey, and Charles Grant
 * PROJECT: MEME
 * COPYRIGHT: 2007, University of Washington
 * DESCRIPTION: Utility functions for reading XML files.
 **************************************************************************/

#include <assert.h>
#include "alphabet.h"
#include "utils.h"
#include "xml-util.h"
#include <libxslt/xslt.h>
#include <libxslt/xsltInternals.h>
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

/***********************************************************************
 * Look up a set of elments using an XPath string.
 * Caller is responsible for freeing the xmlXPathObj.
 ***********************************************************************/
xmlXPathObjectPtr xpath_query(
  xmlXPathContextPtr xpath_ctxt,
  char* path
) {

  xmlXPathObjectPtr xpathObj = NULL;

  xpathObj = xmlXPathEvalExpression(BAD_CAST path, xpath_ctxt);
  if(xpathObj == NULL) {
      die("Error: unable to evaluate xpath expression %s.", path);
  }

  return xpathObj;

}

/***********************************************************************
 * Check if property from an XML node exists.
 ***********************************************************************/
bool check_xml_node_property(
  xmlNodePtr node, 
  char* property_name
) {
  bool exists;
  xmlChar* property = xmlGetProp(node, BAD_CAST property_name);
  exists = (property != NULL);
  free(property);
  return exists;
}

/***********************************************************************
 * Read a property from an XML node.
 * Caller is responsible for freeing the xmlChar string.
 ***********************************************************************/
xmlChar* read_xml_node_property(
  xmlNodePtr node, 
  char* property_name
) {

  xmlChar* property = xmlGetProp(node, BAD_CAST property_name);
  if (property == NULL) {
    die(
      "Error: unable to retreive property %s in tag %s.\n",
      property_name,
      node->name
    );
  }

  return property;

}

/***********************************************************************
 * Read and set the alphabet from a parsed XML document.
 ***********************************************************************/
ALPH_T* read_alphabet_from_xml(xmlXPathContextPtr xpath_ctxt) {
  ALPH_T* alph;

  xmlXPathObjectPtr xpathObj = NULL;
  xmlChar* property = NULL;

  xpathObj = xpath_query(xpath_ctxt, "//*/alphabet");
  property = read_xml_node_property(
    xpathObj->nodesetval->nodeTab[0], 
    "length"
  );
  int alph_size = atoi((char *) property);
  xmlFree(property);
  xmlXPathFreeObject(xpathObj);

  xpathObj = xpath_query(xpath_ctxt, "//*/alphabet/letter");
  // The stated size of the alphabet had better match the
  // number of letter elements in the alphabet.
  assert(alph_size == xpathObj->nodesetval->nodeNr);
  char* buffer = mm_calloc(sizeof(char), alph_size + 1);
  int i = 0;
  xmlNodePtr currLetterNode = NULL;
  for (i = 0; i < alph_size; i++) {
    currLetterNode = xpathObj->nodesetval->nodeTab[i];
    if (currLetterNode == NULL) {
      die("Error: missing letter %d in alphabet.\n", i);
    }
    // Get the letter symbol attribute
    property = read_xml_node_property(currLetterNode, "symbol");
    buffer[i] = *property;
    xmlFree(property);
  }
  buffer[i] = 0;
  alph = alph_type(buffer, alph_size+1);

  myfree(buffer);
  xmlXPathFreeObject(xpathObj);

  return alph;
}

/***********************************************************************
 * Read the background letter frequencies from XML.
 * Caller is responsible for freeing the returned array.
 ***********************************************************************/
ARRAY_T* read_bg_freqs_from_xml(xmlXPathContextPtr xpath_ctxt, ALPH_T* alph) {

  xmlXPathObjectPtr xpathObj = NULL;
  ATYPE    value;
  ARRAY_T* bg_freqs;

  int a_size = alph_size_core(alph);

  // Use XPATH to get the background frequencies from XML
  xpathObj = xpath_query(
    xpath_ctxt, 
    "//*/background_frequencies/alphabet_array/value"
  );
  int num_values = (xpathObj->nodesetval ? xpathObj->nodesetval->nodeNr : 0);
  xmlXPathFreeObject(xpathObj);

  // The number of background frequences should match the alphabet size.
  assert(num_values == a_size);

  // Allocate the array.
  bg_freqs= allocate_array(alph_size_full(alph));

  // XML doesn't enforce any order on the emission probability values,
  // so force reading bg frequency values in alphabet order.
  const int MAX_XPATH_EXPRESSION = 200;
  char xpath_expression[MAX_XPATH_EXPRESSION];
  xmlNodePtr currValueNode = NULL;
  int i_node = 0;
  for (i_node = 0; i_node < a_size; i_node++) {
    // Build the XPATH expression to get bg freq for a character.
    snprintf(
      xpath_expression,
      MAX_XPATH_EXPRESSION,
      "//*/background_frequencies/"
      "alphabet_array/value[@letter_id='letter_%c']",
      alph_char(alph, i_node)
    );
    // Read the selected bg frequency.
    xpathObj = xpath_query(xpath_ctxt, xpath_expression);
    // Should only find one node
    assert(xpathObj->nodesetval->nodeNr == 1);
    // Decode from node set to numeric value for bg freq.
    currValueNode = xpathObj->nodesetval->nodeTab[0];
    xmlXPathFreeObject(xpathObj);
    value = xmlXPathCastNodeToNumber(currValueNode);
    set_array_item(i_node, value, bg_freqs);
  }

  // Make sure the frequencies add up to 1.0. 
  normalize_subarray(0, a_size, 0.0, bg_freqs);

  // Fill in ambiguous characters. 
  calc_ambigs(alph, false, bg_freqs);

  return bg_freqs;

}

/**********************************************************************
print_xml_using_stylesheet

Print the contents of an XML Doc to a file, applying an 
XSLT stylesheet.
**********************************************************************/
bool print_xml_filename_to_file_using_stylesheet(
    char* input_file_path,      /* path to XML input file IN */
    char* stylesheet_file_path, /* path to MEME XSL stylesheet IN */
    FILE* output_file           /* path to HTML output file IN */
) {

  xsltStylesheetPtr stylesheet = NULL;
  xmlDocPtr input_doc = NULL;
  xmlDocPtr output_doc = NULL;
  const int PERFORM_ENTITY_SUBST = 1;

  xmlSubstituteEntitiesDefault(PERFORM_ENTITY_SUBST);
  input_doc = xmlParseFile(input_file_path);
  if (!input_doc) {
    fprintf(stderr, "Unable to parse input file %s.\n", input_file_path);
    return false;
  }
  stylesheet = xsltParseStylesheetFile((const xmlChar *) stylesheet_file_path);
  if (!stylesheet) {
    fprintf(stderr, "Unable to parse stylesheet %s.\n", stylesheet_file_path);
    return false;
  }
  output_doc = xsltApplyStylesheet(stylesheet, input_doc, NULL);
  if (!output_doc) {
    fprintf(
      stderr, 
      "Unable to apply stylsheet %s to to input.\n", 
      stylesheet_file_path
    );
    return false;
  }
  int result = xsltSaveResultToFile(output_file, output_doc, stylesheet);
  if (result == -1) {
    fprintf(stderr, "Unable to save result of applying stylesheet %st.\n", stylesheet_file_path);
    return false;
  }

  xmlFreeDoc(output_doc);
  xsltFreeStylesheet(stylesheet);

  return true;

} /* print_xml_to_file_using_stylesheet_html */

/**********************************************************************
print_xml_filename_to_filename_using_stylesheet

Print the contents of an XML file to another file applying an 
XSLT stylesheet.

Returns true if successful, false otherwise.
**********************************************************************/
bool print_xml_filename_to_filename_using_stylesheet(
    char* input_file_path,        /* path to XML input file IN */
    char* stylesheet_file_path,   /* path to MEME XSL stylesheet IN */
    char* output_file_path        /* path to HTML output file IN */
) {

  xsltStylesheetPtr stylesheet = NULL;
  xmlDocPtr input_doc = NULL;
  xmlDocPtr output_doc = NULL;
  const int PERFORM_ENTITY_SUBST = 1;

  xmlSubstituteEntitiesDefault(PERFORM_ENTITY_SUBST);
  xmlLoadExtDtdDefaultValue = 0;
  stylesheet = xsltParseStylesheetFile((const xmlChar *) stylesheet_file_path);
  if (!stylesheet) {
    fprintf(stderr, "Unable to parse stylesheet %s.\n", stylesheet_file_path);
    return false;
  }
  input_doc = xmlParseFile(input_file_path);
  if (!input_doc) {
    fprintf(stderr, "Unable to parse input file %s.\n", input_file_path);
    return false;
  }
  output_doc = xsltApplyStylesheet(stylesheet, input_doc, NULL);
  if (!output_doc) {
    fprintf(
      stderr, 
      "Unable to apply stylsheet %s to input from file %s.\n", 
      stylesheet_file_path,
      input_file_path
    );
    return false;
  }
  int result = xsltSaveResultToFilename(output_file_path, output_doc, stylesheet, 0);
  if (result == -1) {
    fprintf(
      stderr, 
      "Unable to save result of applying stylesheet %s to %s.\n", 
      stylesheet_file_path, 
      output_file_path
    );
  }

  xsltFreeStylesheet(stylesheet);
  xmlFreeDoc(output_doc);
  xmlFreeDoc(input_doc);
  xsltCleanupGlobals();
  xmlCleanupParser();

  return true;

} /* print_xml_file_html */
