/*
 * This file defines documentation that is repeated in many files and so needs
 * to be centralized. This is the equivalent of the old options.xsl file except
 * it actually works with HTML5 and allows variable substitution, conditional
 * sections and lists.
 */

var shared_doc = {

  /* Fairly common options  */

  "all-desc-opt": ["desc"],
  "all-desc-itm": ["<span class='pdat'>description</span>"],
  "all-desc-desc": [
    "Include the text <span class='pdat'>description</span> in the \"Description\" section of the " +
    "HTML output, taking into account newlines (line breaks) and multiple newlines (paragraph breaks)."
  ],
  "all-desc-def": [
    "No description in the HTML output."
  ],

  "all-dfile-opt": ["dfile"],
  "all-dfile-itm": ["<span class='pdat'>dfile</span>"],
  "all-dfile-desc": [
    "Include the first 500 characters of text from the file named " + 
    "<span class='pdat'>dfile</span> in the \"Description\" secton of the HTML output, " +
    "taking into account newlines (line breaks) and multiple newlines (paragraph breaks). " +
    "This option overrides option <span class='popt'>",
    {"if": "ddash", "yes": "-"},
    "-desc</span>."
  ],
  "all-dfile-def": [
    "No description in the HTML output."
  ],

  "all-help-opt": [{"if": "short", "yes": "h", "no": "help"}],
  "all-help-desc": [
    "Display the usage message and exit."
  ],
  "all-help-def": [
    "Run as normal."
  ],

  "all-version-opt": ["version"],
  "all-version-desc": [
    "Display the program version and exit."
  ],
  "all-version-def": [
    "Run as normal."
  ],

  "all-text-opt": ["text"],
  "all-text-desc": [
    "Limits output to ",
    {"if": "tsv", "yes": "TSV (tab-separated values) formatted results ", "no": "plain text results "},
    "sent to standard output. ",
    {"if": "unsorted", "yes": "The results are unsorted and no q-values are output, allowing very large files to be searched."},
  ],
  "all-text-def": [
    "All formats are output to files in the selected output directory."
  ],

  /* psp */
  "all-psp-opt": ["psp"],
  "all-psp-val": ["file"],
  "all-psp-desc": [
    "File containing position specific priors (PSP) in ",
    "<a href='psp-format.html'>MEME PSP format</a> ",
    "or ",
    "<a href='http://genome.ucsc.edu/goldenPath/help/wiggle.html'>wiggle format</a>. ",
    "This will cause positions with higher prior values to receive larger motif scores. ",
    "This file can be generated ",
    "using the <a href='create-priors.html'>create-priors</a> utility. ",
    "Requires option <span class='popt'>--prior-dist</span>."
  ], 
  "all-psp-def": [
    "A uniform position-specific prior is used."
  ],

  /* prior-dist */
  "all-pdist-opt": ["prior-dist"],
  "all-pdist-val": ["file"],
  "all-pdist-desc": [
     "File containing binned distribution of priors. ",
     "This file can be generated using the ",
     "<a href='create-priors.html'>create-priors</a> utility."
  ],
  "all-pdist-def": [
    "A uniform position-specific prior is used."
  ],


  /* parse-genomic-coord */
  "all-pgc-opt": ["parse-genomic-coord"],
  "all-pgc-desc": [
    "When this option is specified, each FASTA sequence header is ",
    "checked for UCSC style genomic coordinates (e.g., <code>chr1:156887119-156887619</code>). ",
    "The sequence ID in the FASTA header should have the form:",
    "<div style='margin: 5px 0'>",
    "  &gt;<span class='pdat'>sequence name</span>:<span class='pdat'",
    "  >starting position</span>-<span class='pdat'>ending position</span>",
    "</div>",
    "Where",
    "<ul style='margin-top: 0; margin-bottom: 5px'>",
    "  <li><span class='pdat'>sequence name</span> is the name of the sequence,</li>",
    "  <li><span class='pdat'>starting position</span> is the genomic position of the first base and</li>",
    "  <li><span class='pdat'>ending position</span> is the genomic position of the final base.</li>",
    "</ul>",
    "The <span class='pdat'>sequence name</span> may not contain any white space. ",
    "If genomic coordinates are found they is used as the ",
    "coordinates in the output. When no coordinates are found, the default ",
    "behavior is used."
  ],
  "all-pgc-def": [
    "The first position in the sequence is assumed to be 1."
  ],

  /* inc */
  "all-inc-opt": [{"if": "inc", "yes": "inc", "no": "exc"}],
  "all-inc-val": ["pattern"],
  "all-inc-desc": [
    {"if": "inc", "yes": "Use <b>only</b> the ", "no": "Exclude any "},
    {"if": "spamo", "yes": "<b>secondary</b> "},
    {"if": "centrimo", "yes": ""},
    "motifs with names matching the pattern. The pattern can ",
    "contain shell-like wildcards (e.g., \"<code>*</code>\") though they must be escaped ",
    "or quoted to prevent the shell from auto-expanding them. ",
    "This option may be may be repeated, and ",
    {"if": "inc", "yes": "<b>only</b> ", "no": "any "},
    "motifs with names matching at least one pattern is ",
    {"if": "inc", "yes": "used. ", "no": "excluded. "},
    {"if": "centrimo", "yes": [
        {"if": "inc", "yes": "Note: <code>--motif</code> is a synonym for <code>--inc</code>."}
      ]
    },
  ],
  "all-inc-def": [
    {"if": "inc", 
      "yes": ["Unless the <span class='popt'>",
          {"if": "ddash", "yes": "-"},
          "-exc</span> option has been specified, all the motifs are used."
        ],
      "no": ["Unless the <span class='popt'>",
          {"if": "ddash", "yes": "-"},
          "-inc</span> option has been specified, all the motifs is used."
        ]
    }
  ],

  /* bfile/bgfile */
  "all-bg-opt": [{"if": "bfile", "yes": "bfile", "no": "bgfile"}],
  "all-bg-val": ["file"],
  "all-bg-desc": ["Specify the source of a 0-order background model for converting a frequency matrix to a log-odds score matrix ",
    " and for use in estimating the <i>p</i>-values of match scores. ",
    "The background model normalizes for biased distribution of individual letters in the sequences. ",
    "The value of <span class=\"pdat\">file</span> is either the path to a file in ",
    "<a href=\"bfile-format.html\">Markov Background Model Format</a>, or one of the keywords ",
    "<code>--motif--</code>, <code>motif-file</code> or <code>--uniform--</code>. ",
    "The first two keywords cause the 0-order letter frequencies contained in the ",
    {"if": "spamo", "yes": "<b>primary</b> "},
    {"if": "centrimo", "yes": "<b>first</b> "},
    "motif file to be used, ",
    "and <code>--uniform--</code> causes uniform letter frequencies to be used. ",
    {"if": "mast", "yes": "Note: If the motif file contains score matrices, then the background model is used only in estimating <i>p</i>-values of match scores."},
    "If the background model in <span class=\"pdat\">file</span> is higher than 0-order, only the 0-order portion is used. ",
    {"if": "sym", 
      "yes": "The background model is modified by averaging the frequencies of letters and their reverse complements. ",
      "no": "If both strands are being scored, the background model is modified by averaging the frequencies of letters and their reverse complements."
    }
  ],
  "all-bg-def": [
    {"if": "sequences",
      "yes": "The frequencies of the letters in the (primary) <b>sequences</b> are used as the 0-order background model.",
      "no": [ 
        {"if": "dna_only",
          "yes": "Frequencies embedded in the application (from an early version of the non-redundant DNA) are used.", 
          "no": "If the sequence alphabet is DNA or protein, frequencies embedded in the application (from an early version of the non-redundant DNA or protein database, respectively) are used. Otherwise, the letter frequencies contained in the motif file are used."
        }
      ]
    }
  ],

  /* pseudocount */
  "all-pseudo-opt": ["motif-pseudo"],
  "all-pseudo-val": ["pseudocount"],
  "all-pseudo-desc": [
    "Add a this total pseudocount to the counts in each motif column when converting a frequency matrix to a log-odds score matrix. ",
    "The pseudocount added to each count is <span class=\"pdat\">pseudocount</span> times the background frequency of the letter (see option ",
    "<span class=\"popt\">", 
      {"if": "ddash", "yes": "--", "no": "-"},
      {"if": "bfile", "yes": "bfile", "no": "bgfile"}, 
    "</span>, above). ",
    "<b>Note:</b> Counts are computed from <a href=meme-format.html>MEME formatted motifs</a> by multiplying ",
    "the frequency of the letter times the value of <code>nsites</code> given in the motif <code>letter-probability matrix</code> header line. ",
      {"if": "legacy", "yes": "<b>Note:</b> The synonym <span class=\"popt\">--pseudocount</span> is also allowed. "},
  ],
  "all-pseudo-def": [
    "The program applies a pseudocount of 0.1."
  ],

  /* alph */
  "all-alph-opt": ["alph"],
  "all-alph-val": ["alphabet file"],
  "all-alph-desc": ["The <span class=\"pdat\">alphabet file</span> contains the ",
    "<a href=\"alphabet-format.html\">alphabet definition</a> of a (possibly non-standard) alphabet",
    {"if": "meme",
      "yes": " and may contain the additional wildcard symbols <code>*</code> or <code>-</code>. ",
      "no": ". "
    },
    {"if": "dreme", 
      "yes": "Note that DREME works best when there are ambiguous symbols for all likely combinations of core symbols. As DREME is currently implemented it can only access ambiguous symbols that are the union of two other (possibly ambiguous) symbols. ",
      "no": [
        {"if": "fasta-center",
          "yes": "The sequences are verified using the alphabet definition, and sequences consisting of nothing but ambiguous symbols are rejected and optionally written to the reject file. ",
          "no": [
            {"if": "fasta-shuffle",
	      "yes": "Alias symbols are converted to their core representation in the given <a href='alphabet-format.html'>alphabet definition</a> prior to shuffling. Unrecognized symbols are converted to the alphabet's wildcard. ",
              "no": [
		{"if": "getsize",
 	 	  "yes": "Letter frequencies are printed in <a href='bfile-format.html'>background file format</a>. Alias symbols are converted to their core representation in the given <a href='alphabet-format.html'>alphabet definition</a> and ambiguous characters are ignored. ",
 	 	  "no": ""
                }
              ]
            }
          ]
        }
      ]
    },
    "Incompatible with options <span class=\"popt\">", {"if": "ddash", "yes": ["-"]}, "-dna</span>, <span class=\"popt\">",
      {"if": "ddash", "yes": ["-"]}, "-rna</span> and <span class=\"popt\">",
      {"if": "ddash", "yes": ["-"]}, "-protein</span>. "
  ],
  "all-alph-def": [
    {"if": "meme",
      "yes": "The sequences are assumed to be in the protein alphabet and may contain the additional wildcard symbols <code>*</code> or <code>-</code>. ",
      "no": [
        {"if": "fasta-center", 
          "yes": "When no alphabet is specified, the sequence content is ignored. ",
          "no": [
            {"if": "get-markov",
	      "yes": "Autodetect the alphabet (DNA, RNA or protein) from the input sequences. ",
              "no": [
                {"if": "ambig-distinct",
	          "yes": "When no alphabet is specified, symbols (including aliases) are treated as distinct symbols. ",
		  "no": [
		    {"if": "getsize",
		      "yes": "",
	              "no": "The sequences are assumed to be in the DNA alphabet. "
		    }
                  ]
                }
              ]
            }
          ]
        }
      ]
    },
    {"if": "ambig-unknown",
      "yes" : "All ambiguous letters are treated as \"unknown\".",
      "no": "",
    }
  ],

  /* dna, rna, protein */
  "all-stdalph-opt": [
    {"if": "dna", 
      "yes" : "dna",
      "no" : [
        {"if": "rna",
	  "yes" : "rna",
	  "no" : "protein"
        }
      ]
    }
  ],
  "all-stdalph-desc": [
    {"if": "experimental",
      "yes": "<span class='experimental'>Not recommended</span> ",
      "no": ""
    },
    "Same as specifying the <span class='popt'>",
    {"if": "ddash", "yes": ["-"]}, "-alph</span> option with the standard <a href='alphabet-format.html#standard_",
    {"if": "dna",
      "yes" : "DNA'>DNA alphabet definition file</a>. ",
      "no" : [
        {"if": "rna",
          "yes" : "RNA'>RNA alphabet definition file</a>. ",
          "no" : "protein'>protein alphabet definition file</a>. ",
        }
      ]
    },
    {"if": "ambig-unknown",
      "yes": "All ambiguous characters are treated as \"unknown\". ",
      "no": ""
    }, 
    "Incompatible with options ",
    {"if": "dna", 
      "yes" : ["<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-alph</span>, ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-rna</span>, and ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-protein</span>. "],
      "no" : [
        {"if": "rna",
          "yes" : ["<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-alph</span>, ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-dna</span>, and ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-protein</span>. "],
          "no" : ["<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-alph</span>, ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-dna</span>, and ",
               "<span class='popt'>", {"if": "ddash", "yes": ["-"]}, "-rna</span>."]
        }
      ]
    },
    {"if": "experimental",
      "yes": " This option is <b>not</b> recommended for use with this program. ",
      "no": ""
    },
  ],
  "all-stdalph-def": [
    {"if": "streme",
      "yes": "The sequences are assumed to be in the DNA alphabet and all ambiguous letters are treated as \"unknown\". ",
      "no": [
	{"if": "meme",
	  "yes": "The sequences are assumed to be in the protein alphabet and all ambiguous letters are treated as \"unknown\". ",
	  "no": [
	    {"if": "ignore",
	      "yes": "When no alphabet is specified, the sequence content is ignored. ",
	      "no": [
		{"if": "auto",
		  "yes": "Autodetect the alphabet (DNA, RNA or protein) from the input sequences.",
		  "no": [
		    {"if": "ambig-distinct",
		      "yes": "When no alphabet is specified, symbols (including aliases) are treated as distinct symbols while shuffling. ",
		      "no": [
			{"if": "getsize",
			  "yes": "",
			  "no": "The sequences are assumed to be in the DNA alphabet. "
			}
		      ]
		    }
		  ]
		}
	      ]
	    }
          ]
	}
      ]
    }
  ],

  /* [x]alph */
  "all-x_alph-opt": ["[x]alph"],
  "all-x_alph-val": ["alphabet file"],
  "all-x_alph-desc": ["The <span class=\"pdat\">alphabet file</span> contains the ",
    "<a href=\"alphabet-format.html\">alphabet definition</a> of a (possibly non-standard) alphabet. ",
    "The input <b>sequences</b> must be in this alphabet. ",
    "If the input <b>motifs</b> are in a different alphabet than the input sequences, ",
    "specifying <span class='popt'>",
    {"if": "ddash", "yes": ["-"]},
    "-xalph</span> will cause the input motifs to be converted to this new alphabet, ",
    "with the probabilities for the new symbols set to zero prior to applying pseudocounts."
  ],
  "all-x_alph-def": [
     "The sequences are assumed to be in the DNA alphabet and motifs retain the alphabet defined in the motif file.",
  ],

  /* dna2rna */
  "all-dna2rna-opt": ["dna2rna"],
  "all-dna2rna-desc": [
    "The input DNA sequences will be treated as RNA (treating <code>T</code> as representing <code>U</code>), ",
    "and the output motifs will use the RNA alphabet. ",
    "If you provide any known motifs using the ",
    {"if": "memechip",
       "yes": [
         "<span class='popt'>-db</span> "
       ],
       "no": [
         "<span class='popt'>--m</span> "
       ],
    },
    "option, they must use the RNA alphabet or you must also specify ",
    {"if": "memechip",
       "yes": [
         "<span class='popt'>-xalph</span> "
       ],
       "no": [
         "<span class='popt'>--xalph</span> "
       ],
    },
    "with an RNA alphabet file."
  ],
  "all-dna2rna-def": [
    "Treat input DNA sequences as DNA."
  ],

  /* xalph */
  "all-xalph-opt": ["xalph"],
  "all-xalph-val": [
    {"if": "nofile",
      "yes": "",
      "no": ["alphabet file"]
    }
  ],
  "all-xalph-desc": [
    {"if": "spamo",
      "yes": [ 
        "Convert the alphabet of the secondary motif database(s) to the alphabet of the primary ",
        "motif assuming the core symbols of the secondary motif alphabet are a subset. ",
        "The input sequences must be in the alphabet of the primary motif."
      ],
      "no": [ 
        {"if": "tomtom",
          "yes": [
            "Convert the alphabet of the target motif database(s) to the alphabet of the ",
            "query motif database assuming the core symbols of the target motif alphabet are a subset ",
	    "of those used by the query motif(s)."
          ],
	  "no": [
            {"if": "2meme",
              "yes": [
                "Convert all motifs to use the same alphabet as specified in the ", 
                "first motif file. The alphabet specified in the first motif file ", 
                "must contain a superset of all the core symbols in the other ", 
                "motif file alphabets. ", 
              ],
	      "no": [
		"If the input motifs are in a different alphabet than the input sequences, ",
		"and the motif alphabet is a subset of the sequence alphabet, ",
		"you can specify an <span class=\"pdat\">alphabet file</span> containing the sequence ",
		"<a href=\"alphabet-format.html\">alphabet definition</a>. ",
		"The input motifs are converted to this new alphabet, with the probabilities ",
		" for the new symbols set to zero prior to applying pseudocounts."
	      ]
            }
          ]
        }
      ]
    }
  ],
  "all-xalph-def": [
    {"if": "spamo",
      "yes": ["The primary motif alphabet must be identical to the secondary motif alphabet and to that of the input sequences."],
      "no": [
        {"if": "tomtom",
          "yes": ["The query motif alphabet must be identical to the target motif alphabet."],
          "no": ["Motifs retain the alphabet defined in the motif file."]
        }
      ]
    }
  ],

  /* verbosity */
  "all-verbosity-opt": ["verbosity"],
  "all-verbosity-itm": [{"lst": ["1", "2", "3", "4", "5"]}],
  "all-verbosity-desc": ["A number that regulates the verbosity level of the ",
    "output information messages. If set to <span class=\"popt\">1</span> ",
    "(quiet) then it will only output error messages whereas the other ",
    "extreme <span class=\"popt\">5</span> (dump) outputs lots of mostly ",
    "useless information."],
  "all-verbosity-def": ["The verbosity level is set to ",
    "<span class=\"popt\">2</span> (normal verbosity)."],

  /* output without clobber */
  "all-o-opt": ["o"],
  "all-o-val": ["dir"],
  "all-o-desc": [
    "Create a folder called <span class=\"pdat\">dir</span> and write output ",
    "files in it. This option is not compatible with ",
    "<span class=\"popt\">", {"if": "ddash", "yes": ["-"]} ,"-", {"if": "glam2", "yes": "O", "no": "oc"}, "</span> as ",
    "only one output folder is allowed."],
  "all-o-def": [
    {"if": "dir",
      "no": ["The program writes to standard output."], 
      "yes": ["The program behaves as if <code>",
        {"if": "ddash", "yes": ["-"]} ,"-", {"if": "glam2", "yes": "O", "no": "oc"}, " ", {"key": "dir", "def": "dir_out"},
        "</code> had been specified."]
    }],

  /* output with clobber */
  "all-oc-opt": [{"if": "glam2", "yes": "O", "no": "oc"}],
  "all-oc-val": ["dir"],
  "all-oc-desc": [
    "Create a folder called <span class=\"pdat\">dir</span> but if it already ",
    "exists allow overwriting the contents. This option is not compatible ",
    "with <span class=\"popt\">", {"if": "ddash", "yes": ["-"]} ,"-o</span> ",
    "as only one output folder is allowed."],
  "all-oc-def": [
    {"if": "dir",
      "no": ["The program writes to standard output."],
      "yes": ["The program behaves as if <code>",
      {"if": "ddash", "yes": ["-"]}, "-", {"if": "glam2", "yes": "O", "no": "oc"}, " ", {"key": "dir", "def": "dir_out"},
      "</code> had been specified."]
    }],


  /* For format2meme conversion scripts, hence the 2meme prefix */
  "2meme-truncate-opt": ["truncate_names"],
  "2meme-truncate-desc": ["Truncate the motif names (IDs) at the first underscore. (E.g., if the ID is \"Arid3a_primary\", then the motif name will be \"Arid3a\". "],
  "2meme-truncate-def": ["Use the full motif names. "],

  "2meme-numseqs-opt": ["numseqs"],
  "2meme-numseqs-val": ["count"],
  "2meme-numseqs-desc": ["Assume frequencies based on <span class='pdat'>count</span> sequence sites. "],
  "2meme-numseqs-def": ["The motif is created as if it was made from 20 sequence sites."],

  /* output description */
  "2meme-output": [
    "<p>Writes <a href=\"meme-format.html\">MEME motif format</a> to standard ",
    "output.</p>",
    "<p>A probability matrix and optionally a log-odds matrix are output for ",
    "each ", {"if": "sequence", "no": "motif", "yes": "sequence"}, " in the file. ",
    {"if": "meme", "no": [ "The probability matrix is computed using ",
    "pseudo-counts consisting of the background frequency (see ",
    "<span class=\"popt\">-bg</span>, below) multiplied by the total ",
    "pseudocounts (see <span class=\"popt\">-pseudo</span>, below). "]},
    "The log-odds matrix uses the background frequencies in the denominator ",
    "and is log base 2.</p>"
    ],

  /* skip */
  "2meme-skip-opt": ["skip"],
  "2meme-skip-val": ["ID"],
  "2meme-skip-desc": [
    "Skip the motif identified by <span class=\"pdat\">ID</span>. ",
    "This option can be repeated to skip more than one motif."
    ],
  "2meme-skip-def": ["Motifs are not skipped."],
  /* numbers */
  "2meme-numbers-opt": ["numbers"],
  "2meme-numbers-desc": ["Use a number based on the position in the input ",
    "instead of the ", {"key": "id_source", "def": "motif name"}, " as the ",
    "motif identifier."],
  "2meme-numbers-def": ["The ", {"key": "id_source", "def": "motif name"},
    " is used as the motif name."],
  /* bg */
  "2meme-bg-opt": ["bg"],
  "2meme-bg-val": ["background file"],
  "2meme-bg-desc": [
    "The <span class=\"pdat\">background file</span> should be a ",
    "<a href=\"bfile-format.html\">Markov background model</a>. ",
    {"if" : "meme", 
      "no": ["It contains the background frequencies of letters use for assigning pseudocounts. "]
    }, "The background frequencies are included in the resulting MEME file."
    ],
  "2meme-bg-def": [
    {
      "if": "meme",
      "yes": ["Uses background in first specified motif file."],
      "no": ["Uses uniform background frequencies."]
    }],
  /* pseudo */
  "2meme-pseudo-opt": ["pseudo"],
  "2meme-pseudo-val": ["total pseudocounts"],
  "2meme-pseudo-desc": [
    "Add <span class=\"pdat\">total pseudocounts</span> times letter ",
    "background to each frequency."],
  "2meme-pseudo-def": ["No pseudocount is added."],
  /* log odds */
  "2meme-logodds-opt": ["logodds"],
  "2meme-logodds-desc": [
    "Include a log-odds matrix in the output. This is not required for ",
    "versions of the MEME Suite &ge; 4.7.0."],
  "2meme-logodds-def": ["The log-odds matrix is not included in the output."],
  /* url */
  "2meme-url-opt": ["url"],
  "2meme-url-val": ["website"],
  "2meme-url-desc": [
    "The provided website URL is stored with the motif ",
    {"if": "may_have_url", "yes": ["(unless the motif already has a URL) "]},
    "and this can be ",
    "used by MEME Suite programs to provide a direct link to that information ",
    "in their output. If <span class=\"pdat\">website</span> contains the ",
    "keyword MOTIF_NAME the ", {"key": "type", "def": "motif ID"}, " is ",
    "substituted in place of MOTIF_NAME in the output.<br> ",
    "For example if the url is <div class=\"indent cmd\">",
    {"key": "url_start", "def": "http://big-box-of-motifs.com/motifs/"},
    "MOTIF_NAME", {"key": "url_end", "def": ".html"}, "</div> and the ",
    {"key": "type", "def": "motif ID"}, " is <code>", 
    {"key": "example", "def": "motif_id"}, "</code>, the motif will ",
    "contain a link to <div class=\"indent cmd\">",
    {"key": "url_start", "def": "http://big-box-of-motifs.com/motifs/"},
    {"key": "example", "def": "motif_id"}, {"key": "url_end", "def": ".html"},
    "</div>", {"if": "is_transfac", "yes": [" Similarly, the TRANSFAC accession number is substituted for MOTIF_AC and the TRANSFAC ID is substituted for MOTIF_ID in the URL string."]}],
  "2meme-url-def": [
    "The output does not include a URL with the motifs", 
    {"if": "may_have_url", "yes": [" unless the original motif already had a URL"]},"."]

};

function wrdoc2(doc, subs) {
  var i, part;
  if (typeof doc === "string") {
    document.write(doc);
    return;
  }
  for (i = 0; i < doc.length; i++) {
    part = doc[i];
    if (typeof part === "string") {
      document.write(part);
    } else if (typeof part === "object") {
      if (part.hasOwnProperty("key")) {
        if (typeof subs === "object" && typeof subs[part["key"]] === "string") {
          document.write(subs[part["key"]]);
        } else if (part.hasOwnProperty("def")) {
          document.write(part["def"]);
        }
      } else if (part.hasOwnProperty("if")) {
        if (typeof subs === "object" && subs[part["if"]]) {
          if (part.hasOwnProperty("yes")) wrdoc2(part["yes"], subs);
        } else {
          if (part.hasOwnProperty("no")) wrdoc2(part["no"], subs);
        }
      } else if (part.hasOwnProperty("lst")) {
        var list, lside, rside, sep, j;
        // set defaults
        list = [];
        lside = "<span class=\"popt\">";
        rside = "</span>";
        sep = "|&#8203;";
        if (typeof part["lst"] === "string") {
          if (typeof subs === "object" && typeof subs[part["lst"]] === "object") {
            list = subs[part["lst"]];
          }
        } else if (typeof part["lst"] === "object") {
          list = part["lst"];
        }
        if (typeof part["lside"] === "string") {
          if (typeof subs === "object" && typeof subs[part["lside"]] === "object") {
            lside = subs[part["lside"]];
          }
        } else if (typeof part["lside"] === "object") {
          lside = part["lside"];
        }
        if (typeof part["rside"] === "string") {
          if (typeof subs === "object" && typeof subs[part["rside"]] === "object") {
            rside = subs[part["rside"]];
          }
        } else if (typeof part["rside"] === "object") {
          rside = part["rside"];
        }
        if (typeof part["sep"] === "string") {
          if (typeof subs === "object" && typeof subs[part["sep"]] === "object") {
            sep = subs[part["sep"]];
          }
        } else if (typeof part["sep"] === "object") {
          sep = part["sep"];
        }
        for (j = 0; j < list.length; j++) {
          if (j !== 0) wrdoc2(sep, subs);
          wrdoc2(lside, subs);
          wrdoc2(list[j], subs);
          wrdoc2(rside, subs);
        }

      }
    }
  }
}

function wrdoc(doc_key, subs) {
  if (shared_doc.hasOwnProperty(doc_key)) {
    wrdoc2(shared_doc[doc_key], subs);
  }
}

function wropt(opt_key, subs) {
  document.write("<tr><td class=\"popt\">");
  if (typeof subs === "object" && subs["ddash"]) {
    document.write("--");
  } else {
    document.write("-");
  }
  wrdoc(opt_key + "-opt", subs);
  document.write("</td><td>");
  if (shared_doc.hasOwnProperty(opt_key + "-val")) {
    if (typeof subs === "object" && subs["nospan"]) {
    } else {
      document.write("<span class=\"pdat\">");
    }
    wrdoc(opt_key + "-val", subs);
    if (typeof subs === "object" && subs["nospan"]) {
    } else {
      document.write("</span>");
    }
  } else if (shared_doc.hasOwnProperty(opt_key + "-itm")) {
    wrdoc(opt_key + "-itm", subs);
  }
  document.write("</td><td>");
  wrdoc(opt_key + "-desc", subs);
  document.write("</td><td>");
  wrdoc(opt_key + "-def", subs);
  document.write("</td></tr>");
}
