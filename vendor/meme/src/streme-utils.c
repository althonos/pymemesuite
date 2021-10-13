#include "streme-utils.h"

//
// Create a 0-order background model (frequencies) and
// higher-order background model (log conditional probabilities)
// from all the *control* sequences and store them in the multiseq object.
//
void set_multiseq_background(
  STREME_OPTIONS_T *options,
  Multiseq *multiseq,
  ALPH_T *alph,
  BOOL do_rc,
  ARRAY_T *back			// MEME-style background model (or NULL)
) {
  int i;
  Uint alen = alph_size_core(alph);

  int order = options->order;
  if (back != NULL) {
    // SEA: Use background from motifs or uniform.
    DEBUG_FMT(NORMAL_VERBOSE, "# Background model from %s.\n", options->bfile);
  } else if (options->bfile) {
    // Read in the MEME-style background model.
    DEBUG_FMT(NORMAL_VERBOSE, "# Background model by --bfile %s.\n", options->bfile);
    back = load_markov_model(options->alph, &order, options->bfile);
    if (order < options->order) {
      DEBUG_FMT(QUIET_VERBOSE,
        "# Warning: Background order %d requested but the given --bfile is only order %d.\n"
        "#          Setting --order to %d.\n", options->order, order, order);
      options->order = order;
    }
  } else {
    // Create a MEME-style higher-order background from the control sequences.
    // Rereads the negative sequence file (or the posfile if no negfile).
    DEBUG_MSG(NORMAL_VERBOSE, "# Estimating background model from control sequences.\n");
    back = get_markov_from_sequences(
      options->negfile ? options->negfile : options->posfile,	// all the negative sequences
      &order,
      MEME_BG_PRIOR,
      alph,
      options->alph_file,
      options->alphabet_type,
      do_rc
    );
  }

  //
  // Convert tuple probabilities to conditionals: if word ending at i is "wa"
  //   back[i] = Pr(a | w)
  int bsize = get_array_length(back);
  for (i = 0; i < bsize; i += alen) {
    normalize_subarray(i, alen, 1e-7, back);
  }

  //
  // Save 0-order background model.
  //
  multiseq->background = (double *) malloc(alen * sizeof(double));
  for (i = 0; i < alen ; i++) {
    multiseq->background[i] = get_array_item(i, back);
  }

  //
  // Convert from ARRAY to array of double log_2 values for speed.
  //
  multiseq->lcbp = (double *) malloc(bsize * sizeof(double));
  multiseq->bg_order = order;
  for (i = 0; i < bsize; i++) {
    multiseq->lcbp[i] = log(get_array_item(i, back)) / log(2);
  }

  // Print the 0-order (portion of the) background to standard out.
  DEBUG_MSG(NORMAL_VERBOSE, "# Background:");
  for (i=0; i<alen; i++) DEBUG_FMT(NORMAL_VERBOSE, " %c %.3g", I2A(i), multiseq->background[i]);
  DEBUG_MSG(NORMAL_VERBOSE, "\n");
  DEBUG_FMT(NORMAL_VERBOSE, "# Background order: %d Background size: %d\n", multiseq->bg_order, bsize);

  // Free the temporary array.
  free_array(back);
} // set_multiseq_background

//
// Store the log cumulative (higher-order) background probability of each sequence
// in the multisequence object.
//
void set_logcumback(
  STREME_OPTIONS_T *options,
  Multiseq *multiseq
) {
  int i, j, dir, ndir, seqno, pos;
  BOOL do_rc = multiseq->do_rc;
  Uint npos = multiseq->npos;
  Uint nneg = multiseq->nneg;
  Uint ntot = npos + nneg;
  ALPH_T *alph = multiseq->alph; 
  int bg_order = multiseq->bg_order;
  double *background = multiseq->background;
  double *lcbp = multiseq->lcbp;		// (kmer-1)-order log_2 background conditional probabilites Pr(a | w)
  Uint alen = alph_size_core(alph);		// length of alphabet
  double *logcumback = (double *) malloc(sizeof(double) * multiseq->totallength);
  multiseq->logcumback = logcumback;

  // Compute cumulative probability for forward (and possibly reverse) sequences.
  ndir = do_rc ? 2 : 1;
  for (dir=0; dir<ndir; dir++) {
    for (seqno=0; seqno<ntot; seqno++) {
      Uint seqstart = multiseq->seqstarts[seqno] + (dir==0 ? 0 : multiseq->totallength/2 + 1);
      Uint seqlen = multiseq->seqlengths[seqno];
      Uchar *seq = multiseq->sequence + seqstart;
      double lcb = 0;			// logcumback up to previous character
      int w = 0;			// width of valid word preceding current character
					// capped at bg_order
      for (pos=0; pos<seqlen; pos++) {
        int lcbp_index=0, a2n=1, offset=0;
	double log_pwa=0;
	Uint aindex = A2I(seq[pos]);	// current character
	if (aindex < alen) {		// not a SEPARATOR
	  for (j=0; j<=w; j++) {
	    if (j > 0) {
	      a2n *= alen;
	      offset += a2n;
	      lcbp_index *= alen;
	    }
	    lcbp_index += A2I(seq[pos-w+j]);
	  }
	  lcbp_index += offset;
	  log_pwa = lcbp[lcbp_index];
          if (++w > bg_order) w = bg_order;
        } else {			// SEPARATOR
          w = 0;
	  log_pwa = 0;
        }
	lcb = logcumback[seqstart+pos] = lcb + log_pwa;
      } // pos
    } // seqno
  } // dir

} // set_logcumback

//
// Shuffle the letters of each of the sequences
// preserving k-mer frequencies and the positions of
// the fixed character.
//
void shuffle_multiseq(
  Multiseq *multiseq,			// multiseq to shuffle
  int kmer,				// preserve frequencies of k-mers
  Uchar fixed				// the character to fix locations of
) {
  int seqno, i, j;

  // Get arrays of sequence starts and lengths and record the maximum length.
  Uint *seqstarts = (Uint *) malloc(multiseq->numofsequences * sizeof(Uint));
  Uint *seqlengths = (Uint *) malloc(multiseq->numofsequences * sizeof(Uint));
  Uint maxlength = 0;
  for (seqno=0; seqno < multiseq->numofsequences; seqno++) {
    // Get the a pointer to the sequence and its length.
    seqstarts[seqno] = (seqno == 0) ? 0 : *(PEEKARRAY(&multiseq->markpos, seqno-1, Uint)) + 1;
    Uint *markposptr = PEEKARRAY(&multiseq->markpos, seqno, Uint);
    Uint seqend = (markposptr == NULL) ? multiseq->totallength : *markposptr;
    seqlengths[seqno] = seqend - seqstarts[seqno];
    if (seqlengths[seqno] > maxlength) maxlength = seqlengths[seqno];
  }

  // Shuffle each sequence.
  Uchar *shuffled = (Uchar *) malloc(maxlength * sizeof(Uchar));
  for (seqno=0; seqno < multiseq->numofsequences; seqno++) {
    Uint seqstart = seqstarts[seqno];
    Uint seqlen = seqlengths[seqno];
    Uchar *original = multiseq->sequence + seqstart;
    ushuffle((char *) original, (char *) shuffled, seqlen, kmer);
    for (i=0, j=0; i<seqlen; i++) {
      // Advance to non-fixed character in shuffled seq.
      while (j<seqlen && shuffled[j]==fixed) j++;
      if (original[i] == fixed) {
        /* NOOP */ 			// fix the position of the fixed character
      } else {
        original[i] = shuffled[j++];	// copy shuffled character
      }
    }
  }

  // Free space.
  free(shuffled);
  free(seqstarts);
  free(seqlengths);
} // shuffle_multiseq

//
// Read in positive and negative sequences.
// Non-core letters are converted to SEPARATOR.
// Create negative sequences by shuffling in none given.
// Read or create the background model.
//
// Returns a Multiseq object containing:
//	Pos
//	Neg
//	[ RC(Pos)
//	  RC(Neg)
//      ]
//
Multiseq *read_pos_neg_seqs(
  STREME_OPTIONS_T *options,
  BOOL do_rc,				// append reverse complement of pos+neg sequences
					// negative sequences by shuffling positive sequences
  BOOL allow_ambigs,			// don't convert ambiguous characters to SEPARATOR
  BOOL no_trim,				// don't trim sequences to average length
  BOOL set_back,			// set the background in the multiseq objects
  ARRAY_T *background,                  // IN MEME-style background model
  Multiseq **test_multiseq_ptr		// OUT hold-out set
) {
  char *posfile = options->posfile;	// name of FASTA file with positive sequences
  char *negfile = options->negfile; 	// name of FASTA file with negative sequences or NULL
  ALPH_T *alph = options->alph;		// the sequence alphabet
  ALPHABET_T alphabet_type = options->alphabet_type;
  int order = options->order;		// shuffle using m-order shuffle
  double hofract = options->hofract;	// size of test set
  int minwidth = options->minwidth;	// mininum motif width
  int maxtotallength = options->totallength;	// truncate each each sequence set (+/-) to this length
  Multiseq *multiseq=NULL, *test_multiseq=NULL, *negmultiseq=NULL, *test_negmultiseq=NULL;

  // Convert DNA to RNA?
  BOOL is_rna = (alphabet_type == Rna || (alphabet_type == Custom && alph_check(alph, RNA)));
  if (is_rna) DEBUG_MSG(NORMAL_VERBOSE, "# NOTE: Will convert any DNA sequences to RNA.\n");

  // Read in the positive sequences.
  if (read_fasta_to_multiseqs(&multiseq, &test_multiseq, hofract, MIN_HO_SIZE, do_rc, is_rna, allow_ambigs, posfile, alph, COMPALPH, minwidth, maxtotallength) != 0) { 
    DEBUG_FMT(QUIET_VERBOSE, "ERROR: There was a problem reading the primary sequence file '%s'.\n", posfile); 
    exit(EXIT_FAILURE);
  }

  // Check that there are enough positive sequences.
  if (multiseq->numofsequences < MIN_SEQUENCES) {
    DEBUG_FMT(QUIET_VERBOSE, "ERROR: There are too few (valid) primary sequences (%d < %d).\n", 
      multiseq->numofsequences, MIN_SEQUENCES); 
    exit(EXIT_FAILURE);
  }

  // Check that the positive hold-out set was created if needed.
  if (hofract > 0 && test_multiseq == NULL) {
    DEBUG_FMT(QUIET_VERBOSE, 
      "# Warning: No hold-out set was created because the primary hold-out set \n"
      "#          would have had fewer than %d sequences.\n", MIN_HO_SIZE);
  }

  // Read in the negative sequences if needed.
  if (options->objfun == DE || options->objfun == NO_OBJFUN) {
    if (read_fasta_to_multiseqs(&negmultiseq, &test_negmultiseq, (test_multiseq==NULL ? 0 : hofract), MIN_HO_SIZE, do_rc, is_rna, allow_ambigs, negfile, alph, COMPALPH, minwidth, maxtotallength) != 0) { 
      DEBUG_FMT(QUIET_VERBOSE, "ERROR: There was a problem reading the control sequence file '%s'.\n", negfile); 
      exit(EXIT_FAILURE);
    }
    // Check that the negative hold-out set was created if needed and
    // redo without hold-out sets if not.
    if (hofract > 0 && test_multiseq != NULL && test_negmultiseq == NULL) {
      DEBUG_FMT(QUIET_VERBOSE, 
        "# Warning: No hold-out set was created because the control hold-out set\n"
	"#          would have had fewer than %d sequences.\n", MIN_HO_SIZE);
      // Reread both sets with hofract=0.
      freemultiseq(multiseq);
      freemultiseq(test_multiseq);
      freemultiseq(negmultiseq);
      (void) read_fasta_to_multiseqs(&multiseq, &test_multiseq, 0, MIN_HO_SIZE, do_rc, is_rna, allow_ambigs, posfile, alph, COMPALPH, minwidth, maxtotallength);
      (void) read_fasta_to_multiseqs(&negmultiseq, &test_negmultiseq, (test_multiseq==NULL ? 0 : 0), MIN_HO_SIZE, do_rc, is_rna, allow_ambigs, negfile, alph, COMPALPH, minwidth, maxtotallength);
    }
  }

  // Finish up the multiseq objects.
  DEBUG_FMT(NORMAL_VERBOSE, "# Positive sequences \"%s\" - training: %d hold-out: %d\n", 
    posfile, multiseq->numofsequences, test_multiseq ? test_multiseq->numofsequences : 0);
  multiseq->npos = multiseq->numofsequences;
  multiseq->pos_length = multiseq->totallength;
  if (negmultiseq) multiseq->nneg = negmultiseq->numofsequences;
  // Check that there are enough negative sequences.
  if (options->objfun != CD && multiseq->nneg < MIN_SEQUENCES) {
    DEBUG_FMT(QUIET_VERBOSE, "ERROR: There are too few (valid) control sequences (%d < %d).\n", 
      multiseq->nneg, MIN_SEQUENCES); 
    exit(EXIT_FAILURE);
  }
  if (test_multiseq) {
    test_multiseq->npos = test_multiseq->numofsequences;
    test_multiseq->pos_length = test_multiseq->totallength;
    if (test_negmultiseq) test_multiseq->nneg = test_negmultiseq->numofsequences;;
  }
  if (negmultiseq) {
    if (negfile != posfile) {
      DEBUG_FMT(NORMAL_VERBOSE, "# Negative sequences \"%s\" - training: %d hold-out: %d\n", 
	negfile, negmultiseq->numofsequences, test_negmultiseq ? test_negmultiseq->numofsequences : 0);
    } else {
      DEBUG_FMT(NORMAL_VERBOSE, "# Negative sequences are shuffled primary sequences (%d-order) - training: %d hold-out: %d\n",
	order, multiseq->numofsequences, test_multiseq ? test_multiseq->numofsequences : 0);
    }
  }
  
  // Shuffle the negative sequences if they are the positive sequences.
  if (negfile == posfile) {
    shuffle_multiseq(negmultiseq, order+1, SEPARATOR);
    if (test_multiseq) shuffle_multiseq(test_negmultiseq, order+1, SEPARATOR);
  }

  // Create background model from ALL the control sequences unless it was passed in.
  if (set_back) {
    set_multiseq_background(options, multiseq, alph, do_rc, background);
    if (test_multiseq) {
      test_multiseq->bg_order = multiseq->bg_order;
      test_multiseq->background = multiseq->background;
      test_multiseq->lcbp = multiseq->lcbp;
    } 
  }

  // Append the training negative sequences to the training positive sequences
  // trimming the negatives if they have longer average length.
  if (negmultiseq) {
    append_to_multiseq(multiseq, negmultiseq, no_trim, false);
    freemultiseq(negmultiseq);
    // Append the hold-out negative sequences to the hold-out training positive sequences.
    if (test_multiseq) {
      append_to_multiseq(test_multiseq, test_negmultiseq, no_trim, true);
      freemultiseq(test_negmultiseq);
    }
  } else {
    multiseq->nneg = 0;
    multiseq->neg_length = 0;
    if (test_multiseq) {
      test_multiseq->nneg = 0;
      test_multiseq->neg_length = 0;
    }
  }

  // Ignore <pvt> if there is no hold-out set for computing unbiased p-values.
  if (test_multiseq == NULL && options->nmotifs == 0) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Ignoring <thresh> (%g) and setting <nmotifs> to %d.\n", options->thresh, DEFAULT_NMOTIFS);
    options->nmotifs = DEFAULT_NMOTIFS;
    options->thresh = DEFAULT_EVT;
  }
  
  // Append the reverse complement sequences to the positive (and negative sequences).
  if (do_rc) {
    append_rc_to_multiseq(multiseq);
    multiseq->do_rc = true;
    if (test_multiseq) {
      append_rc_to_multiseq(test_multiseq);
      test_multiseq->do_rc = true;
    }
  }

  // Initialize the lookup table of starts, lengths and 
  // minimum and maximum lengths.
  Uint seqno;
  int i;
  multiseq->use_binomial = false;
  multiseq->bernoulli = -1;
  for (i=0; i<2; i++) {
    Multiseq *mseq = (i==0) ? multiseq : test_multiseq;
    if (! mseq) continue;
    mseq->min_poslen = mseq->totallength;
    mseq->max_poslen = 0;
    mseq->min_neglen = mseq->totallength;
    mseq->max_neglen = 0;
    mseq->seqstarts = (Uint *) malloc(mseq->numofsequences * sizeof(Uint));
    mseq->seqlengths = (Uint *) malloc(mseq->numofsequences * sizeof(Uint));
    for (seqno = 0; seqno < mseq->numofsequences; seqno++) {
      // Get the a pointer to the sequence and its length.
      mseq->seqstarts[seqno] = (seqno == 0) ? 0 : *(PEEKARRAY(&mseq->markpos, seqno-1, Uint)) + 1;
      Uint *markposptr = PEEKARRAY(&mseq->markpos, seqno, Uint);
      Uint seqend = (markposptr == NULL) ? mseq->totallength : *markposptr;
      Uint seqlen = seqend - mseq->seqstarts[seqno];
      mseq->seqlengths[seqno] = seqlen;
      if (seqno < mseq->npos) {
        if (seqlen < mseq->min_poslen) mseq->min_poslen = seqlen;
        if (seqlen > mseq->max_poslen) mseq->max_poslen = seqlen;
      } else {
        if (seqlen < mseq->min_neglen) mseq->min_neglen = seqlen;
        if (seqlen > mseq->max_neglen) mseq->max_neglen = seqlen;
      }
    }
    // Get the average lengths of positive and negative sequences.
    // Account for separators between sequences in computing average sequence length.
    mseq->avg_poslen = (mseq->pos_length - mseq->npos + 1.0) / mseq->npos;
    mseq->avg_neglen = mseq->nneg ? (mseq->neg_length - mseq->nneg + 1.0) / mseq->nneg : 0;

    // Use Binomial Test or Fisher Exact Test with DE objective function?
    // We use the Fisher Exact Test unless the average lengths of
    // the main positive and negative sequences differ by >= 0.01% of the length
    // of the main positives.
    if (options->objfun == DE) {
      if (mseq == multiseq) {
	double diff = fabs(mseq->avg_poslen - mseq->avg_neglen)/mseq->avg_poslen;
	if (diff >= 0.0001) {
          multiseq->use_binomial = true;
	  // Note: this is a compromise over all motif widths.
	  double avg_w = (options->minwidth + options->maxwidth) / 2.0;
	  double pos_sites = mseq->npos * MAX(1, (mseq->avg_poslen - avg_w + 1));
	  double neg_sites = mseq->nneg * MAX(1, (mseq->avg_neglen - avg_w + 1));
	  multiseq->bernoulli = pos_sites / (pos_sites + neg_sites);
        }
      } else {
        // Set the test sequences to use the same test as the main sequences.
        test_multiseq->use_binomial = multiseq->use_binomial;
	test_multiseq->bernoulli = multiseq->bernoulli;
      }
    }
  }
 
  // Check that sequence lengths are OK for objective function.
  if (options->objfun == CD && 
    (
      multiseq->min_poslen != multiseq->max_poslen ||
      (test_multiseq && test_multiseq->min_poslen != test_multiseq->max_poslen)
    )
  ) {
    DEBUG_MSG(QUIET_VERBOSE, "ERROR: All sequences must be the same length with --objfun cd.\n");
    exit(EXIT_FAILURE);
  }

  // Announce how p-values are computed.
  if (options->objfun == DE) {
    if (multiseq->use_binomial) {
      DEBUG_FMT(NORMAL_VERBOSE, 
        "# Using Binomial test for p-values because primary and control sequences\n"
        "# have different average lengths: %g vs. %g. Bernoulli = %f\n",
        multiseq->avg_poslen, multiseq->avg_neglen, multiseq->bernoulli);
      DEBUG_MSG(QUIET_VERBOSE, 
        "# Warning: p-values will be inaccurate if primary and control\n"
        "#          sequences have different length distributions.\n");
    } else {
      DEBUG_MSG(NORMAL_VERBOSE, "# Using Fisher Exact test for p-values.\n");
      if (test_multiseq && (test_multiseq->avg_poslen > test_multiseq->avg_neglen)) {
        DEBUG_FMT(QUIET_VERBOSE, 
          "# Warning: p-values may be too small because primary hold-out sequences\n"
          "#          are longer on average than control hold-out sequences: %g > %g.\n",
          test_multiseq->avg_poslen, test_multiseq->avg_neglen);
      }
    }
  } else if (options->objfun == CD) {
    DEBUG_MSG(NORMAL_VERBOSE, "# Using Cumulative Bates distribution for p-values.\n");
  }

  *test_multiseq_ptr = test_multiseq;
  return(multiseq);
} // read_pos_neg_seqs

//
// Score the given sequences using the model's PSPM.
// Find the optimum score threshold if requested.
// Set the enrichment p-value in the model if requested.
// Save the matches in the model if requested.
//
void score_model_pssm(
  STREME_OPTIONS_T *options,		// STREME options
  Multiseq *multiseq,			// the positive and negative sequences
  Model *model,				// the model
  BOOL find_score_threshold, 		// find the best score threshold and set in model
  BOOL set_pvalue,			// compute the enrichment p-value of the model
  SAVE_MATCHES_T save_matches, 		// ZOOPS: save best match in every sequence (positive or negative) for nesting
					// PASSING: save all matches with score >= model->score_threshold for erasing
					// NONE: don't save matches
  BOOL is_ho, 				// data is the hold-out set
  BOOL use_cache  			// use the p-value cache
) {
  Uint i, j, seqno, pos, rc_pos;
  OBJFUN_T objfun = options->objfun;
  double log_pvalue;
  double score_threshold = model->score_threshold;
  BOOL do_rc = multiseq->do_rc;
  double avg_poslen = multiseq->avg_poslen;
  double avg_neglen = multiseq->avg_neglen;
  Uint npos = multiseq->npos;
  Uint nneg = multiseq->nneg;
  Uint ntot = npos + nneg;
  BOOL use_binomial = multiseq->use_binomial;
  int bg_order = multiseq->bg_order;
  double *background = multiseq->background;
  double *lcbp = multiseq->lcbp;
  double *logcumback = multiseq->logcumback;
  Uint w = model->width;
  Uint alen = model->alen;
  // Sites for computing threshold and/or p-value.
  BOOL save_sites = (find_score_threshold || set_pvalue);
  Site *sites = save_sites ? (Site *) malloc((ntot+1) * sizeof(Site)) : NULL;
  Uint nsites = 0; 
  // Matches to be returned in model.
  Uint nmatches = 0;
  Site *matches = NULL;
  double pos_sites = npos * MAX(1, (avg_poslen - w + 1));
  double neg_sites = nneg * MAX(1, (avg_neglen - w + 1));
  double bernoulli = use_binomial ? pos_sites / (pos_sites + neg_sites) : -1;

  // Convert PSPM to PSSM.  
  INIT_PSSM_FROM_PROBS(model, background, bg_order);

  // Find the matching sites in each input sequence on either strand.
  // Seqno is always an input sequence, not a reverse complement.
  Uint rc_seqstart = 0;
  Uchar *rc_word = NULL;
  Uchar strand;
  for (seqno=0; seqno<ntot; seqno++) {
    Uint seqstart = multiseq->seqstarts[seqno];
    Uint seqlen = multiseq->seqlengths[seqno];
    if (seqlen < w) continue;
    if (do_rc) rc_seqstart = seqstart + multiseq->totallength/2 + 1;
    double best_score = model->min_possible_score - 1;	// make sure there is one site at least
    Uint best_pos = -1;
    Uchar best_strand = '+';
    for (pos=0; pos<seqlen-w+1; pos++) {
      Uchar *word = multiseq->sequence + seqstart + pos;
      if (do_rc) {
        rc_pos = (seqlen - pos) - w;
        rc_word = multiseq->sequence + rc_seqstart + rc_pos;
      }
      double score=0, rc_score=0;
      double score_adj=0, rc_score_adj=0;
      for (i=0; i<w; i++) {
        Uint aindex = A2I(word[i]);
        if (aindex == alen) break;		// Check for SEPARATOR.
        score += model->pssm[aindex][i];
        if (bg_order > 0 && !logcumback) SUBLCBP(word, i+1, alen, lcbp, bg_order, score_adj);
        if (do_rc) {
	  Uint rc_aindex = A2I(rc_word[i]);
          if (rc_aindex == alen) {		// Check for SEPARATOR.
            i = w - i - 1;
            break;
          }
	  rc_score += model->pssm[rc_aindex][i];
          if (bg_order > 0 && !logcumback) SUBLCBP(rc_word, i+1, alen, lcbp, bg_order, rc_score_adj);
        }
      } // w
      // See if there was a separator character in the word.
      if (i < w) {
        // Move position to the separator character and skip this site.
        pos += i;
        continue;
      } else {
        // Finish computing log likelihood ratio(s).
        if (bg_order > 0) {
          if (logcumback) {
            score_adj -= logcumback[seqstart+pos+w-1] - (pos > 0 ? logcumback[seqstart+pos-1] : 0);
            if (do_rc) rc_score_adj -= logcumback[rc_seqstart+rc_pos+w-1] - (rc_pos > 0 ? logcumback[rc_seqstart+rc_pos-1] : 0);
          }
          score += score_adj;
          if (do_rc) rc_score += rc_score_adj;
        }
        // Round the score(s) so sorting will be consistent across platforms.
        RND(score, RNDDIG, score);
        if (do_rc) RND(rc_score, RNDDIG, rc_score);
        // Save position if best for sequence so far.
        if (do_rc && rc_score > score) {
          score = rc_score;
          strand = '-';
        } else {
          strand = '+';
        }
        if (score > best_score) {
          best_pos = pos;
          best_score = score;
          best_strand = strand;
        }
        // Save the site if it is a MATCH and we are returning all PASSING matches.
	if (save_matches==PASSING && score >= score_threshold) {
          if (nmatches % RCHUNK == 0) Resize(matches, nmatches+RCHUNK, Site);
	  matches[nmatches].seqno = seqno;
	  matches[nmatches].pos = pos;
          matches[nmatches].is_positive = (seqno < npos);
          matches[nmatches].score = score;
          matches[nmatches].strand = strand;
	  nmatches++;
	}
      }
    } // pos
    // Record the BEST SITE in the current sequence UNLESS
    // there were no valid sites in it.
    if (save_sites && best_pos != -1) {
      sites[nsites].seqno = seqno;
      sites[nsites].pos = best_pos;
      sites[nsites].score = best_score;
      sites[nsites].strand = best_strand;
      sites[nsites].is_positive = (seqno < npos);
      sites[nsites].dtc = fabs(best_pos - (avg_poslen + w)/2.0) + 0.5;
      nsites++;
    }
    // Save the best site as a MATCH if we are returning just ZOOPS matches.
    if (save_matches==ZOOPS) {
      if (nmatches % RCHUNK == 0) Resize(matches, nmatches+RCHUNK, Site);
      matches[nmatches].seqno = seqno;
      matches[nmatches].pos = best_pos;
      matches[nmatches].score = best_score;
      matches[nmatches].strand = best_strand;
      matches[nmatches].is_positive = (seqno < npos);
      nmatches++;
    }
  } // seqno

  // Save the matches in the model if requested.
  if (save_matches==ZOOPS || save_matches==PASSING) {
    model->matches = matches;
    model->nmatches = nmatches;
  }

  // Done if not computing the threshold or p-value.
  if (!find_score_threshold && !set_pvalue) {
    return;
  }

  // Sort the sites by decreasing score
  // breaking ties by placing negative sites first
  // so that p-values will be conservative.
  qsort(sites, nsites, sizeof(Site), compare_site_score);
  
  // Find the best split point in terms of the objective function
  // or use the given score threshold.
  Uint pos_count=0, neg_count=0, best_pos_count=0, best_neg_count=0;
  Uint n_eff_tests=0;
  double dtc_sum=0, best_dtc_sum=0;
  double best_log_pvalue=0, best_score_threshold=0;
  for (i=0; i<nsites; i++) {
    // The score check below favors more specific motifs by requiring
    // that the log-odds score be greater than DEFAULT_OPTSCORE.
    if (sites[i].score <= DEFAULT_OPTSCORE) break;
    // See if we have passed the given score threshold.
    if (!find_score_threshold && sites[i].score < score_threshold) break;
    // Update the counts.
    if (sites[i].is_positive) {
      pos_count++;
      if (objfun==CD) dtc_sum += sites[i].dtc;
    } else {
      neg_count++;
    }
    // Compute the p-value and save the state if it is the best so far.
    if (find_score_threshold) {
      log_pvalue = 0;
      if (objfun == CD) {
	// CD:
	// The p-value can't be minimum unless the next site has higher dtc than the current mean.
	if (
	  i==nsites-1 || // last site
	  sites[i+1].dtc > dtc_sum/pos_count	// next site's dtc > mean dtc
	) {
	  GET_PVALUE(log_pvalue, objfun, use_binomial, pos_count, 0, 0, 0, 0.0, dtc_sum, w, avg_poslen, 
	    "score_model_pssm", options->pv_cache, options->cache_length);
	}
      } else if (objfun == DE) {
	// DE:
	// The p-value can't improve unless the site was a positive, and
	// no need to compute the p-value until the last in a run of positives.
	// Assume that sites are sorted so that runs of the same score are
	// sorted with negative sites first. That way p-value is only
	// computed at the end of a run of positives.
	if (sites[i].is_positive && 
	  (
	    i==nsites-1 || // last site
	    sites[i+1].is_positive == false || // end of positive run
	    sites[i+1].score <= 0 // next site's log-odds is not greater than 0
	  )
	) {
	  double **cache = (use_cache) ? options->pv_cache : NULL;
	  int cache_length = (use_cache) ? options->cache_length : 0;
	  GET_PVALUE(log_pvalue, objfun, use_binomial, pos_count, neg_count, npos, nneg, bernoulli, 0, 0, 0, 
	    "score_model_pssm", cache, cache_length);
	}
      } // p-value

      // Update score threshold if p-value improved.
      if (log_pvalue < best_log_pvalue) {
	best_score_threshold = sites[i].score;
	best_log_pvalue = log_pvalue; 
	best_pos_count = pos_count;
	best_neg_count = neg_count;
        best_dtc_sum = dtc_sum;
        n_eff_tests++;
      }
    } // find_score_threshold
  } // site

  // Set the best counts if we were not finding the score threshold.
  if (find_score_threshold) {
    // Nudge score down in case of roundoff later.
    model->score_threshold = best_score_threshold - FLOAT_EPS;
  } else {
    best_pos_count = pos_count;
    best_neg_count = neg_count;
    best_dtc_sum = dtc_sum;
  }

  // Compute the p-value.
  GET_PVALUE(log_pvalue, objfun, use_binomial, 
    best_pos_count, best_neg_count, npos, nneg, bernoulli,
    best_dtc_sum, w, avg_poslen, "score_model_pssm", 
    (use_cache ? options->pv_cache : NULL), 
    (use_cache) ? options->cache_length : 0);
  // Compute the DTC.
  double dtc = (objfun == CD && best_pos_count > 0) ? best_dtc_sum/best_pos_count : -1;
  // Compute the enrichment ratio.
  double enr_ratio = (objfun == DE) ? 
    ((best_pos_count+1.0)/(npos+1.0)) / ((best_neg_count+1.0)/(nneg+1.0)) :
    -1;

  //
  // Record the counts and the p-value in the model.
  //
  if (is_ho) {
    model->test_pos_count = best_pos_count;
    model->test_neg_count = best_neg_count;
    model->test_dtc = dtc;
    model->test_ratio = enr_ratio;
    model->test_log_pvalue = log_pvalue;
  } else {
    model->train_pos_count = best_pos_count;
    model->train_neg_count = best_neg_count;
    model->train_dtc = dtc;
    model->train_ratio = enr_ratio;
    model->train_log_pvalue = log_pvalue;
  }
  model->n_eff_tests = MAX(1, n_eff_tests);

  // Free space.
  free(sites);
  
} // score_model_pssm

//
// For sorting sites by decreasing score.
// Break ties by placing negative sites first.
//
int compare_site_score(
  const void *v1,
  const void *v2
)
{
  const Site *s1 = (const Site *) v1;
  const Site *s2 = (const Site *) v2;

  if (fabs(s1->score - s2->score) > FLOAT_EPS) {
    if (s1->score < s2->score) {
      return(+1);
    } else {
      return(-1);
    }
  } else if (s1->is_positive && !s2->is_positive) {
   return(+1);
  } else if (!s1->is_positive && s2->is_positive) {
   return(-1);
  } else if (s1->seqno > s2->seqno) {
    return(+1);
  } else if (s1->seqno < s2->seqno) {
    return(-1);
  } else if (s1->pos > s2->pos) {
    return(+1);
  } else {
    return(-1);
  }
} // compare_site_score

//
// For sorting sequences by decreasing score.
//
int compare_sequence_score(
  const void *v1,
  const void *v2
)
{
  const Passing_seq *s1 = (const Passing_seq *) v1;
  const Passing_seq *s2 = (const Passing_seq *) v2;

  if (fabs(s1->score - s2->score) > FLOAT_EPS) {
    if (s1->score < s2->score) {
      return(+1);
    } else {
      return(-1);
    }
  } else if (s1->dbno != s2->dbno) {
   return(s1->dbno - s2->dbno);
  } else {
   return(s1->seqno - s2->seqno);
  } 
} // compare_sequence_score

//
// Initialize the alphabet lookup table.
// Letters in alphabet have index in [0,alen).
// Non-core letters will have index = alen and value SEPARATOR.
// Returns the alphabet length.
//
Uint initialize_st_alphabet(
  ALPH_T *alph                  // the MEME-style alphabet
) {
  int i;
  int alen = alph_size_core(alph);
  for (i=0; i<MAX_ALENGTH; i++) {
   ALPHINDEX[i] = alen;
   COMPINDEX[i] = alen;
   ALPH[i] = COMPALPH[i] = SEPARATOR;
  }
  // Initialize the core symbols in the letter-to-index lookup table.
  for (i=0; i<alen; i++) {
   char c = ALPH[i] = alph_char(alph, i);
   ALPHINDEX[(Uint) c] = alph_index(alph, c);
  }
  // Initialize the core symbols in the complement letter-to-index lookup table.
  if (alph_has_complement(alph)) {
    for (i=0; i<alen; i++) {
      Uint c = ALPH[i];
      Uint ic = COMPINDEX[c] = alph_complement(alph, i);
      COMPALPH[c] = ALPH[ic];
    }
  }
  return(alen);
} // initialize_st_alphabet

//
// Append the reverse complement sequences to the end
// of the sequences in a multisequence object.
//
void append_rc_to_multiseq (Multiseq *multiseq)
{
  Uint i, j;
  Uint numofsequences = 2 * multiseq->numofsequences;
  Uint totallength = 2 * multiseq->totallength + 1;
  // Make space for the description pointers
  multiseq->startdesc = ALLOCSPACE_TLB(multiseq->startdesc,Uint,numofsequences+1);
  // Double the size of the sequence.
  multiseq->sequence = realloc(multiseq->sequence, totallength * sizeof(Uchar));
  Uchar *sequence = multiseq->sequence;
  // Append the reverse complement of the sequences.
  Uchar *original = sequence;
  Uchar *revcomp = sequence + multiseq->totallength;
  for (i=0; i<multiseq->numofsequences; i++) {
    if (i != 0) original++;			// skip separator
    Uint *markposptr = PEEKARRAY(&multiseq->markpos, i, Uint);
    Uint seqstart = original - sequence;
    Uint seqend = (markposptr == NULL) ? multiseq->totallength : *markposptr;
    Uint seqlen = seqend-seqstart;
    // Put in the separator and mark the position.
    PUSHARRAY (&multiseq->markpos, Uint, 128, (Uint) (revcomp - multiseq->sequence));
    *revcomp++ = SEPARATOR;
    for (j=0; j<seqlen; j++) {
      revcomp[seqlen-j-1] = COMP(original[j]);
    }
    original += seqlen;
    revcomp += seqlen;
    // Copy the pointer to the description of the sequence.
    multiseq->startdesc[multiseq->numofsequences + i] = multiseq->startdesc[i];
  }
  multiseq->totallength = totallength;
  multiseq->numofsequences *= 2;
} // append_rc_to_multiseq

//
// Append the negative sequences to the end
// of the positive sequences in a multisequence object.
// If the average length of the negative sequences is longer,
// they are centrally trimmed to have the same average length.
//
void append_to_multiseq(
  Multiseq *multiseq, 		// old sequences
  Multiseq *multiseq2,		// sequences to append
  BOOL no_trim,			// don't trim sequences to average length
  BOOL is_ho			// this is the hold-out set
) {
  int i;
  Uint numofsequences = multiseq->numofsequences + multiseq2->numofsequences;
  Uint totallength = multiseq->totallength + multiseq2->totallength + 1;
  Uint neglength;

  // Trim the second set of sequences if they are longer (on average).
  double avglen = (multiseq->totallength - multiseq->numofsequences + 1.0) / multiseq->numofsequences;
  double avglen2 = (multiseq2->totallength - multiseq2->numofsequences + 1.0) / multiseq2->numofsequences;
  double trim_fraction = no_trim ? 1 : ((avglen2 > avglen) ? avglen/avglen2 : 1);
  if (trim_fraction != 1) 
    DEBUG_FMT(NORMAL_VERBOSE, "# Trimming control %ssequences by %.2f%% to average primary sequence length (%.1f).\n", 
      is_ho ? "holdout " : "", 100*(1-trim_fraction), avglen);

  // Append the new sequence descriptions.
  Uint index_offset = TOPINDEXARRAY(&(multiseq->descspace), Uchar) + 1;
  Uint new_desc_len = TOPINDEXARRAY(&(multiseq2->descspace), Uchar) + 1;
  for (i=0; i<new_desc_len; i++) {
    Uchar *new_char_ptr = PEEKARRAY(&multiseq2->descspace, i, Uchar);
    PUSHARRAY(&multiseq->descspace, Uchar, 128, *new_char_ptr);
  }
  multiseq->startdesc = ALLOCSPACE_TLB(multiseq->startdesc,Uint,numofsequences+1);
  Uint *startdesc = multiseq->startdesc;
  startdesc += multiseq->numofsequences;
  for (i=0; i<=multiseq2->numofsequences; i++) {
    *startdesc++ = multiseq2->startdesc[i] + index_offset;
  }

  // Append the new sequences.
  // Increase the size of the sequence.
  Uchar *sequence = (Uchar *) malloc((totallength) * sizeof(Uchar));
  // Copy the old sequence into the new space.
  memcpy(sequence, multiseq->sequence, multiseq->totallength*sizeof(Uchar));
  free(multiseq->sequence);
  
  multiseq->sequence = sequence;
  sequence += multiseq->totallength;

  // Append a separator after the old sequences and mark it.
  *sequence++ = SEPARATOR;
  PUSHARRAY(&multiseq->markpos, Uint, 128, multiseq->totallength);
  
  // Append the new sequence to the new space, trimming if requested.
  if (trim_fraction == 1) {
    // Not trimming the new sequences.
    memcpy(sequence, multiseq2->sequence, multiseq2->totallength*sizeof(Uchar));
    // Append the markpos for each of the new sequences.
    for (i=0; i<multiseq2->numofsequences-1; i++) {
      PUSHARRAY(&multiseq->markpos, Uint, 128, multiseq->totallength + *(PEEKARRAY(&multiseq2->markpos, i, Uint)) + 1);
    }
    multiseq->neg_length = multiseq2->totallength;
  } else {
    // Trimming the new sequences.
    totallength = multiseq->totallength + 1;
    Uint seqno;
    for (seqno = 0; seqno < multiseq2->numofsequences; seqno++) {
      // Get the a pointer to the sequence and its length.
      Uint seqstart = (seqno == 0) ? 0 : *(PEEKARRAY(&multiseq2->markpos, seqno-1, Uint)) + 1;
      Uint *markposptr = PEEKARRAY(&multiseq2->markpos, seqno, Uint);
      Uint seqend = (markposptr == NULL) ? multiseq2->totallength : *markposptr;
      Uint seqlen = seqend - seqstart;
      Uint trimlen = MAX(1, trim_fraction * seqlen + 0.5);
      // Trim keeping the center of the new sequence and append it.
      Uint offset = (seqlen - trimlen)/2.0;
      memcpy(sequence, multiseq2->sequence + seqstart + offset, trimlen*sizeof(Uchar));
      sequence += trimlen;
      totallength += trimlen;
      // Add a separator and mark it unless this is the last sequence.
      if (seqno < multiseq2->numofsequences-1) {
	*sequence++ = SEPARATOR;
	PUSHARRAY(&multiseq->markpos, Uint, 128, totallength++);
      }
    }
    multiseq->neg_length = totallength - multiseq->totallength - 1;
  }

  multiseq->totallength = totallength;
  multiseq->numofsequences = numofsequences;
} // append_to_multiseq

//
//  Get a single-letter consensus from a STREME model.
//
char *get_single_letter_consensus(
  Model *model
) {
  int r, c;
  int w = model->width;
  int alength = model->alen;
  ALPH_T *alphabet = model->alph;

  // Create a MEME Suite motif from the STREME model.
  MATRIX_T *probs = allocate_matrix(w, alength);
  if (probs != NULL) {
    for(r=0; r<alength; r++)
      for(c=0; c<w; c++)
        set_matrix_cell(c, r, model->probs[r][c], probs);
  } else {
    DEBUG_MSG(QUIET_VERBOSE, "ERROR: Unable to allocate matrix get_single_letter_consensus.");
    exit(EXIT_FAILURE);
  }
  MOTIF_T *motif = allocate_motif("", "", alphabet, probs, NULL);

  // Get a string representing the motif consensus.
  STR_T *id_buf = str_create(10);
  str_clear(id_buf);
  motif2consensus(motif, id_buf, true);
  char *motif_id = str_internal(id_buf);
  str_destroy(id_buf, true);

  // Free space.
  free_matrix(probs);
  destroy_motif(motif);

  return(motif_id);
} // get_single_letter_consensus

//
// Read the contents of a description file.
//
char *get_description_file(
  char *filename			// name of description file
) {
  int i, j;
  FILE *dfile = NULL;
  char *description = NULL;

  if ((dfile = fopen(filename, "r")) == NULL) {
    DEBUG_FMT(QUIET_VERBOSE, "# Warning: Unable to open description file \"%s\" for reading.\n", filename);
    exit(EXIT_FAILURE);
  } else {
    char *text = malloc(502);			// Room for 500 characters.
    description = malloc(502);			// Room for 500 characters.
    int nread = fread(text, sizeof(char), 500, dfile);        // Read up to 500 char.
    char *c = description;
    // Convert to Unix newlines and condense multiple newlines.
    char old_c = '\0';
    int n_newlines = 0;
    for (i=0, j=0; i<nread; i++) {
      char c = text[i];
      old_c = (j==0) ? '\0' : description[j-1];
      // convert Microsoft \r\n to Unix \n.
      if (c == '\r') {
	description[j++] = '\n';		// convert '\r' to '\n'
        n_newlines++;
      } else if (c == '\n') {
        if (old_c != '\r') {			// condense \r\n
	  description[j++] = '\n';
	  n_newlines++;
        }
      } else {
	description[j++] = c;
        n_newlines = 0;
      }
      if (n_newlines > 2) j--;			// condense more than 2 newlines
    }
    // Remove trailing newline.
    if (j > 0 && description[j-1] == '\n') j--;
    description[j] = '\0';   // Null terminate string.
    free(text);
  }
  return(description);
} // get_description_file

//
// Get the sequences with sites passing the threshold and store in model.
//
void get_passing_sequences(
  STREME_OPTIONS_T *options,	// STREME options
  Model *model,			// the model, including its PASSING sites (see score_model_pssm)
  Multiseq *multiseq, 		// the sequences
  BOOL append,			// append passing sites to model if true
  BOOL sort			// sort the passing sequences by score if true
) {
  int i, j, seqno;

  // Get the best score for each sequence from the matches in the model.
  int npos = multiseq->npos;
  int nneg = multiseq->nneg;
  int nseqs = npos + nneg;
  double *best_score = malloc(nseqs * sizeof(double));
  double score_threshold = model->score_threshold;
  int npassing = 0;
  for (seqno=0; seqno<nseqs; seqno++) best_score[seqno] = score_threshold-1;
  for (i=0; i<model->nmatches; i++) {
    seqno = model->matches[i].seqno;
    double score = model->matches[i].score;
    if (score >= score_threshold && score > best_score[seqno]) {
      if (best_score[seqno] < score_threshold) npassing++;
      best_score[seqno] = score;
    }
  }

  // Add each sequence with a passing site to the list of passing sequences
  Passing_seq *sequences = append ? 
    (Passing_seq *) realloc(model->passing_seqs, (npassing+model->npassing) * sizeof(Passing_seq)) :
    (Passing_seq *) malloc(npassing * sizeof(Passing_seq));
  npassing = append ? model->npassing : 0;
  for (seqno=0; seqno<nseqs; seqno++) {
    if (best_score[seqno] >= score_threshold) {
      sequences[npassing].dbno = append ? 1 : 0;	// present to allow unique sort order
      sequences[npassing].seqno = seqno;		// present to allow unique sort order
      sequences[npassing].desc = DESCRIPTIONPTR(multiseq, seqno);
      sequences[npassing].score = best_score[seqno];
      sequences[npassing].is_tp = (seqno < npos);
      npassing++;
    }
  }

  // Resize the array before sorting.
  sequences = (Passing_seq *) realloc(sequences, npassing * sizeof(Passing_seq));

  // Sort the passing sequences.
  if (sort) qsort(sequences, npassing, sizeof(Passing_seq), compare_sequence_score);

  // Set the passing sequences in the model.
  model->passing_seqs = sequences;
  model->npassing = npassing;
} // get_passing_sequences

//
// Get the positional site distribution and site count histogram 
// in the positive sequences and save in model.
//
void get_site_distr(
  Model *model,			// The model, including its PASSING sites (see score_model_pssm)
  Multiseq *multiseq,		// The sequences
  SEQ_ALIGN_T align		// Align sequences left, center or right.
) {
  int i, offset=0;
  int w = model->width;		// width of motif
  int npos = multiseq->npos;	// number of positive sequences
  int maxlen = multiseq->max_poslen; 	// maximum length of positive sequences
  double score_threshold = model->score_threshold;

  // Create the array of best score for each seqno, set to score_threshold-1.
  double *best_score = (double *) malloc(npos * sizeof(double));
  for (i=0; i<npos; i++) best_score[i] = score_threshold - 1;

  // Create the array of the numbers of sites and numbers of
  // of best sites for each seqno, set to 0.
  int *n_sites = (int *) calloc(npos, sizeof(int));
  int *n_best_sites = (int *) calloc(npos, sizeof(int));

  // Set the best score and number of passing sites for each seqno.
  for (i=0; i<model->nmatches; i++) {
    // Sequence is positive?
    int seqno = model->matches[i].seqno;
    if (seqno >= npos) continue;
    // Site is passing?
    double score = model->matches[i].score;
    if (score < score_threshold) continue;		
    // Increment the number of passing sites in this sequence.
    n_sites[seqno]++;
    // Increment the number of (tied) best sites in this sequence.
    if (score > best_score[seqno]) {
      best_score[seqno] = score;
      n_best_sites[seqno] = 1;
    } else if (score == best_score[seqno]) {
      n_best_sites[seqno]++;
    }
  }

  // Create the array of positional site counts, set to zero.
  model->site_distr = (int *) calloc(maxlen, sizeof(int));

  // Update the positional site counts, spreading ties.
  for (i=0; i<model->nmatches; i++) {
    // Sequence is positive?
    int seqno = model->matches[i].seqno;
    if (seqno >= npos) continue;
    // Sequence has a passing site?
    if (n_best_sites[seqno] == 0) continue;
    // Site has the best score for sequence?
    double score = model->matches[i].score;
    if (score < best_score[seqno]) continue;
    int pos = model->matches[i].pos;
    int seqlen = multiseq->seqlengths[seqno];
    switch (align) {
      case LEFT:
	// align sequence left edges; position site at its left edge
	offset = 0;
	break;
      case CENTER:
	// align sequence centers; position site at its left edge
	offset = (maxlen - seqlen)/2.0;
	break;
      case RIGHT:
	// align sequence right edges; position site at its left edge
	offset = (maxlen - seqlen);
	break;
    }
    // Update the count of sites at this position.
    model->site_distr[pos+offset] += 1.0/n_best_sites[seqno];
  }

  // Get the total number of sequences with passing sites.
  model->total_sites = 0;
  for (i=0; i<npos; i++) if (n_best_sites[i] > 0) model->total_sites++;

  // Create the histogram of the number of sequences with each
  // number of matching sites.
  model->max_sites = 0;
  model->site_hist = (int *) calloc(maxlen+1, sizeof(int));
  for (i=0; i<npos; i++) { 
    int count = n_sites[i];
    if (count > 0) {
      model->site_hist[count]++;
      if (count > model->max_sites) model->max_sites = count;
    }
  }
  model->site_hist = realloc(model->site_hist, (model->max_sites + 1) * sizeof(int));

} // get_site_distr
