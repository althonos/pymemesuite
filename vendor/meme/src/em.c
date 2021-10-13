/***********************************************************************
*                                                                      *
*       MEME                                                           *
*       Copyright 1994-2017, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
// 6-22-99 tlb; combine OOPS and ZOOPS e_step 
// 6-28-99 tlb; add weight on prior on nsites, wnsites 
// 5-17-99 tlb; fix bug--strcpy cons to model->cons 
// em.c 

#include "meme.h"
#include "psp.h"
 
static bool check_convergence(
  THETA old_theta, // before EM iteration 
  THETA new_theta, // after EM iteration 
  int w,           // width of motif 
  double distance, // convergence radius 
  int alength,     // alphabet length 
  int iter,        // current iteration number 
  int maxiter      // maximum iterations 
);

/**********************************************************************/
/*
        em

        Uses a version of the EM algorithm (expectation maximization).

*/
/**********************************************************************/
void em(
  MODEL *model,                 // the model 
  DATASET *dataset,             // the dataset 
  DATASET *neg_dataset,         // the control dataset
  S_POINT *s_point              // starting point in case EM wanders
)
{
  MOTYPE mtype = model->mtype;                 // type of model 
  int max_w = model->w;                        // width of motif 
  int alength = alph_size_wild(dataset->alph); // length of alphabet 
  PRIORS *priors = dataset->priors;            // the priors 
  double wnsites = dataset->wnsites;           // weight on prior on nsites 
  int maxiter = dataset->maxiter;              // maximum number of iterations 
  double distance = dataset->distance;         // stopping criterion 
  THETA theta_save;
  int iter;                                    // iteration number 
  double (*E_STEP)(MODEL *, DATASET *);        // expectation step 
  double (*E_STEP0)(MODEL *, DATASET *);       // expectation step 
  // maximization step function 
  void (*M_STEP)(MODEL *, DATASET *, PRIORS *, double);
  bool converged = false; // EM has converged 

  // create a place to save old value of theta 
  create_2array(theta_save, double, max_w, alength);

  // set up the correct type of EM to run 
  M_STEP = m_step;
  E_STEP = e_step;
  E_STEP0 = e_step;
  switch (mtype) {
    case Oops:
    case Zoops:
      E_STEP = e_step;
      break;
    case Tcm:
      E_STEP = tcm_e_step;
      break;
    default:
      fprintf(stderr, "Unknown model type in em()! \n");
      exit(1);
      break;
  }

  /* get the probability that a site starting at position x_ij would
     NOT overlap a previously found motif; used in E_STEP.
  */
  get_not_o(dataset, model->w);

  // renormalize the PSP to the current motif width.
  if (dataset->pspfile && model->mtype != Tcm) {
    psp_renormalize(dataset, model->w, model->invcomp, model->mtype);
    if (neg_dataset) { // set PSP to uniform for control dataset
      psp_renormalize(neg_dataset, model->w, model->invcomp, model->mtype);
    }
  }

  // Perform EM for number of iterations or until no improvement 
  char *cons;
  if (PRINT_LL) {
    cons = get_single_letter_consensus(model, dataset->alph);
    printf("iter=0 w=%d consensus=%s\n", model->w, cons);
    myfree(cons);
  }
  for (iter=0; iter < maxiter; iter++) {
    int w = model->w;           // width of motif 
    THETA theta = model->theta; // final theta of motif 

    if (PRINTALL) { printf("\niter %d\n", iter); }
    if ((!NO_STATUS) && ((iter % 10) == 0)) {
      fprintf(stderr, 
        "\rem: w=%4d, psites=%4.0f, iter=%4d ", w, model->psites, iter);
    }

    // save current contents of theta 
    copy_theta(theta, theta_save, w, alength);

    // expectation step 
    model->ll = E_STEP(model, dataset);
    if (PRINT_Z) print_zij(dataset, model);

    // maximization step 
    M_STEP(model, dataset, priors, wnsites);

    // print status if requested 
    if (PRINT_LL) {
      double nsites_obs = model->lambda_obs * wps(dataset, w, 1);
      double ll0 = model->mll_0;
      double ll1 = model->mll_1;
      double llr = ll1 - ll0;
      cons = get_single_letter_consensus(model, dataset->alph);
      printf("iter=%d w=%d llr=%8.2f nsites_obs=%6.1f consensus=%s\n",
        iter+1, model->w, llr, nsites_obs, cons);
      myfree(cons);
      printf("w %d ll1 = %f ll0 = %f\n", model->w, ll1, ll0);
    }
    if (PRINTALL) {
      int n = model->nsites_dis; // number of sites 
      printf("lambda= %8.6f (theta and obs)\n", model->lambda);
      printf("obs: \n");
      print_theta(0, NULL, 1, n, model->obs, model->w, 0, "", dataset, stdout);  
      printf("freq: \n");
      print_theta(0, NULL, 1, n, model->theta, model->w, 0, "", dataset, stdout);
    }

    // see if EM has converged 
    converged = check_convergence(theta_save, theta, w, distance, alength,
      iter, maxiter);

    if (converged) {iter++; break;} // done 
  }

  // save the number of iterations (counting from zero)
  model->iter += iter;

  // Discretize
  (void) discretize(model, dataset, neg_dataset, E_STEP);            

#ifdef FIXME
  // Check that EM didn't wander off; starting point should be among sites.
  // FIXME: Test if Classic would be improved by this step as well.
  if (dataset->objfun != Classic && s_point) {
    int iseq = s_point->iseq;
    int ioff = s_point->iseq;
    int imotif = model->imotif;
    char seq[MAXSITE+1];
    int i;
    for (i = 0; i < s_point->w0; i++) {
      seq[i] = alph_char(dataset->alph, s_point->e_cons0[i]);
    }
    seq[s_point->w0] = '\0';
    bool start_in_motif = false;
    bool cons_given = iseq < 0;		// true if -cons was given
    int w = cons_given ? strlen(s_point->cons0) : 0;
    for (i=0; i<model->nsites; i++) {
      if (cons_given) {
        // Check if start word matches a site.
        iseq = model->maxima[i].x;
        ioff = model->maxima[i].y;
//FIXME: handle wild cards
        if (!strncmp(s_point->cons0, dataset->samples[iseq]->seq+ioff, w)) {
	  start_in_motif = true;
	  break;
        }
      } else {
        // Check if start position is one of sites.
	if (model->maxima[i].x == iseq && model->maxima[i].y == ioff) {
	  start_in_motif = true;
	  break;
	}
      }
    }
    if (!start_in_motif) {		// EM wandered off)
      fprintf(stderr, "EM wandered from %s--using starting point as motif (w0 %d w %d)\n", seq, s_point->w0, model->w);
      // Reinitialize the model from the starting point.
      init_model(s_point, model, dataset, imotif);
      // Discretize the starting point.
      (void) discretize(model, dataset, neg_dataset, E_STEP);            
    }
  }
#endif

  /* Run 1 m_step, get relative entropy.  
     Include only group 0 primary samples in this final M_STEP if doing NZ.
     Otherwise, include all primary samples in this final M_STEP.
     Use b=0 if using MegaP heuristic.
  */
  if (priors->ptype == MegaP) SWAP(PriorLib *, priors->plib, priors->plib0);
  set_seq_groups_to_include(dataset, dataset->primary_groups.nsites);
  M_STEP(model, dataset, priors, wnsites);
  restore_seq_groups_to_include(dataset);
  if (priors->ptype == MegaP) SWAP(PriorLib *, priors->plib, priors->plib0);
  
  // get the consensus of the model 
  {
    THETA theta = model->theta;
    int w = model->w;
    char *cons;
    cons = get_consensus(theta, w, dataset, 1, MINCONS); 
    strcpy(model->cons, cons);
    myfree(cons);   
  }

  free_2array(theta_save, max_w);
} // em 

/**********************************************************************/
/*
        check_convergence
*/
/**********************************************************************/
static bool check_convergence(
  THETA old_theta, // before EM iteration 
  THETA new_theta, // after EM iteration 
  int w,           // width of motif 
  double distance, // convergence radius 
  int alength,     // alphabet length 
  int iter,        // current iteration number 
  int maxiter      // maximum iterations 
)
{
  int i, j;
  double euclid;   // distance between old_theta and new_theta 
  bool converged;

  // calculate the euclidean change in theta 
  euclid = 0;
  for(i=0; i<w; i++) {
    for(j=0; j<alength; j++) {
      double diff = theta_ref(old_theta, i, j) - theta_ref(new_theta, i, j);
      euclid += diff * diff;
    }
  }
  euclid = sqrt(euclid);
  if (PRINTALL || PRINT_LL) { printf(" d_theta = %f\n", euclid); }

  if (euclid < distance) { // converged? 
    if (TRACE) printf("Converged to motif (< %g change) after %d iterations\n",
      distance, iter+1);
    converged = true;
  } else if (maxiter > 1 && iter == maxiter - 1) {
    // Use fprintf to print from all nodes in parallel. 
    if (TRACE) {
      fprintf(stdout, "Failed to converge after %d iterations!\n", maxiter);
    }
    converged = false;
  } else {
    converged = false;
  }

  return converged;
}
