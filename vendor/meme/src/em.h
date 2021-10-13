/*
 * $Id: em.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:37:45  nadya
 * Initial revision
 *
 */

#ifndef EM_H
#  define EM_H

/* 6-28-99 tlb; remove zoops and mcm e_steps */
/* 6-28-99 tlb; add wnsitesw */

extern void em(
  MODEL *model,			/* the model */
  DATASET *dataset, 		/* the dataset */
  DATASET *neg_dataset,		/* the control dataset */
  S_POINT *s_point		// starting point in case EM wanders
);
extern void m_step(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  PRIORS *priors,		/* the priors */
  double wnsitesw                /* weight on prior on nsites */
);
extern double e_step(
  MODEL *model,                 		/* the model */
  DATASET *dataset              		/* the dataset */
);
extern double tcm_e_step(
  MODEL *model,					/* the model */
  DATASET *dataset				/* the dataset */
);
extern double like_e_step(
  MODEL *model,                 	 	/* the model */
  DATASET *dataset              		/* the dataset */
);
extern void discretize(
  MODEL *model,					/* the model */
  DATASET *dataset,				/* the dataset */
  DATASET *neg_dataset,				/* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *)		// E_STEP function
);
extern void set_z (
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
);
extern void convert_theta_to_log(
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
);
int estep_maxima_nsites(
  MODEL *model,                         /* the model */
  DATASET *primary,                     /* the primary dataset */
  DATASET *control,                     /* the control dataset */
  double (*E_STEP)(MODEL *, DATASET *), // E_STEP function
  bool primary_groups[3],               // which primary groups to include
  bool control_groups[3],               // which control groups to include
  int min_nsites,                       /* minimum nsites */
  int max_nsites,                       /* maximum nsites */
  double thresh_in,                     // restrict to maxima
                                        // with prob >= thresh_in
                                        // set to -BIG to ignore
  double *log_pv,                       /* log p-value of score */
  double *llr,                          // LLR of best p-value
  double *thresh_out                    // optimal site prob threshold
);
#endif
