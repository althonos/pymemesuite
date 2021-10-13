/********************************************************************
 * FILE: ame.h
 * AUTHOR: Robert McLeay
 * CREATE DATE: 19/08/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * rsClover is a yet unpublished algorithm that seeks to assist in
 * determining whether a given transcription factor regulates a set
 * of genes that have been found by an experimental method that
 * provides ranks, such as microarray or ChIp-chip.
 *
 *
 * rsClover is a code name. The name will change prior to publication.
 ********************************************************************/
#ifndef __AME_H__
#define __AME_H__

#include <stdbool.h>
#include "html-monolith.h"
#include "motif.h"
#include "pssm.h"
#include "seq.h"

/*
 * ame-specific macros
 */

#define AVG_ODDS 0
#define MAX_ODDS 1
#define SUM_ODDS 2
#define TOTAL_HITS 3

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a>b)?a:b

#define MAX_SEQ_LENGTH 1e6

/*
 * Define ame constants
 */

#define MEME_FORMAT 1   // input format is meme

#define UNIFORM_BG 0    //uniform.
#define MOTIF_BG 1      // background frequencies taken from motif file
#define FILE_BG 2       // background frequencies taken from specified file

#define QUICK_RS 0      //use my bodgy way of doing the ranksum
#define BETTER_RS 1     //use fabian's correct way of doing the ranksum

#define POS_FASTA 0  //use the fluorescence sort as positive indicator
#define POS_PWM 1 //use the PWM sort as positive indicator

#define RANKSUM_METHOD 0        //use the ranksum test
#define FISHER_METHOD  1        //use Fisher's exact test
#define MULTIHG_METHOD  2       //use Fisher's exact test modified to label with 0,1, or 2.
#define LONG_MULTIHG_METHOD  3  //use Fisher's exact test modified to label with 0,1,2, or 3.
#define PEARSON_METHOD 4        //use the significance of the pearson correlation coefficient
#define SPEARMAN_METHOD 5       //use the significance of the Spearman rank correlation coefficient

#else
#endif
