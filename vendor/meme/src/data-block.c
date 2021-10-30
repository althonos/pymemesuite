/******************************************************************************
 * FILE: data-block.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-16
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the implementation of data structures and functions used 
 * to store and access a block of data. The meaning of a "block of data" is
 * determined by the constructor function, and the implementation of the 
 * functions providing the data_block_reader interface. Current implementations
 * include priors from MEME PSP files and sequence segments from FASTA files.
 *****************************************************************************/

#include <stdlib.h>
#include "data-block.h"
#include "utils.h"

struct data_block {
  size_t position;
  size_t block_size;
  size_t num_read;
  void *data;
};

/******************************************************************************
 * This function creates and initializes a data block structure
 * which will contain one prior as a double.
 *****************************************************************************/
DATA_BLOCK_T *new_prior_block() {
  DATA_BLOCK_T *data_block = NULL;
  data_block = mm_malloc(sizeof(DATA_BLOCK_T) * 1);
  data_block->position = 0;
  data_block->block_size = 1;
  data_block->num_read = 0;
  double *prior = mm_malloc(sizeof(double) * 1);
  data_block->data = prior;
  return data_block;
}

/******************************************************************************
 * This function creates and initializes a data block structure
 * which will contain a sequence segment of the given length.
 *****************************************************************************/
DATA_BLOCK_T *new_sequence_block(size_t block_size) {
  DATA_BLOCK_T *data_block = NULL;
  data_block = mm_malloc(sizeof(DATA_BLOCK_T) * 1);
  data_block->position = 0;
  data_block->block_size = block_size;
  data_block->num_read = 0;
  data_block->data = NULL;
  return data_block;
}

/******************************************************************************
 * This function frees a data block structure
 *****************************************************************************/
void free_data_block(DATA_BLOCK_T *data_block) {
  data_block->data = NULL;
  data_block->block_size = 0;
  data_block->num_read = 0;
  data_block->position = 0;
  myfree(data_block);
}

/******************************************************************************
 * This function returns the coordinate of the first data element in 
 * this block of data.
 *****************************************************************************/
size_t get_start_pos_for_data_block(DATA_BLOCK_T *data_block) {
  return data_block->position;
}

/******************************************************************************
 * This function sets the coordinate of the first prior in 
 * this block of priors.
 *****************************************************************************/
void set_start_pos_for_data_block(DATA_BLOCK_T *data_block, size_t start_pos) {
  data_block->position = start_pos;
}

/******************************************************************************
 * This function returns the number of data items allowed in this data block
 *****************************************************************************/
size_t get_block_size_from_data_block(DATA_BLOCK_T *data_block) {
  return data_block->block_size;
}

/******************************************************************************
 * This function returns the number of data items obtained in the last read
 *****************************************************************************/
size_t get_num_read_into_data_block(DATA_BLOCK_T *data_block) {
  return data_block->num_read;
}

/******************************************************************************
 * This function the number of data items read in the last read
 *****************************************************************************/
void set_num_read_into_data_block(DATA_BLOCK_T *data_block, size_t num_read) {
  data_block->num_read = num_read;
}

/******************************************************************************
 * This function returns an array of characters. The size of the array
 * can be found by calling get_num_read_into_data_block()
 *
 * The user should not free this array. It will be freed when the data reader
 * is freed.
 *****************************************************************************/
char *get_sequence_from_data_block(DATA_BLOCK_T *data_block) {
  return (char *) data_block->data;
}

/******************************************************************************
  * This function sets the sequence string in a data_block
 *****************************************************************************/
void set_sequence_in_data_block(DATA_BLOCK_T *data_block, char* seq) {
  data_block->data = seq;
}

/******************************************************************************
 * This function returns a single prior from a data block.
 *****************************************************************************/
void set_prior_in_data_block(DATA_BLOCK_T *data_block, double prior) {
  *((double *) data_block->data) = prior;
}

/******************************************************************************
 * This function set a single prior from a data block.
 *****************************************************************************/
double get_prior_from_data_block(DATA_BLOCK_T *data_block) {
  return *((double *) data_block->data);
}
