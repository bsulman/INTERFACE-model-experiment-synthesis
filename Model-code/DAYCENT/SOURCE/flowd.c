
/*              Copyright 1993 Colorado State University                    */
/*                      All Rights Reserved                                 */

/*****************************************************************************
**
**  FILE:      flowd.c
**
**  PURPOSE:   Global data file for flow routines in century.
**
**  GLOBAL VARIABLES:
**    flowstack[] - variables used in Century flow routines
**    LENFST      - maximum number of elements in flowstack array
**    nflows      - indicates number of unflowed events stored in flowstack
**
*****************************************************************************/

#define LENFST 500  /* Length of flowstack */

struct stack{
  double *from;      /* Source */
  double *to;        /* Destination */
  double when;       /* Time to flow */
  double amt;        /* Amount */
} flowstack[LENFST+1];

int nflows;         /* Number of flows */
