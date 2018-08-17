
/*              Copyright 1993 Colorado State University     */
/*                      All Rights Reserved                  */

/* flow.h */
/* Header file for flow routines in century. */

#define LENFST 500  /* Length of flowstack */

extern struct stack{
  double *from;  /* Source */
  double *to;    /* Destination */
  double when;   /* Time to flow */
  double amt;    /* Amount */
} flowstack[LENFST+1];

extern int nflows;  /* Number of flows */

/* global function */

void flow_err(int error_num, double when);
