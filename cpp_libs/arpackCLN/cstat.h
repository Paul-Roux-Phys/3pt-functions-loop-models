#ifndef TIMING_H
#define TIMING_H

struct timingstruct
{
  float t0,t1,t2,t3,t4,t5;
  
  int nopx, nbx, nrorth, nitref, nrstrt;
  float tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv;
  float tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv;
  float tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv;
  float tmvopx, tmvbx, tgetv0, titref, trvec;
};

extern struct timingstruct timing;

#endif