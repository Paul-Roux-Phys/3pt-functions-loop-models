#ifndef DEBUGHEADER_H
#define DEBUGHEADER_H
struct debugstruct
{
  int logfil, ndigit, mgetv0;
  int msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  int mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd;
  int mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
};

extern struct debugstruct debug;

#endif