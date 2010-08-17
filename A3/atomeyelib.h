/* atomeye C api to be called from Python */

#ifndef _ATOMEYELIB_H
#define _ATOMEYELIB_H

#define ATOMEYELIB_REDRAW 1
#define ATOMEYELIB_RUN_COMMAND 2
#define ATOMEYELIB_LOAD_ATOMS 3

#include "libatoms.h"

typedef struct {
  fortran_t params[12];
  fortran_t properties[12];
  double lattice[3][3];
  int n_atom;
} Atomeyelib_atoms;

#define ATOMEYELIB_MAX_EVENTS 10
#define ATOMEYELIB_STR_LEN 255

typedef struct {
  int event;
  char instr[ATOMEYELIB_STR_LEN];
  Atomeyelib_atoms *data;
} Atomeyelib_events;

#define ATOMEYELIB_MAX_EVENT_ID 5


extern void (*atomeyelib_on_click_atom)(int iw, int atom);
extern void (*atomeyelib_on_close)(int iw);
extern void (*atomeyelib_on_advance)(int iw, char *instr);
extern void (*atomeyelib_on_new)(int iw);

int atomeyelib_init(int argc, char* argv[], Atomeyelib_atoms* data);
void atomeyelib_set_handlers(void (*on_click)(int iw, int atom), 
			     void (*on_close)(int iw), 
			     void (*on_advance)(int iw, char *instr),
			     void (*on_new)(int iw));
int atomeyelib_open_window(int icopy);
void atomeyelib_set_title(int iw, char *title);
int atomeyelib_treatevent(int iw);
int atomeyelib_queueevent(int iw, int event, char *instr, Atomeyelib_atoms *data, char *outstr);
void atomeyelib_wait(int iw);


#endif
