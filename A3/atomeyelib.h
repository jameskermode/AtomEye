/* atomeye C api to be called from Python */

#ifndef _ATOMEYELIB_H
#define _ATOMEYELIB_H

#define ATOMEYELIB_REDRAW 1
#define ATOMEYELIB_RUN_COMMAND 2

#if HAVE_LIBATOMS
#include "libatoms.h"
#endif

typedef struct {
  fortran_t params[SIZEOF_FORTRAN_T];
  fortran_t properties[SIZEOF_FORTRAN_T];
  double lattice[3][3];
  int n_atom;
  int allocated;
  char title[255];
} Atomeyelib_atoms;

#define ATOMEYELIB_MAX_EVENTS 10
#define ATOMEYELIB_STR_LEN 255

typedef struct {
  int event;
  char instr[ATOMEYELIB_STR_LEN];
  Atomeyelib_atoms *data;
} Atomeyelib_events;

#define ATOMEYELIB_MAX_EVENT_ID 5


extern void (*atomeyelib_on_click_atom)(int mod_iw, int iw, int atom);
extern void (*atomeyelib_on_close)(int mod_id, int iw);
extern void (*atomeyelib_on_advance)(int mod_id, int iw, char *instr);
extern void (*atomeyelib_on_new)(int mod_id, int iw);
extern int (*atomeyelib_on_redraw)(int mod_id, int iw, Atomeyelib_atoms *atoms);

extern int atomeyelib_mod_id;

int atomeyelib_init(int argc, char* argv[], Atomeyelib_atoms* data);
void atomeyelib_set_handlers(void (*on_click)(int mod_id, int iw, int atom), 
			     void (*on_close)(int mod_id, int iw), 
			     void (*on_advance)(int mod_id, int iw, char *instr),
			     void (*on_new)(int mod_id, int iw),
			     int (*on_redraw)(int mod_id, int iw, Atomeyelib_atoms *atoms));
int atomeyelib_open_window(int mod_id, int icopy);
int atomeyelib_treatevent(int iw);
int atomeyelib_queueevent(int iw, int event, char *instr, Atomeyelib_atoms *data, char *outstr);
void atomeyelib_wait(int iw);
void atomeyelib_get_visible(int *n_shown, int **idx);

#endif
