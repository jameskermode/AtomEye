/* atomeye C api to be called from Python */

#ifndef _ATOMEYELIB_H
#define _ATOMEYELIB_H

#define ATOMEYELIB_REDRAW 1
#define ATOMEYELIB_RUN_COMMAND 2
#define ATOMEYELIB_LOAD_ATOMS 3

extern void (*atomeyelib_on_click_atom)(int iw, int atom);
extern void (*atomeyelib_on_close)(int iw);
extern void (*atomeyelib_on_advance)(int iw, char *instr);
extern void (*atomeyelib_on_new)(int iw);

int atomeyelib_init(int argc, char* argv[], void* data);
void atomeyelib_set_handlers(void (*on_click)(int iw, int atom), 
			     void (*on_close)(int iw), 
			     void (*on_advance)(int iw, char *instr),
			     void (*on_new)(int iw));
int atomeyelib_open_window(int icopy);
void atomeyelib_set_title(int iw, char *title);
int atomeyelib_treatevent(int iw);
int atomeyelib_queueevent(int iw, int event, char *instr, void *data, char *outstr);


#endif
