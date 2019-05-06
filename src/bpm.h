#ifndef BPM_H
#define BPM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"

/* Must be called before bpm_256!!!!  */
extern void set_broadcast_mask(void);

extern uint8_t bpm_256(const uint8_t* t,const uint8_t* p,int n,int m);
extern uint8_t bpm(const uint8_t* t,const uint8_t* p,int n,int m);


#endif
