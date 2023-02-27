#ifndef ECE5210_H
#define ECE5210_H

#include <stdint.h>

#define PASSTHROUGH_RIGHT

void init_firwin(void);
float polyphase_decimation(float x, float *h, uint16_t M);
float polyphase_interpolation(float x, float *h, uint16_t L);
int16_t process_sample(int16_t sample_in);

#endif
