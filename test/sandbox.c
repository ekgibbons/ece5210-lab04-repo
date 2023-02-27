#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "ece5210.h"

#define DECIMATE 10

typedef struct array_f32
{
    uint32_t size;
    float data[];
} array_f32;

extern float h_aa[];


array_f32 * initialize_array_f32(uint32_t size)
{
    array_f32 *out = malloc(sizeof(*out) +
                            sizeof(float)*(size_t)size);

    out->size = size;
    return out;
}


array_f32 * rand_uniform_f32(float start, float end, uint32_t size)
{
    
    array_f32 *out = initialize_array_f32((uint32_t)size);

    srand((unsigned int)time(NULL));
    
    float a = end - start;
    float offset = a/2;
    for (uint32_t i = 0; i < out->size; i++)
    {
        out->data[i] = (float)rand()/(float)(RAND_MAX/a) - offset;
    }

    return out;
}


int main(int argc, char *argv[])
{

    init_firwin();

    array_f32 * x = rand_uniform_f32(-15000, 15000, 1000000);

    clock_t start_time = clock();
    for (uint32_t i = 0; i < x->size; i++)
    {
        polyphase_decimation(x->data[i], h_aa, DECIMATE);
    }
    clock_t end_time = clock();

    double poly_time = ((double)(end_time - start_time)/
                        CLOCKS_PER_SEC);

    printf("%f\n",poly_time);

    free(x);
    x = NULL;    
    return 0;
}
