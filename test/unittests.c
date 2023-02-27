#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "test_utils/utest.h"

#include "ece5210.h"

#define NUM_TAPS 900
#define DECIMATE 100
#define INTERPOLATE 100

typedef struct array
{
    uint32_t size;
    int16_t data[];
} array;

typedef struct array_f32
{
    uint32_t size;
    float data[];
} array_f32;

extern float h_aa[];
extern float h_poly[];

static float interp_state[NUM_TAPS];
static float dec_state[NUM_TAPS];

double mean_absolute_error(array *y, array *x);
double mean_absolute_error_f32(array_f32 *y, array_f32 *x);
array * initialize_array(uint32_t size);
array_f32 * initialize_array_f32(uint32_t size);
void skip_lines(FILE *fp, int n_lines);
array * read_1d_mtx(char filename[]);
array_f32 * read_1d_f32_mtx(char filename[]);
array_f32 * rand_uniform_f32(float start, float end, uint32_t size);
float naive_interpolation(float x, float *h,
                          uint16_t L);
float naive_decimation(float x, float *h,
                       uint16_t M, uint16_t *count);


UTEST(ece5210_lab04, polyphase_components)
{
    init_firwin();

    int ret = system("python test_utils/ece5210_lab04_poly.py");
    char filename_h[] = "h_out.mtx";
    
    array_f32 *h_poly_sol = read_1d_f32_mtx(filename_h);

    double diff = 0.;
    
    for (uint16_t i = 0; i < h_poly_sol->size; i++)
    {

        diff += fabs((double)h_poly[i] -
                     (double)h_poly_sol->data[i]);
    }
    
    double error =  diff/h_poly_sol->size;

    ASSERT_LT(error, 0.001);

    ret = system("rm h_out.mtx");
    (void)ret;

    free(h_poly_sol);

    h_poly_sol = NULL;
    
}

UTEST(ece5210_lab04, decimate_accuracy)
{
    init_firwin();

    int ret = system("python test_utils/ece5210_lab04_dec.py");
    char filename_x[] = "x_out.mtx";
    char filename_y[] = "y_out.mtx";

    array_f32 *x = read_1d_f32_mtx(filename_x);
    array_f32 *y_sol = read_1d_f32_mtx(filename_y);

    array_f32 *y_my = initialize_array_f32(x->size);

    uint16_t M = (uint16_t)x->size/y_sol->size;

    
    // perform the filtering
    for (uint16_t i = 0; i < x->size; i++)
    {
        y_my->data[i] = polyphase_decimation(x->data[i],
                                             h_poly, M);

    }

    double diff = 0;

    for (uint16_t i = 0; i < y_sol->size; i++)
    {

        diff += fabs((double)y_my->data[M*i] -
                     (double)y_sol->data[i]);
    }

    
    double error =  diff/y_sol->size;

    ASSERT_LT(error, 0.001);

    ret = system("rm x_out.mtx y_out.mtx");
    (void)ret;

    free(x);
    free(y_sol);
    free(y_my);

    x = NULL;
    y_sol = NULL;
    y_my = NULL;
    
}

UTEST(ece5210_lab04, decimate_time)
{
    init_firwin();
    uint16_t count = 0;

    array_f32 * x = rand_uniform_f32(-15000, 15000, 1000000);

    clock_t start_time = clock();
    for (uint32_t i = 0; i < x->size; i++)
    {
        naive_decimation(x->data[i], h_aa, DECIMATE, &count);
    }
    clock_t end_time = clock();

    double naive_time = ((double)(end_time - start_time)/
                         CLOCKS_PER_SEC);
    
    count = 0;
    
    start_time = clock();
    for (uint32_t i = 0; i < x->size; i++)
    {
        polyphase_decimation(x->data[i], h_aa, DECIMATE);
    }
    end_time = clock();

    double poly_time = ((double)(end_time - start_time)/
                        CLOCKS_PER_SEC);

    ASSERT_LT(poly_time/naive_time,0.1);

   /* printf("%f %f\n",naive_time, poly_time); */

    free(x);
    x = NULL;
}


UTEST(ece5210_lab04, interpolate_accuracy)
{
    init_firwin();
    
    int ret = system("python test_utils/ece5210_lab04_int.py");
    char filename_x[] = "x_out.mtx";
    char filename_y[] = "y_out.mtx";

    array_f32 *x = read_1d_f32_mtx(filename_x);
    array_f32 *y_sol = read_1d_f32_mtx(filename_y);

    array_f32 *y_my = initialize_array_f32(y_sol->size);

    uint16_t L = (uint16_t)y_sol->size/x->size;

    // perform the filtering
    for (uint16_t i = 0; i < y_sol->size; i++)
    {
        if (i % L == 0)
        {
            y_my->data[i] = polyphase_interpolation(
                x->data[i/L],
                h_poly, L
                );            
        }
        else
        {
            y_my->data[i] = polyphase_interpolation(
                0,
                h_poly, L
                );            
        }
    }

    double diff = 0;

    for (uint16_t i = 0; i < y_sol->size; i++)
    {
        diff += fabs((double)y_my->data[i] -
                     (double)y_sol->data[i]);
    }
    
    double error =  diff/y_sol->size;

    ASSERT_LT(error, 0.003);

    ret = system("rm x_out.mtx y_out.mtx");
    (void)ret;

    free(x);
    free(y_sol);
    free(y_my);

    x = NULL;
    y_sol = NULL;
    y_my = NULL;
    
}


UTEST(ece5210_lab04, interpolate_time)
{
    init_firwin();
    
    array_f32 * x = rand_uniform_f32(-15000, 15000, 1000000);

    clock_t start_time = clock();
    for (uint32_t i = 0; i < x->size; i++)
    {
        naive_interpolation(x->data[i], h_aa, DECIMATE);
    }
    clock_t end_time = clock();

    double naive_time = ((double)(end_time - start_time)/
                         CLOCKS_PER_SEC);
    
    start_time = clock();
    for (uint32_t i = 0; i < x->size; i++)
    {
        polyphase_interpolation(x->data[i], h_aa, DECIMATE);
    }
    end_time = clock();

    double poly_time = ((double)(end_time - start_time)/
                        CLOCKS_PER_SEC);

    ASSERT_LT(poly_time/naive_time,0.1);
    
    free(x);
    x = NULL;
}

UTEST_MAIN();

/************************************* 
         SUPPORT FUNCTIONS 
*************************************/

double mean_absolute_error(array *y, array *x)
{
    double diff = 0;
    for (int n = 0; n < y->size; n++)
    {
        diff += fabs((double)y->data[n] - (double)x->data[n]);
    }

    return diff/y->size;
}

double mean_absolute_error_f32(array_f32 *y, array_f32 *x)
{
    double diff = 0;
    for (int n = 0; n < y->size; n++)
    {
        diff += fabs((double)y->data[n] - (double)x->data[n]);
    }

    return diff/y->size;
}


void skip_lines(FILE *fp, int n_lines)
{
    // Skip the first n lines
    for (int i = 0; i < n_lines; ++i)
    {
        if (fscanf(fp, "%*[^\n]\n") == -1)
        {
            printf("ERROR: fscanf() failed in %s on line %i\n",
                   __FILE__,__LINE__);
        }
    }
}

array * initialize_array(uint32_t size)
{
    array *out = malloc(sizeof(*out) +
                        sizeof(int16_t)*(size_t)size);

    out->size = size;
    return out;
}

array_f32 * initialize_array_f32(uint32_t size)
{
    array_f32 *out = malloc(sizeof(*out) +
                            sizeof(float)*(size_t)size);

    out->size = size;
    return out;
}



array * read_1d_mtx(char filename[])
{
    int temp;
    
    FILE *fp;
    fp = fopen(filename,"r");

    int n_lines = 2;
    skip_lines(fp, n_lines);


    // read in data and find max value
    if (fscanf(fp,"%i %*i", &temp) == -1)
    {
        printf("ERROR: fscanf in %s in line %i failed\n",
               __FILE__,__LINE__);
        exit(1);
    }

    array *out = initialize_array((uint32_t)temp);

    for (uint32_t m = 0; m < out->size; m++)
    {
        if ((fscanf(fp,"%i", &temp) == 1) &&
            !feof(fp))
        {
            out->data[m] = (int16_t)temp;		
        }
        else
        {
            printf("ERROR:  fscanf() failed\n");
        }
    }

    fclose(fp);
    
    return out;
}


array_f32 * read_1d_f32_mtx(char filename[])
{

    
    int size;
    float temp;
    
    FILE *fp;
    fp = fopen(filename,"r");

    int n_lines = 2;
    skip_lines(fp, n_lines);

    // read in the size
    if (fscanf(fp,"%i %*f", &size) == -1)
    {
        printf("ERROR: fscanf in %s in line %i failed\n",
               __FILE__,__LINE__);
        exit(1);
    }

    array_f32 *out = initialize_array_f32((uint32_t)size);

    for (uint32_t m = 0; m < out->size; m++)
    {
        if ((fscanf(fp,"%f", &temp) == 1) &&
            !feof(fp))
        {
            out->data[m] = temp;		
        }
        else
        {
            printf("ERROR:  fscanf() failed\n");
        }
    }

    fclose(fp);
    
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

float fir(float x, float *h, float *state, uint16_t size)
{
    float out = 0.f;

    for (uint16_t i = size - 1; i > 0; i--)
    {
        state[i] = state[i-1];
    }
    state[0] = x;

    
    for (uint16_t k = 0; k < size; k++)
    {
        out += h[k]*state[k];
    }

    return out;
}


float naive_interpolation(float x, float *h,
                          uint16_t L)
{
    float out = fir(x, h, interp_state, NUM_TAPS);
   
    return out;
}

float naive_decimation(float x, float *h, uint16_t M,
                       uint16_t *count)
{
    float out = fir(x, h, dec_state, NUM_TAPS);

    (*count)++;
    if (*count == M)
    {
        *count = 0;
    }
    
    if (*count == 1)
    {
        return out;
    }
    else
    {
        return 0; 
    }
}

