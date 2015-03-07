/*
 Copyright (c) 2015 Ricardo Iván Vieitez Parra

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include "utilities.h"

#define POLYNOMIAL_ERROR_POINT_MULTIPLIER 524288.0L

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795028841971L
#endif

/* Utility functions */
struct sampled_function* sample_values(struct function const * const function, long double const start, long double const end, long double const sampling_interval) {
    if(start >= end) {
        return NULL;
    }
    long double nl_samples = floorl((end - start) / sampling_interval) + 1.0L;
    size_t n_samples = (size_t)nl_samples;
    struct sampled_function * const sampled_function = malloc(sizeof(struct sampled_function));
    if(sampled_function == NULL) {
        fprintf(stderr, "Unable to allocate memory for sampling struct of %s\n", function->name);
        return NULL;
    }
    long double * const samples = malloc(sizeof(long double)*n_samples);
    if(samples == NULL) {
        free(sampled_function);
        fprintf(stderr, "Unable to allocate memory for %ld samples of %s\n", n_samples, function->name);
        return NULL;
    }
    if(function->name != NULL) {
        if((sampled_function->name = malloc(sizeof(char) * (strlen(function->name) + 1))) != NULL) {
            strcpy((char*) sampled_function->name, function->name);
        }
    } else {
        sampled_function->name = NULL;
    }
    sampled_function->start = start;
    sampled_function->end = start + nl_samples*sampling_interval;
    sampled_function->sampling_interval = sampling_interval;
    sampled_function->n_samples = n_samples;
    sampled_function->samples = samples;
    for(size_t i = 0; i != n_samples; i++) {
        samples[i] = function->f(start + (sampling_interval * ((long double)i)), function->arg);
    }
    return sampled_function;
}

void destroy_sample(struct sampled_function * const sample)
{
    free((void*)(sample->samples));
    if(sample->name != NULL) {
        free((void*)(sample->name));
    }
    free(sample);
}

struct sampled_function* sample_derivative(struct sampled_function const * const sampled_function) {
    struct sampled_function * const sampled_derivative = malloc(sizeof(struct sampled_function));
    char* name = NULL;
    if(sampled_function->n_samples < 2) {
        return NULL;
    }
    if(sampled_function->name != NULL) {
        size_t len = strlen(sampled_function->name) + 3;
        name = malloc(len * sizeof(char));
        snprintf(name, len, "(%s)\'", sampled_function->name);
    }
    size_t n_samples = sampled_function->n_samples - 1;
    long double sampling_interval = sampled_derivative->sampling_interval;
    if(sampled_derivative == NULL) {
        fprintf(stderr, "Unable to allocate memory for sampling struct of %s\n", name);
        return NULL;
    }
    long double * const samples = malloc(sizeof(long double) * (sampled_function->n_samples - 1));
    if(samples == NULL) {
        free(sampled_derivative);
        fprintf(stderr, "Unable to allocate memory for %ld samples of %s\n", n_samples, name);
        return NULL;
    }
    sampled_derivative->name = name;
    sampled_derivative->start = sampled_function->start;
    sampled_derivative->end = sampled_function->end - sampled_function->sampling_interval;
    sampled_derivative->sampling_interval = sampling_interval;
    sampled_derivative->n_samples = n_samples;
    sampled_derivative->samples = samples;
    for(size_t i = 0; i != n_samples; i++) {
        samples[i] = (sampled_function->samples[i + 1] - sampled_function->samples[i]) / sampling_interval;
    }
    return sampled_derivative;
}

// TODO: report rate of convergence
void report_result(struct result const * const result)
{
    size_t precision = 8;
    if(result->error != 0.0L && isfinite(result->error)) {
        precision = 1L + -1L * (size_t)(ceill(log10l(fabsl(result->error))));
    }
    printf("result: %.*Lf ± %.1LE, iterations: %ld\n", (int)precision, result->value, result->error, result->iterations);
}

void gnuplot(char const * const base, size_t const n_functions, long double const start, long double const end, unsigned long points, ...)
{
    char const **function_names = malloc(sizeof(char*)*n_functions);
    struct function const *t_f;
    char filename_buffer[1024];
    va_list ap;
    FILE *fp;
    va_start(ap, points);
    long double const sampling_interval = (end - start) / ((long double)points);
    size_t const base_len = strlen(base);
    memcpy(filename_buffer, base, sizeof(char)*base_len);
    memcpy(filename_buffer + base_len, "___", sizeof(char) * 3);
    for(size_t i = 0; i < n_functions; i++) {
        snprintf(filename_buffer + base_len + 3, 1024 - base_len - 3, "d%ld.csv", i);
        fp = fopen(filename_buffer, "w");
        t_f = va_arg(ap, struct function const *);
        function_names[i] = t_f->name;
        for(long double x = start; x < end; x += sampling_interval) {
            fprintf(fp, "%Lf,%Lf\n", x, t_f->f(x, t_f->arg));
        }
        fclose(fp);
    }
    memcpy(filename_buffer + base_len, ".gnuplot", sizeof(char) * 9);
    fp = fopen(filename_buffer, "w");
    fprintf(fp, "set datafile separator \",\";");
    fprintf(fp, "plot ");
    for(size_t i = 0; i < n_functions; i++) {
        fprintf(fp, "\"%s___d%ld.csv\" using 1:2 title \'%s\' with lines,", base, i, function_names[i]);
    }
    fprintf(fp, "\n");
    fclose(fp);
    va_end(ap);
    free(function_names);
}

static long double function_error(struct function const* const function1, struct function const* const function2, long double const start, long double const end, unsigned long const points)
{
    if(end <= start) {
        return NAN;
    }
    long double sampling_interval = (end - start) / (long double)points;
    struct sampled_function* const sampled_function1 = sample_values(function1, start, end, sampling_interval);
    if(sampled_function1 == NULL) {
        fprintf(stderr, "function_error(): Unable to take samples.\n");
        return NAN;
    }
    struct sampled_function* const sampled_function2 = sample_values(function2, start, end, sampling_interval);
    if(sampled_function2 == NULL) {
        destroy_sample(sampled_function1);
        fprintf(stderr, "function_error(): Unable to take samples.\n");
        return NAN;
    }
    /* Let's be on the safe side */
    size_t const actual_points = (sampled_function1->n_samples < sampled_function2->n_samples) ? sampled_function1->n_samples : sampled_function2->n_samples;
    long double difference2 = 0.0L;
    long double f2 = 0.0L;
    long double error;

    for(size_t i = 0; i < actual_points; i++) {
        difference2 += powl(sampled_function1->samples[i] - sampled_function2->samples[i], 2.0L);
        f2 += powl(sampled_function1->samples[i], 2.0L);
    }
    error = sqrtl(difference2 / f2);
    destroy_sample(sampled_function1);
    destroy_sample(sampled_function2);
    return error;
}

struct interpolation* allocate_interpolation(struct function const * const function, long double const start, long double const end, size_t const order) {
    struct interpolation *interpolation = malloc(sizeof(struct interpolation));
    if(interpolation != NULL) {
        interpolation->coefficients = malloc(sizeof(long double) * (order + 1));
        if(interpolation->coefficients == NULL) {
            free(interpolation);
            return NULL;
        }
        interpolation->function = function;
        interpolation->start = start;
        interpolation->end = end;
        interpolation->order = order;
        interpolation->name = NULL;
    }
    return interpolation;
}

void destroy_interpolation(struct interpolation* const interpolation)
{
    free(interpolation->name);
    free(interpolation->coefficients);
    free(interpolation);
}

long double polynomial_value(long double const x, struct interpolation const * const interpolation)
{
    long double y = interpolation->coefficients[0];
    long double xn = x;
    for(size_t i = 1; i <= interpolation->order; i++) {
        y += interpolation->coefficients[i] * xn;
        xn *= x;
    }
    return y;
}

long double polynomial_error(struct interpolation const * const interpolation)
{
    struct function const function2 = {
        interpolation->name,
        (long double(*)(long double, void const*))polynomial_value,
        interpolation
    };
    return function_error(interpolation->function, &function2, interpolation->start, interpolation->end, interpolation->order * POLYNOMIAL_ERROR_POINT_MULTIPLIER + 1);
}

long double piecewise_linear_value(long double const x, struct interpolation const * const interpolation)
{
    if(x < interpolation->start || x > interpolation->end) {
        return NAN;
    }
    long double const sampling_interval = interpolation->sampling_interval;
    long double const x_sample = x - interpolation->start;
    size_t const x0 = (size_t)floorl(x_sample / sampling_interval);
    size_t const x1 = (size_t)ceill(x_sample / sampling_interval);
    long double const a = sampling_interval * ((long double)x0 + 1.0L) - x_sample, b = sampling_interval * ((long double)x0) - x_sample;
    return interpolation->coefficients[x0] * a - interpolation->coefficients[x1] * b;
}

long double piecewise_linear_error(struct interpolation const * const interpolation)
{
    struct function const function2 = {
        interpolation->name,
        (long double(*)(long double, void const*))piecewise_linear_value,
        interpolation
    };
    return function_error(interpolation->function, &function2, interpolation->start, interpolation->end, interpolation->order * POLYNOMIAL_ERROR_POINT_MULTIPLIER + 1);
}

long double raised_cosine_value(long double const x, struct interpolation const * const interpolation)
{
    if(x < interpolation->start || x > interpolation->end) {
        return NAN;
    }
    long double const sampling_interval = interpolation->sampling_interval;
    long double const x_sample = x - interpolation->start;
    size_t const x0 = (size_t)floorl(x_sample / sampling_interval);
    size_t const x1 = (size_t)ceill(x_sample / sampling_interval);
    long double const a = x_sample / sampling_interval - ((long double)x0);
    long double const b = x_sample / sampling_interval - ((long double)x0 + 1.0L);
    return interpolation->coefficients[x0] * (1.0L + cosl(M_PI * a)) + interpolation->coefficients[x1] * (1.0L + cosl(M_PI * b));
}

long double raised_cosine_error(struct interpolation const * const interpolation)
{
    struct function const function2 = {
        interpolation->name,
        (long double(*)(long double, void const*))raised_cosine_value,
        interpolation
    };
    return function_error(interpolation->function, &function2, interpolation->start, interpolation->end, interpolation->order * POLYNOMIAL_ERROR_POINT_MULTIPLIER + 1);
}
