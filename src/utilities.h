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

#include <string.h>
#include "matrix.h"

struct function {
    char const* name;
    long double(*f)(long double, void const*);
    void const* arg;
};

struct sampled_function {
    char const* name;
    long double start;
    long double end;
    long double sampling_interval;
    size_t n_samples;
    long double const* samples;
};

struct result {
    long double value;
    long double error;
    unsigned long iterations;
    unsigned char convergence_rate;
};

struct interpolation {
    struct function const* function;
    char *name;
    long double start;
    long double end;
    size_t order;
    long double *coefficients;
    long double sampling_interval;
};

// struct linear_interpolation
// struct raised_cosine_interpolation
// struct least_squares_interpolation

struct sampled_function* sample_values(struct function const*, long double, long double, long double);
void destroy_sample(struct sampled_function*);

struct sampled_function* sample_derivative(struct sampled_function const*);

void report_result(struct result const*);

void gnuplot(char const * base, size_t n_functions, long double start, long double end, unsigned long points, ...);

struct interpolation* allocate_interpolation(struct function const*, long double start, long double end, size_t order);
void destroy_interpolation(struct interpolation*);
long double polynomial_value(long double x, struct interpolation const *interpolation);
long double polynomial_error(struct interpolation const*);
long double piecewise_linear_value(long double x, struct interpolation const *interpolation);
long double piecewise_linear_error(struct interpolation const*);
long double raised_cosine_value(long double x, struct interpolation const *interpolation);
long double raised_cosine_error(struct interpolation const*);