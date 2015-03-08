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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include "utilities.h"
#include "project1.h"

#define LEAST_SQUARES_POINTS 524288.0L
#define SQUARE_ROOT_TOLERANCE 1E-7L

struct result bisection_method(struct function const* function, long double x0, long double x1, long double tolerance) {
    struct result result;
    result.iterations = 0L;
    if(x0 >= x1) {
        result.value = NAN;
        result.error = NAN;
        return result;
    }

    long double y0, y1, ym;

    y0 = function->f(x0, function->arg);
    y1 = function->f(x1, function->arg);
    /* Error is halved to report it in +/- form */
    result.error = (x1 - x0) / 2.0L;
    result.value = (x0 + x1) / 2.0L;

    if(y0 * y1 > 0) {
        result.value = NAN;
        result.error = NAN;
        return result;
    }

    /* Tolerence is halved because it is compared to the error in +/- form, also halved */
    tolerance /= 2.0L;

    while(result.error > tolerance) {
        result.error /= 2.0L;
        result.iterations++;
        if(y0 == 0.0L) {
            result.error = 0.0L;
            result.value = x0;
            break;
        } else if(y1 == 0.0L) {
            result.error = 0.0L;
            result.value = x1;
            break;
        } else {
            ym = function->f(result.value, function->arg);
            if(y0 * ym < 0.0L) {
                x1 = result.value;
                y1 = ym;
            } else {
                x0 = result.value;
                y0 = ym;
            }
            result.value = (x0 + x1) / 2.0L;
        }
    }
    result.convergence_rate = 1;
    return result;
}

struct result newtons_method(struct function const* function, struct function const* derivative, long double x0, unsigned long max_iterations, long double tolerance) {
    long double errors[] = {0.0L, 0.0L, 0.0L};
    struct result result;

    for(result.iterations = 0L; result.iterations != max_iterations; result.iterations++) {
        result.value = x0 - function->f(x0, function->arg) / derivative->f(x0, derivative->arg);
        if((result.error = fabsl((result.value - x0) / result.value)) < tolerance) {
            break;
        }
	errors[0] = errors[1];
	errors[1] = errors[2];
	errors[2] = result.error;
        x0 = result.value;
    }
    result.error /= 2.0L;
    result.convergence_rate = (result.iterations < 3)?NAN:roundl(logl(errors[2]/errors[1])/logl(errors[1]/errors[0]));
    return result;
}

struct result altered_newtons_method(struct function const* const function, struct function const* const derivative, struct function const* const secondderivative, long double x0,  unsigned long const max_iterations, long double const tolerance) {
    long double errors[] = {0.0L, 0.0L, 0.0L};
    struct result result;
    double long temp_f, temp_fd, temp_fdd;
    for(result.iterations = 0L; result.iterations != max_iterations; result.iterations++) {
        temp_f = function->f(x0, function->arg);
        temp_fd = derivative->f(x0, function->arg);
        temp_fdd = secondderivative->f(x0, function->arg);
        result.value = x0 - (temp_f * temp_fd) / (temp_fd * temp_fd - temp_f * temp_fdd);
        if((result.error = fabsl((result.value - x0) / result.value)) < tolerance) {
            break;
        }
	errors[0] = errors[1];
	errors[1] = errors[2];
	errors[2] = result.error;
        x0 = result.value;
    }
    result.error /= 2.0L;
    result.convergence_rate = (result.iterations < 3)?NAN:roundl(logl(errors[2]/errors[1])/logl(errors[1]/errors[0]));
    return result;
}

struct interpolation const* lagrange_interpolation(struct function const *const function, long double const x0, long double const x1, unsigned long const order) {
    if(x0 >= x1) {
        return NULL;
    }
    long double const sampling_interval = (x1 - x0) / ((long double)order);

    struct sampled_function * const sampled_function = sample_values(function, x0, x1, sampling_interval);
    if(sampled_function == NULL) {
        fprintf(stderr, "lagrange_interpolation(): Unable to take samples.\n");
        return NULL;
    }

    struct interpolation * const lagrange = allocate_interpolation(function, x0, x1, sampled_function->n_samples - 1);
    if(lagrange == NULL) {
        fprintf(stderr, "lagrange_interpolation(): Unable to allocate memory.\n");
        destroy_sample(sampled_function);
        return NULL;
    }

    if(sampled_function->name != NULL) {
        size_t const len = strlen(sampled_function->name) + 50L;
        lagrange->name = malloc(len * sizeof(char));
        if(lagrange->name != NULL) {
            snprintf(lagrange->name, len, "Lagrange Interpolation of %s (order %ld)", sampled_function->name, order);
        }
    }

    long double * const temp_coefficients = malloc(sizeof(long double) * (order + 1));
    if(temp_coefficients == NULL) {
        fprintf(stderr, "lagrange_interpolation(): Unable to allocate memory.\n");
        destroy_interpolation(lagrange);
        destroy_sample(sampled_function);
        return NULL;
    }
    long double * const coefficients = lagrange->coefficients;

    long double coefficient_denominator;
    size_t temp_poly_degree;

    /* Initialise coefficients */
    for(size_t i = 0; i <= order; i++) {
        coefficients[i] = 0.0L;
    }

    for(size_t i = 0; i <= order; i++) {
        coefficient_denominator = powl(sampling_interval, order);
        temp_coefficients[0] = sampled_function->samples[i];
        for(size_t j = 0; j <= order; j++) {
            /* Compute denominator */
            if(i != j) {
                coefficient_denominator *= (long double)i - (long double)j;
                temp_poly_degree = j - ((j > i) ? 1 : 0);
                /* Distribute multiplication to get coefficient numerators */
                long double xi = -(x0 + sampling_interval * (long double)j);
                temp_coefficients[temp_poly_degree + 1] = temp_coefficients[temp_poly_degree];
                for(size_t l = temp_poly_degree; l != 0; l--) {
                    temp_coefficients[l] *= xi;
                    temp_coefficients[l] += temp_coefficients[l - 1];
                }
                temp_coefficients[0] *= xi;
            }
        }
        /* Add coefficients to overall coefficients */
        for(size_t j = 0; j <= order; j++) {
            coefficients[j] += temp_coefficients[j] / coefficient_denominator;
        }
    }

    destroy_sample(sampled_function);
    free(temp_coefficients);

    return lagrange;
}

struct interpolation const* piecewise_linear_interpolation(struct function const *const function, long double const x0, long double const x1, unsigned long const order) {
    if(x0 >= x1) {
        return NULL;
    }
    long double const sampling_interval = (x1 - x0) / ((long double)order);

    struct sampled_function * const sampled_function = sample_values(function, x0, x1, sampling_interval);
    if(sampled_function == NULL) {
        fprintf(stderr, "piecewise_linear_interpolation(): Unable to take samples.\n");
        return NULL;
    }

    struct interpolation * const piecewise_linear = allocate_interpolation(function, x0, x1, sampled_function->n_samples - 1);
    if(piecewise_linear == NULL) {
        fprintf(stderr, "piecewise_linear_interpolation(): Unable to allocate memory.\n");
        destroy_sample(sampled_function);
        return NULL;
    }

    if(sampled_function->name != NULL) {
        size_t const len = strlen(sampled_function->name) + 60L;
        piecewise_linear->name = malloc(len * sizeof(char));
        if(piecewise_linear->name != NULL) {
            snprintf(piecewise_linear->name, len, "Piecewise Linear Interpolation of %s (order %ld)", sampled_function->name, order);
        }
    }

    piecewise_linear->sampling_interval = sampled_function->sampling_interval;
    long double * const coefficients = piecewise_linear->coefficients;

    /* Initialise coefficients */
    for(size_t i = 0; i < sampled_function->n_samples; i++) {
        coefficients[i] = sampled_function->samples[i] / sampled_function->sampling_interval;
    }

    destroy_sample(sampled_function);

    return piecewise_linear;
}

struct interpolation const* raised_cosine_interpolation(struct function const *const function, long double const x0, long double const x1, unsigned long const order) {
    if(x0 >= x1) {
        return NULL;
    }
    long double const sampling_interval = (x1 - x0) / ((long double)order);

    struct sampled_function * const sampled_function = sample_values(function, x0, x1, sampling_interval);
    if(sampled_function == NULL) {
        fprintf(stderr, "raised_cosine_interpolation(): Unable to take samples.\n");
        return NULL;
    }

    struct interpolation * const raised_cosine = allocate_interpolation(function, x0, x1, sampled_function->n_samples - 1);
    if(raised_cosine == NULL) {
        fprintf(stderr, "raised_cosine_interpolation(): Unable to allocate memory.\n");
        destroy_sample(sampled_function);
        return NULL;
    }

    if(sampled_function->name != NULL) {
        size_t const len = strlen(sampled_function->name) + 60L;
        raised_cosine->name = malloc(len * sizeof(char));
        if(raised_cosine->name != NULL) {
            snprintf(raised_cosine->name, len, "Raised Cosine Interpolation of %s (order %ld)", sampled_function->name, order);
        }
    }

    raised_cosine->sampling_interval = sampled_function->sampling_interval;
    long double * const coefficients = raised_cosine->coefficients;

    /* Initialise coefficients */
    for(size_t i = 0; i < sampled_function->n_samples; i++) {
        coefficients[i] = sampled_function->samples[i] / 2.0L;
    }

    destroy_sample(sampled_function);

    return raised_cosine;
}

struct interpolation const* least_squares_interpolation(struct function const *const function, long double const x0, long double const x1, unsigned long const order) {
    long double const sampling_interval = (x1 - x0) / (LEAST_SQUARES_POINTS);
    struct interpolation * least_squares = NULL;
    struct sampled_function * const sampled_function = sample_values(function, x0, x1, sampling_interval);
    if(sampled_function == NULL) {
        printf("error\n");
        goto error1;
    }
    struct matrix const sample_vector_transpose = {
        .elements = (long double*)sampled_function->samples,
        .cols = 1L,
        .rows = sampled_function->n_samples
    };
    struct matrix * const matrix_transpose = create_matrix(sampled_function->n_samples, order + 1L);
    if(matrix_transpose == NULL) {
        printf("error\n");
        goto error2;
    }
    matrix_set_row(matrix_transpose, 1L, 1.0L);
    long double * const xs = malloc(sizeof(long double)*sampled_function->n_samples);
    if(xs == NULL) {
        printf("error\n");
        goto error3;
    }
    for(size_t i = 0; i != sampled_function->n_samples; i++) {
        xs[i] = sampled_function->start + ((long double)i) * sampled_function->sampling_interval;
    }
    for(size_t i = 2; i <= (order + 1); i++) {
        matrix_set_row_vector_power(matrix_transpose, i, xs, (long double)i - 1.0L);
    }
    struct matrix * const matrix = transpose_matrix(matrix_transpose);
    if(matrix == NULL) {
        printf("error\n");
        goto error4;
    }
    struct matrix * const matrix_transpose_by_matrix_product = matrix_multiply(matrix_transpose, matrix);
    if(matrix_transpose_by_matrix_product == NULL) {
        printf("error\n");
        goto error5;
    }
    struct matrix * const matrix_transpose_by_matrix_product_inverse = matrix_inverse(matrix_transpose_by_matrix_product);
    if(matrix_transpose_by_matrix_product_inverse == NULL) {
        printf("error\n");
        goto error6;
    }
    struct matrix * const matrix_transpose_by_matrix_product_inverse_by_matrix_transpose = matrix_multiply(matrix_transpose_by_matrix_product_inverse, matrix_transpose);
    if(matrix_transpose_by_matrix_product_inverse_by_matrix_transpose == NULL) {
        printf("error\n");
        goto error7;
    }
    struct matrix *const matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose = matrix_multiply(matrix_transpose_by_matrix_product_inverse_by_matrix_transpose, &sample_vector_transpose);
    if(matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose == NULL || matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose->cols != 1L) {
        printf("error\n");
        goto error8;
    }
    least_squares = allocate_interpolation(function, x0, x1, order);
    if(least_squares == NULL) {
        goto error9;
    }
    if(sampled_function->name != NULL) {
        size_t const len = strlen(sampled_function->name) + 60L;
        least_squares->name = malloc(len * sizeof(char));
        if(least_squares->name != NULL) {
            snprintf(least_squares->name, len, "Least Squares Interpolation of %s (order %ld)", sampled_function->name, order);
        }
    }
    for(size_t i = 0L; i != matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose->rows; i++) {
        least_squares->coefficients[i] = matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose->elements[i];
    }

error9:
    destroy_matrix(matrix_transpose_by_matrix_product_inverse_by_matrix_transpose_by_sample_vector_transpose);
error8:
    destroy_matrix(matrix_transpose_by_matrix_product_inverse_by_matrix_transpose);
error7:
    destroy_matrix(matrix_transpose_by_matrix_product_inverse);
error6:
    destroy_matrix(matrix_transpose_by_matrix_product);
error5:
    destroy_matrix(matrix);
error4:
    free(xs);
error3:
    destroy_matrix(matrix_transpose);
error2:
    destroy_sample(sampled_function);
error1:
    return least_squares;
}

static long double square_root_helper(double long const x, double long *k)
{
    return x * x - *k;
}

static long double square_root_helper_derivative(double long const x, double long *k)
{
    return 2 * x;
}

struct result square_root_calculator(double long const k) {
    struct result result = {
        NAN,
        0,
        0
    };
    if(k < 0) {
        return result;
    } else if(k == 0.0L) {
        result.value = 0.0L;
        return result;
    } else if(k == 1.0L) {
        result.value = 1.0L;
        return result;
    }

    struct function const f[] = {
        {
            NULL,
            (long double(*)(long double, void const*))square_root_helper,
            &k
        },
        {
            NULL,
            (long double(*)(long double, void const*))square_root_helper_derivative,
            &k
        }
    };
    result = bisection_method(&f[0], 0.0L, k, k / 16);
    if(result.error != 0.0L) {
        unsigned long iterations = result.iterations;
        struct result newtons_result = newtons_method(&f[0], &f[1], result.value, 256, SQUARE_ROOT_TOLERANCE);
        newtons_result.iterations += iterations;
        return newtons_result;
    } else {
        return result;
    }
}

struct result adjusting_newtons_method(struct function const* function, struct function const* derivative, long double x0, unsigned long max_iterations, long double tolerance) {
    struct result result;
    char adjusting = 1;
    long double errors[] = {0.0L, 0.0L, 0.0L};
    long double m = 1.0L;
    for(result.iterations = 0L; result.iterations != max_iterations; result.iterations++) {
        result.value = x0 - m * function->f(x0, function->arg) / derivative->f(x0, derivative->arg);
        if((result.error = fabsl((result.value - x0) / result.value)) < tolerance) {
            break;
        }
	errors[0] = errors[1];
	errors[1] = errors[2];
	errors[2] = result.error;
	if (adjusting && result.iterations > 2 && (result.iterations % 3 == 0)) {
	  if (roundl(logl(errors[2]/errors[1])/logl(errors[1]/errors[0])) < 2.0L) {
	    m++;
	  } else {
	    adjusting = 0;
	  }
	}
        x0 = result.value;
    }
    result.error /= 2.0L;
    result.convergence_rate = (result.iterations < 3)?NAN:roundl(logl(errors[2]/errors[1])/logl(errors[1]/errors[0]));
    return result;
}
