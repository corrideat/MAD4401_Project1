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

#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

struct matrix* create_matrix(size_t const cols, size_t const rows) {
    if(cols == 0 || rows == 0) {
        return NULL;
    }
    struct matrix * const matrix = malloc(sizeof(struct matrix));
    if(matrix != NULL) {
        matrix->elements = malloc(sizeof(long double) * rows * cols);
        if(matrix->elements == NULL) {
            free(matrix);
            return NULL;
        }
        matrix->rows = rows;
        matrix->cols = cols;
    }
    return matrix;
}

void destroy_matrix(struct matrix * const matrix)
{
    free(matrix->elements);
    free(matrix);
}

void matrix_set_row(struct matrix * const matrix, size_t const row, long double const value)
{
    if(row == 0 || row > matrix->rows) {
        return;
    }
    long double * const rowp = matrix->elements + (row - 1L) * matrix->cols;
    for(size_t i = 0; i != matrix->cols; i++) {
        rowp[i] = value;
    }
}

void matrix_set_row_vector(struct matrix * const matrix, size_t const row, long double const* const vector)
{
    if(row == 0 || row > matrix->rows) {
        return;
    }
    long double * const rowp = matrix->elements + (row - 1L) * matrix->cols;
    for(size_t i = 0; i != matrix->cols; i++) {
        rowp[i] = vector[i];
    }
}

void matrix_set_row_vector_power(struct matrix* const matrix, size_t const row, long double const* const vector, long double const power)
{
    if(row == 0 || row > matrix->rows) {
        return;
    }
    long double * const rowp = matrix->elements + (row - 1L) * matrix->cols;
    for(size_t i = 0; i != matrix->cols; i++) {
        rowp[i] = powl(vector[i], power);
    }
}

struct matrix* transpose_matrix(struct matrix const * const matrix) {
    struct matrix * const transpose = create_matrix(matrix->rows, matrix->cols);
    if(transpose != NULL) {
        for(size_t i = 0; i < matrix->rows; i++) {
            for(size_t j = 0; j < matrix->cols; j++) {
                transpose->elements[i + j * (matrix->rows)] = matrix->elements[j + i * (matrix->cols)];
            }
        }
    }
    return transpose;
}

struct matrix* matrix_multiply(struct matrix const * const A, struct matrix const * const B) {
    if(A->cols != B->rows) {
        printf("Dimensions mismatch: A: (cols: %ld, rows: %ld), B: (cols: %ld, rows: %ld)\n", A->cols, A->rows, B->cols, B->rows);
        return NULL;
    }
    struct matrix * const result = create_matrix(B->cols, A->rows);
    if(result != NULL) {
        for(size_t i = 0; i < A->rows; i++) {
            for(size_t j = 0; j < B->cols; j++) {
                result->elements[j + i * B->cols] = 0.0L;
                for(size_t k = 0; k < A->cols; k++) {
                    result->elements[j + i * B->cols] += A->elements[i * A->cols + k] * B->elements[k * B->cols + j];
                }
            }
        }
    }
    return result;
}

/* Helper for Gaussian reduction */
static char matrix_scale_row(struct matrix * const matrix, struct matrix * const secondary_matrix, size_t const row, size_t const col, long double const scale_to)
{
    if(row == 0L || row > matrix->rows || col == 0L || col > matrix->cols) {
        return 1;
    }
    long double element = matrix->elements[(row - 1L) * matrix->cols + col - 1L];
    if(element == 0.0L) {
        for(size_t i = 0L; i != matrix->rows; i++) {
            element = (matrix->elements[i * matrix->cols + col - 1L]);
            if(element != 0.0L) {
                for(size_t j = 0L; j < matrix->cols; j++) {
                    matrix->elements[(row - 1L)*matrix->cols + i] += matrix->elements[i * matrix->cols + col];
                    secondary_matrix->elements[(row - 1L)*matrix->cols + i] += secondary_matrix->elements[i * matrix->cols + col];
                }
            }
        }
        if(element == 0.0L) {
            return 1;
        }
        element = matrix->elements[(row - 1L) * matrix->cols + col - 1L];
    }
    long double const scale = scale_to / element;
    if(!isfinite(scale)) {
        return 1;
    }
    if(scale != 1.0L) {
        for(size_t i = 0; i < matrix->cols; i++) {
            matrix->elements[(row - 1L)*matrix->cols + i] *= scale;
            secondary_matrix->elements[(row - 1L)*matrix->cols + i] *= scale;
        }
    }
    return 0;
}

/* Helper for Gaussian reduction */
static char matrix_reduce_row(struct matrix * const matrix, struct matrix * const secondary_matrix, size_t const row_dst, size_t const row_src, size_t const col)
{
    if(row_src == row_dst || row_dst == 0L || row_dst > matrix->rows || row_src == 0L || row_src > matrix->rows || col == 0L || col > matrix->cols) {
        return 1;
    }
    if(matrix->elements[(row_src - 1L)*matrix->cols + col - 1L] == 0.0L) {
        return 1;
    }
    long double const scale = matrix->elements[(row_dst - 1L) * matrix->cols + col - 1L];
    if(scale != 0.0L) {
        for(size_t i = 0L; i < matrix->cols; i++) {
            matrix->elements[(row_dst - 1L)*matrix->cols + i] -= scale * matrix->elements[(row_src - 1L) * matrix->cols + i];
            secondary_matrix->elements[(row_dst - 1L)*matrix->cols + i] -= scale * secondary_matrix->elements[(row_src - 1L) * matrix->cols + i];
        }
    }
    return 0;
}

struct matrix* matrix_inverse(struct matrix const * const matrix) {
    if(matrix->rows != matrix->cols) {
        return NULL;
    }
    struct matrix *const temp = create_matrix(matrix->rows, matrix->cols);
    if(temp != NULL) {
        struct matrix *const inverse = create_matrix(matrix->rows, matrix->cols);
        if(inverse != NULL) {
            char result = 0;
            /* Duplicate matrix into temp */
            memcpy(temp->elements, matrix->elements, matrix->rows * matrix->cols * sizeof(long double));
            /* Create identify matrix */
            for(size_t i = 0; i < matrix->rows; i++) {
                for(size_t j = 0; j < matrix->cols; j++) {
                    inverse->elements[j + i * matrix->cols] = (i == j) ? 1.0L : 0.0L;
                }
            }
            /* Do Gaussian reduction */
            for(size_t i = 1; i <= matrix->rows && result == 0; i++) {
                result = matrix_scale_row(temp, inverse, i, i, 1.0L);
                for(size_t j = i + 1; j <= matrix->rows && result == 0; j++) {
                    result = matrix_reduce_row(temp, inverse, j, i, i);
                }
            }
            for(size_t i = matrix->rows; i >= 1 && result == 0; i--) {
                for(size_t j = i - 1; j >= 1 && result == 0; j--) {
                    result = matrix_reduce_row(temp, inverse, j, i, i);
                }
            }
            if(result != 0) {
                destroy_matrix(temp);
                destroy_matrix(inverse);
                return NULL;
            }
            destroy_matrix(temp);
            return inverse;
        } else {
            destroy_matrix(temp);
        }
    }
    return NULL;
}

void print_matrix(struct matrix const * const matrix)
{
    printf("Matrix: %lu × %lu\n", matrix->rows, matrix->cols);
    printf("----BEGIN MATRIX----\n");
    for(size_t i = 0; i < matrix->rows; i++) {
        printf("[ ");
        for(size_t j = 0; j < matrix->cols; j++) {
            printf("%.3Lf ", matrix->elements[i * matrix->cols + j]);
        }
        printf("]\n");
    }
    printf("-----END MATRIX-----\n");
}
