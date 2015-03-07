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

struct matrix {
    long double *elements;
    size_t cols;
    size_t rows;
};

struct matrix* create_matrix(size_t cols, size_t rows);
void destroy_matrix(struct matrix * const);
void matrix_set_row(struct matrix * matrix, size_t row, long double);
void matrix_set_row_vector(struct matrix *matrix, size_t row, long double const* vector);
void matrix_set_row_vector_power(struct matrix * matrix, size_t row, long double const* vector, long double power);
struct matrix* transpose_matrix(struct matrix const* matrix);
struct matrix* matrix_multiply(struct matrix const* A, struct matrix const* B);
struct matrix* matrix_inverse(struct matrix const *);
void print_matrix(struct matrix const * const matrix);
