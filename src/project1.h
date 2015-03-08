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


struct result bisection_method(struct function const* function, long double x0, long double x1, long double tolerance);

struct result newtons_method(struct function const* function, struct function const* derivative, long double x0, unsigned long max_iterations, long double tolerance);
struct result altered_newtons_method(struct function const* function, struct function const* derivative, struct function const* secondderivative, long double x0, unsigned long max_iterations, long double tolerance);
struct interpolation const* lagrange_interpolation(struct function const* function, long double x0, long double x1, unsigned long order);
struct interpolation const* piecewise_linear_interpolation(struct function const* function, long double x0, long double x1, unsigned long order);
struct interpolation const* raised_cosine_interpolation(struct function const* function, long double x0, long double x1, unsigned long order);
struct interpolation const* least_squares_interpolation(struct function const* function, long double x0, long double x1, unsigned long order);
struct result square_root_calculator(double long const k);