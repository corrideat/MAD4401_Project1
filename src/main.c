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

#include <locale.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "utilities.h"
#include "project1.h"

#define TOLERANCE 1E-7L
#define TOLERANCE_3 1.0L/35184372088832.0L
#define EXPORT_POINTS 524288L

/* The function we are interested in for this project (1-3) */
static long double f(long double x)
{
    return expl(-x / 5.0L) - sinl(x);
}

static long double f_derivative(long double x)
{
    return -expl(-x / 5.0L) / -5.0L - cosl(x);
}

static long double f_secondderivative(long double x)
{
    return -expl(-x / 5.0L) / 25.0L + sinl(x);
}

/* The function we are interested in for this project (4) */
static long double g(long double x)
{
    return powl(x - 3.0L, 4.0L) * sinl(x);
}

static long double g_derivative(long double x)
{
    return powl(x - 3.0L, 3.0L) * (4.0L * sinl(x) + (x - 3.0L) * cosl(x));
}

/* The function we are interested in for this project (bounus) */
static long double f_bonus(long double x)
{
    return powl(x - 4.0L, 2.0L) * sinl(x);
}

static long double f_bonus_derivative(long double x)
{
    return (x - 4.0L) * (2.0L * sinl(x) + (x - 4.0L) * cosl(x));
}

static long double g_bonus(long double x)
{
    return powl(x - 4.0L, 3.0L) * sinl(x);
}

static long double g_bonus_derivative(long double x)
{
    return powl(x - 4.0L, 2.0L) * (3.0L * sinl(x) + (x - 4.0L) * cosl(x));
}

/* The function we are interested in interpolating */
static long double h(long double x)
{
    return 1.0L / (powl(x, 2.0L) + 1.0L);
}

struct function const study_functions[] = {
    {
        "e**(-x/5)/sin(x)",
        (long double(*)(long double, void const*))f,
        NULL
    },
    {
        "(x-3)**4*sin(x)",
        (long double(*)(long double, void const*))g,
        NULL
    }
};

struct function const study_function_derivatives[] = {
    {
        "(e**(-x/5)/sin(x))'",
        (long double(*)(long double, void const*))f_derivative,
        NULL
    },
    {
        "(x-3)**4*sin(x)'",
        (long double(*)(long double, void const*))g_derivative,
        NULL
    },
};

struct function const study_function_secondderivatives[] = {
    {
        NULL,
        NULL,
        NULL
    },
    {
        "(e**(-x/5)/sin(x))''",
        (long double(*)(long double, void const*))f_secondderivative,
        NULL
    },
};

struct function const bonus_functions[] = {
    {
        "(x-4)**2*sin(x)",
        (long double(*)(long double, void const*))f_bonus,
        NULL
    },
    {
        "(x-4)**3*sin(x)",
        (long double(*)(long double, void const*))g_bonus,
        NULL
    }
};

struct function const bonus_function_derivatives[] = {
    {
        "(x-4)**2*sin(x)'",
        (long double(*)(long double, void const*))f_bonus_derivative,
        NULL
    },
    {
        "(x-4)**3*sin(x)",
        (long double(*)(long double, void const*))g_bonus_derivative,
        NULL
    },
};

struct function const interpolation_function = {
    "1/(1+x**2)",
    (long double(*)(long double, void const*))h,
    NULL
};

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "");

    gnuplot("1_visual_inspection", 1, 0.0L, 10.0L, EXPORT_POINTS, &study_functions[0]);

    struct result const bisection_result[] = {
        bisection_method(&study_functions[0], 0.5L, 1.5L, TOLERANCE),
        bisection_method(&study_functions[0], 2.0L, 3.0L, TOLERANCE),
        bisection_method(&study_functions[0], 6.0L, 7.0L, TOLERANCE),
        bisection_method(&study_functions[0], 9.0L, 10.0L, TOLERANCE),
    };

    printf("Bisection Method: %s\n", study_functions[0].name);
    for(int i = 0; i != 4; i++) {
        report_result(&bisection_result[i]);
    }

    struct result const newtons_result[] = {
        newtons_method(&study_functions[0], &study_function_derivatives[0], 1.0L, bisection_result[0].iterations * 4, TOLERANCE),
        newtons_method(&study_functions[0], &study_function_derivatives[0], 2.5L, bisection_result[1].iterations * 4, TOLERANCE),
        newtons_method(&study_functions[0], &study_function_derivatives[0], 6.5L, bisection_result[2].iterations * 4, TOLERANCE),
        newtons_method(&study_functions[0], &study_function_derivatives[0], 9.9L, bisection_result[3].iterations * 4, TOLERANCE)
    };

    printf("Newton's Method: %s\n", study_functions[0].name);
    for(int i = 0; i != 4; i++) {
        report_result(&newtons_result[i]);
    }

    struct result const newtons_result_3 = newtons_method(&study_functions[1], &study_function_derivatives[1], 2.0L, 256, TOLERANCE_3);
    printf("Newton's Method (part 3): %s\n", study_functions[1].name);
    report_result(&newtons_result_3);

    struct result const altered_newtons_result_3 = altered_newtons_method(&study_functions[1], &study_function_derivatives[1], &study_function_secondderivatives[1], 2.0L, 256, TOLERANCE_3);
    printf("Altered Newton's Method (part 3): %s\n", study_functions[1].name);
    report_result(&altered_newtons_result_3);

    struct interpolation const* lagrange[] = {
        lagrange_interpolation(&interpolation_function, -5.0L, 5.0L, 5),
        lagrange_interpolation(&interpolation_function, -5.0L, 5.0L, 10),
        lagrange_interpolation(&interpolation_function, -5.0L, 5.0L, 20)
    };
    struct function const lagrange_functions[] = {
        {
            lagrange[0]->name,
            (long double(*)(long double, void const*))polynomial_value,
            lagrange[0]
        },
        {
            lagrange[1]->name,
            (long double(*)(long double, void const*))polynomial_value,
            lagrange[1]
        },
        {
            lagrange[2]->name,
            (long double(*)(long double, void const*))polynomial_value,
            lagrange[2]
        }
    };
    gnuplot("lagrange", 4, -5.0L, 5.0L, EXPORT_POINTS, &interpolation_function, &lagrange_functions[0], &lagrange_functions[1], &lagrange_functions[2]);

    printf("Lagrange interpolation coefficients for %s\n", interpolation_function.name);
    for(int i = 0; i != 3; i++) {
        for(ssize_t j = lagrange[i]->order; j >= 0; j--) {
            printf("%.4LE x**%ld%s", fabsl(lagrange[i]->coefficients[j]), j, j == 0 ? "" : (lagrange[i]->coefficients[j - 1] < 0.0L ? " - " : " + "));
        }
        printf(" order: %ld, error: %.2LE\n", lagrange[i]->order, polynomial_error(lagrange[i]));
        destroy_interpolation((struct interpolation*)lagrange[i]);
    }

    struct interpolation const* piecewise_linear[] = {
        piecewise_linear_interpolation(&interpolation_function, -5.0L, 5.0L, 5),
        piecewise_linear_interpolation(&interpolation_function, -5.0L, 5.0L, 10),
        piecewise_linear_interpolation(&interpolation_function, -5.0L, 5.0L, 20)
    };
    struct function const piecewise_linear_functions[] = {
        {
            piecewise_linear[0]->name,
            (long double(*)(long double, void const*))piecewise_linear_value,
            piecewise_linear[0]
        },
        {
            piecewise_linear[1]->name,
            (long double(*)(long double, void const*))piecewise_linear_value,
            piecewise_linear[1]
        },
        {
            piecewise_linear[2]->name,
            (long double(*)(long double, void const*))piecewise_linear_value,
            piecewise_linear[2]
        }
    };
    gnuplot("piecewise_linear", 4, -5.0L, 5.0L, EXPORT_POINTS, &interpolation_function, &piecewise_linear_functions[0], &piecewise_linear_functions[1], &piecewise_linear_functions[2]);

    printf("Piecewise linear interpolation coefficients for %s\n", interpolation_function.name);
    for(int i = 0; i != 3; i++) {
        printf(" order: %ld, error: %.2LE\n", piecewise_linear[i]->order, piecewise_linear_error(piecewise_linear[i]));
        destroy_interpolation((struct interpolation*)piecewise_linear[i]);
    }

    struct interpolation const* raised_cosine[] = {
        raised_cosine_interpolation(&interpolation_function, -5.0L, 5.0L, 5),
        raised_cosine_interpolation(&interpolation_function, -5.0L, 5.0L, 10),
        raised_cosine_interpolation(&interpolation_function, -5.0L, 5.0L, 20)
    };
    struct function const raised_cosine_functions[] = {
        {
            raised_cosine[0]->name,
            (long double(*)(long double, void const*))raised_cosine_value,
            raised_cosine[0]
        },
        {
            raised_cosine[1]->name,
            (long double(*)(long double, void const*))raised_cosine_value,
            raised_cosine[1]
        },
        {
            raised_cosine[2]->name,
            (long double(*)(long double, void const*))raised_cosine_value,
            raised_cosine[2]
        }
    };
    gnuplot("raised_cosine", 4, -5.0L, 5.0L, EXPORT_POINTS, &interpolation_function, &raised_cosine_functions[0], &raised_cosine_functions[1], &raised_cosine_functions[2]);

    printf("Raised cosine interpolation coefficients for %s\n", interpolation_function.name);
    for(int i = 0; i != 3; i++) {
        printf(" order: %ld, error: %.2LE\n", raised_cosine[i]->order, raised_cosine_error(raised_cosine[i]));
        destroy_interpolation((struct interpolation*)raised_cosine[i]);
    }

    struct interpolation const* least_squares[] = {
        least_squares_interpolation(&interpolation_function, -5.0L, 5.0L, 5),
        least_squares_interpolation(&interpolation_function, -5.0L, 5.0L, 10),
        least_squares_interpolation(&interpolation_function, -5.0L, 5.0L, 20)
    };

    struct function const least_squares_functions[] = {
        {
            least_squares[0]->name,
            (long double(*)(long double, void const*))polynomial_value,
            least_squares[0]
        },
        {
            least_squares[1]->name,
            (long double(*)(long double, void const*))polynomial_value,
            least_squares[1]
        },
        {
            least_squares[2]->name,
            (long double(*)(long double, void const*))polynomial_value,
            least_squares[2]
        }
    };
    gnuplot("least_squares", 4, -5.0L, 5.0L, EXPORT_POINTS, &interpolation_function, &least_squares_functions[0], &least_squares_functions[1], &least_squares_functions[2]);
    printf("Least squares interpolation coefficients for %s\n", interpolation_function.name);
    for(int i = 0; i != 3; i++) {
        for(ssize_t j = least_squares[i]->order; j >= 0; j--) {
            printf("%.4LE x**%ld%s", fabsl(least_squares[i]->coefficients[j]), j, j == 0 ? "" : (least_squares[i]->coefficients[j - 1] < 0.0L ? " - " : " + "));
        }
        printf(" order: %ld, error: %.2LE\n", least_squares[i]->order, polynomial_error(least_squares[i]));
        destroy_interpolation((struct interpolation*)least_squares[i]);
    }

    /* Bonus Problem 1 */

    char square_root_errors = 0;
    for(long double i = 10.0L; i <= 10000.0L; i += 1.0L / 8192.0L) {
        struct result squareroot = square_root_calculator(i);
        double long real_sqrt = sqrtl(i);
        if(labs(squareroot.value - real_sqrt) > squareroot.error) {
            square_root_errors++;
            printf("Error for square root of %.0Lf (%Lf ± %LE not %Lf) \n", i, squareroot.value, squareroot.error, real_sqrt);
        }
    }
    if(square_root_errors == 0) {
        printf("Success for square root\n");
    }

    /* Bonus Problem 2 */
    struct result const bonus_newtons_result[] = {
        newtons_method(&bonus_functions[0], &bonus_function_derivatives[0], 5.0L, 256, TOLERANCE),
        newtons_method(&bonus_functions[1], &bonus_function_derivatives[1], 5.0L, 256, TOLERANCE)
    };
    struct result const adjusting_bonus_newtons_result[] = {
        adjusting_newtons_method(&bonus_functions[0], &bonus_function_derivatives[0], 5.0L, 256, TOLERANCE),
        adjusting_newtons_method(&bonus_functions[1], &bonus_function_derivatives[1], 5.0L, 256, TOLERANCE)
    };

    printf("Bonus Problem 2: Adjusting Newton's Method\n");
    for(int i = 0; i != 2; i++) {
        printf("function %s:\n", bonus_functions[i].name);
        report_result(&bonus_newtons_result[i]);
        report_result(&adjusting_bonus_newtons_result[i]);
    }

    return EXIT_SUCCESS;
}
