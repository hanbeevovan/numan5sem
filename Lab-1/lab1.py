from math import cos, sin, exp, atan, pi

n = 20 
eps = 0.05  
a = -2
b = 2

def f(x):
    return -x * (1 - exp(cos(x) + 2 * x * cos(x)))

# def f(x):
#     return exp(x / 10) * sin(x) + x ** 3 + cos(x)

def funcA(left, right, n):

    return [left + i * (right - left) / n for i in range(n + 1)]

def funcB(left, right, n):

    return [0.5 * (right + left - (right - left) * cos((2 * k + 1) / (2 * (n + 1)) * pi)) for k in range(n + 1)]


def Polynom_Lagrange(x,x_values, func):

    y_values = [func(x) for x in x_values]
    Ln = 0.0

    for i in range(n + 1):
        f_i = y_values[i]
        l_i = 1.0
        
        for k in range(n + 1):
            if i == k:
                continue
            l_i = l_i * (x - x_values[k]) / (x_values[i] - x_values[k])

        Ln = Ln + f_i * l_i

    return Ln

def divided_difference_forward(x_values, y_values, k):
    result = 0
    
    for j in range(k + 1):
        mul = 1
        for i in range(k + 1):
            if i != j:
                mul *= x_values[j] - x_values[i]
        result += y_values[j] / mul
    return result

def divided_difference_backward(x_values, y_values, k):
    result = 0

    for j in range(k + 1):
        mul = 1
        for i in range(k + 1):
            if i != j:
                mul *= x_values[n - j] - x_values[n - i]
        result += y_values[n - j] / mul
    return result

def Polynom_Newton_Forward(x, x_values, func):
    differences = []
    y_values = [func(x) for x in x_values]

    for i in range(1, n + 1):
        differences.append(divided_difference_forward(x_values, y_values, i))
    result = y_values[0]
    for k in range(1, n + 1):
        mul = 1.0
        for j in range(k):
            mul *= (x - x_values[j])
        result += differences[k - 1] * mul
    return result

def Polynom_Newton_Backward(x, x_values, func):
    differences = []
    y_values = [func(x) for x in x_values]

    for i in range(1, n + 1):
        differences.append(divided_difference_backward(x_values, y_values, i))
    result = y_values[n]
    for k in range(1, n + 1):
        mul = 1.0
        for j in range(k):
            mul *= (x - x_values[n - j])
        result += differences[k - 1] * mul
    return result


def show(INTERPOL, FUNC, eps):

    x_values = FUNC(a, b, n)
    h = (b - a) / n
    
    Lx = [INTERPOL(x, x_values, f) for x in x_values]
    fx = [f(x) for x in x_values]
    diff = [abs(Lx[i] - fx[i]) for i in range(n + 1)]
    
    x_values_moved = [x + eps * h for x in x_values[:-1]]
    Lx_moved = [INTERPOL(x, x_values, f) for x in x_values_moved]
    fx_moved = [f(x) for x in x_values_moved]
    diff_moved = [abs(Lx_moved[i] - fx_moved[i]) for i in range(n)]
    
    print("At the nodes")
    print(f"{'x':<15}{'f(x)':<15}{'L(x)':<15}{'|f(x)-L(x)|':<30}")
    for i in range(n + 1):
        print(f"{x_values[i]:<15.10f}{fx[i]:<15.10f}{Lx[i]:<15.10f}{diff[i]:<30.26f}")
    
    print("With offset")
    print(f"{'x':<15}{'f(x)':<15}{'L(x)':<15}{'|f(x)-L(x)|':<30}")
    for i in range(n):
        print(f"{x_values_moved[i]:<15.10f}{fx_moved[i]:<15.10f}{Lx_moved[i]:<15.10f}{diff_moved[i]:<30.26f}")

print("Lagrange with uniform nodes")
show(Polynom_Lagrange, funcA, eps)
print("\nLagrange with Chebyshev nodes")
show(Polynom_Lagrange, funcB, eps)  
print("\nNewton Forward with uniform nodes")
show(Polynom_Newton_Forward, funcA, eps)
print("\nNewton Forward with Chebyshev nodes")
show(Polynom_Newton_Forward, funcB, eps)        
print("\nNewton Backward with uniform nodes")
show(Polynom_Newton_Backward, funcA, eps)   
print("\nNewton Backward with Chebyshev nodes")
show(Polynom_Newton_Backward, funcB, eps)
