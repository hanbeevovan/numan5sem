import numpy as np
from math import sin, cos, exp, pi, atan

n = 20 
eps = 0.1
a = 0
b = 1

def f(x):
    return -x * (1 - exp(-cos(x)) + 2 * x * cos(x))


def funcA(left, right, n):
    return [left + i * (right - left) / n for i in range(n + 1)]

def funcB(left, right, n):
    return [0.5 * (right + left - (right - left) * cos((2 * k + 1) / (2 * (n + 1)) * pi)) for k in range(n + 1)]

def spline(x, x_nodes, f_nodes, m_coefs, n_local):
    if x <= x_nodes[0]:
        i = 1
    elif x >= x_nodes[n_local]:
        i = n_local
    else: 
        i= 1 
        while x > x_nodes[i]:
            i += 1
    
    h_i = x_nodes[i] - x_nodes[i - 1]
    if h_i == 0:
        return f_nodes[i]
    
    x_i = x_nodes[i]
    x_i_prev = x_nodes[i - 1]
    f_i = f_nodes[i]
    f_i_prev = f_nodes[i - 1]
    m_i = m_coefs[i]
    m_i_prev = m_coefs[i - 1]

    term1 = m_i_prev * ((x_i - x) ** 3) / (6 * h_i)
    term2 = m_i * ((x - x_i_prev) ** 3) / (6 * h_i)
    term3 = (f_i_prev - (h_i ** 2 / 6) * m_i_prev) * (x_i - x) / h_i
    term4 = (f_i - (h_i ** 2 / 6) * m_i) * (x - x_i_prev) / h_i

    S_x = term1 + term2 + term3 + term4

    return S_x

def progonka(l,m,diag,d):
    n_local = len(m) - 1 

    c_prime = list(diag) + [0.0]
    d_prime = list(d)
    m_coefs = [0.0] * (n_local + 1)

    c_prime[0] = c_prime[0] / m[0]
    d_prime[0] = d_prime[0] / m[0]

    for i in range(1, n_local + 1):
        temp = m[i] - l[i - 1] * c_prime[i - 1]

        if i < n_local:
            c_prime[i] = c_prime[i] / temp

        d_prime[i] = (d_prime[i] - l[i - 1] * d_prime[i - 1]) / temp
    
    m_coefs[n_local] = d_prime[n_local]

    for i in range(n_local - 1, -1, -1):
        m_coefs[i] = d_prime[i] - c_prime[i] * m_coefs[i + 1]
    
    return m_coefs


def show(S, INTERPOL, FUNC, e=0.0):
    print("=" * 70)
    print(f"ПОБУДОВА СПЛАЙНА: {S}")
    print(f"Відрізок [a, b] = [{a}, {b}], Кількість вузлів n = {n}, Зміщення e = {e}")
    print("-" * 70)
    
    x_nodes = FUNC(a, b, n)
    
    if e != 0.0:
        x_nodes = [x_nodes[0]] + [x + e for x in x_nodes[1:-1]] + [x_nodes[-1]]
        
    f_nodes = [f(x) for x in x_nodes]

    h = [0.0] * (n + 1)
    for i in range(1, n + 1):
        h[i] = x_nodes[i] - x_nodes[i - 1]

    l = [0.0] * n
    m = [0.0] * (n + 1)
    diag = [0.0] * n
    d = [0.0] * (n + 1)

    m[0] = 1.0
    diag[0] = 0.0
    d[0] = 0.0

    for i in range(1, n):
        l[i-1] = h[i] / 6.0
        m[i] = (h[i] + h[i+1]) / 3.0
        diag[i] = h[i+1] / 6.0
        d[i] = (f_nodes[i+1] - f_nodes[i]) / h[i+1] - (f_nodes[i] - f_nodes[i-1]) / h[i]

    l[n-1] = 0.0
    m[n] = 1.0
    d[n] = 0.0

    m_coeffs = progonka(l, m, diag, d)
    

    
    alpha = 0.5
    h_test = (b - a) / n
    x_nodes_uniform = funcA(a, b, n)
    test_points = [x_nodes_uniform[i] + alpha * h_test for i in range(n)]

    print(f"{'x_tilde':>15} | {'f(x_tilde)':>15} | {'S(x_tilde)':>15} | {'|f - S|':>15}")
    print("-" * 70)
    
    max_error = 0.0
    for xt in test_points:
        f_xt = f(xt)
        

        S_xt = INTERPOL(xt, x_nodes, f_nodes, m_coeffs, n)
        
        error = abs(f_xt - S_xt)
        
        if error > max_error:
            max_error = error
            
        print(f"{xt:15.8f} | {f_xt:15.8f} | {S_xt:15.8f} | {error:15.8e}")
        
    print("-" * 70)
    print(f"Максимальна похибка на тестових точках: {max_error:15.8e}")
    print("=" * 70 + "\n")

show("S(X):A", spline, funcA)
show("S(X):B", spline, funcB)
show("S(x):A+e", spline, funcA, eps)
show("S(x):B+e", spline, funcB, eps)