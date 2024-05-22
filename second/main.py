import math

STEPS = 1024
LEFT_BORDER = 0
EPSILON = 1e-6
E = 2.71828182846

def Tfunc(x):
    n = 0
    node_0 = x
    ans = x
    while abs(node_0) > 1e-6:
        q = ((-1) * x * x * (2 * n + 1)) / ((2 * n + 2) * (2 * n + 3) * (2 * n + 3))
        node_0 *= q
        ans += node_0
        n += 1
    return ans

def func(t):
    if t != 0:
        return math.sin(t) / t
    return 1

def leftRectangles(func, a, b, steps):
    result = 0.0
    stepSize = (b - a) / steps
    for i in range(steps):
        x_i = a + stepSize * i
        result += stepSize * func(x_i)
    return result

def rightRectangles(func, a, b, steps):
    result = 0.0
    stepSize = (b - a) / steps
    for i in range(1, steps + 1):
        x_i_1 = a + stepSize * i
        result += stepSize * func(x_i_1)
    return result

def middleRectangles(func, a, b, steps):
    result = 0.0
    stepSize = (b - a) / steps
    for i in range(steps):
        x_i = a + stepSize * i
        x_i_1 = a + stepSize * (i + 1)
        result += stepSize * func((x_i + x_i_1) / 2)
    return result

def trapezeFormula(func, a, b, steps):
    result = func(a) + func(b)
    stepSize = (b - a) / steps
    for i in range(1, steps):
        x_i_1 = a + stepSize * i
        result += 2 * func(x_i_1)
    result *= stepSize / 2
    return result

def SypmsonsFormula(func, a, b, steps):
    stepSize = (b - a) / steps
    result = 0
    x = a
    for i in range(steps):
        result += (func(x) + 4 * func(x + stepSize / 2) + func(x + stepSize)) * stepSize / 6
        x += stepSize
    return result

def GaussFormula(func, a, b, steps):
    stepSize = (b - a) / steps
    ad1 = (1 - 1.0 / math.sqrt(3)) * stepSize / 2
    ad2 = (1 + 1.0 / math.sqrt(3)) * stepSize / 2
    result = 0
    x = a
    for i in range(steps):
        result += (func(x + ad1) + func(x + ad2)) * stepSize / 2
        x += stepSize
    return result

def calculateFunc(points, function):
    for point in points:
        steps = 1
        lastResult = 0.0
        currentResult = 0.0
        while True:
            steps *= 2
            lastResult = currentResult
            currentResult = function(func, LEFT_BORDER, point, steps)
            if abs(lastResult - currentResult) <= EPSILON or steps >= STEPS:
                break

        difference = abs(Tfunc(point) - currentResult)

        print(f"x_i = {point:.1f} | J_o = {Tfunc(point):.10f} | J_n = {currentResult:.10f} | |J_o - J_n| = {difference:.10f} | N = {steps}")

if __name__ == "__main__":
    points = [0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0]
    print("Правые прямоугольники")
    calculateFunc(points, rightRectangles)
    print("Левые прямоугольники")
    calculateFunc(points, leftRectangles)
    print("Центральные прямоугольники")
    calculateFunc(points, middleRectangles)
    print("Трапеции")
    calculateFunc(points, trapezeFormula)
    print("Симпсон")
    calculateFunc(points, SypmsonsFormula)
    print("Гаус")
    calculateFunc(points, GaussFormula)
