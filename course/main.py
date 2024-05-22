import numpy as np
import matplotlib.pyplot as plt

def fx_test(t, x, y):
    return x / (2 + 2 * t) - 2 * t * y #x тестовое

def fy_test(t, x, y):
    return y / (2 + 2 * t) + 2 * t * x #у тестовое

def fx_dynamic(t, x, y, alfa): # x в безмерной форме
    eps = 0.1
    return (1 - eps * x) * x - (x * y) / (1 + alfa * x) 

def fy_dynamic(t, x, y, alfa): # y в безмерной форме
    gamma = 1
    return gamma * (x / (1 + alfa * x) - 1) * y

def solve_test(t1, t2, h, x0, y0):
    n = int((t2 - t1) / h + 1)#вычислем количество шагов
    t = np.linspace(t1, t2, n)#создаем массив временных значений t
    x = np.zeros(n)#создаем массивы с нулями
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0

    for i in range(n - 1):
        xk1 = fx_test(t[i], x[i], y[i])
        yk1 = fy_test(t[i], x[i], y[i])
        xk2 = fx_test(t[i] + h / 3, x[i] + h * xk1 / 3, y[i] + h * yk1 / 3)
        yk2 = fy_test(t[i] + h / 3, x[i] + h * xk1 / 3, y[i] + h * yk1 / 3)
        xk3 = fx_test(t[i] + 2 * h / 3, x[i] - h * xk1 / 3 + h * xk2, y[i] - h * yk1 / 3 + h * yk2)
        yk3 = fy_test(t[i] + 2 * h / 3, x[i] - h * xk1 / 3 + h * xk2, y[i] - h * yk1 / 3 + h * yk2)
        xk4 = fx_test(t[i] + h, x[i] + h * (xk1 - xk2 + xk3), y[i] + h * (yk1 - yk2 + yk3))
        yk4 = fy_test(t[i] + h, x[i] + h * (xk1 - xk2 + xk3), y[i] + h * (yk1 - yk2 + yk3))
        x[i + 1] = x[i] + h * (xk1 + 3 * xk2 + 3 * xk3 + xk4) / 8
        y[i + 1] = y[i] + h * (yk1 + 3 * yk2 + 3 * yk3 + yk4) / 8

    return t, x, y

def solve_dynamic(t1, t2, h, x0, y0, alfa):
    n = int((t2 - t1) / h + 1)
    t = np.linspace(t1, t2, n)
    x = np.zeros(n)
    y = np.zeros(n)
    x[0] = x0
    y[0] = y0

    for i in range(n - 1):
        xk1 = fx_dynamic(t[i], x[i], y[i], alfa)
        yk1 = fy_dynamic(t[i], x[i], y[i], alfa)
        xk2 = fx_dynamic(t[i] + h / 3, x[i] + h * xk1 / 3, y[i] + h * yk1 / 3, alfa)
        yk2 = fy_dynamic(t[i] + h / 3, x[i] + h * xk1 / 3, y[i] + h * yk1 / 3, alfa)
        xk3 = fx_dynamic(t[i] + 2 * h / 3, x[i] - h * xk1 / 3 + h * xk2, y[i] - h * yk1 / 3 + h * yk2, alfa)
        yk3 = fy_dynamic(t[i] + 2 * h / 3, x[i] - h * xk1 / 3 + h * xk2, y[i] - h * yk1 / 3 + h * yk2, alfa)
        xk4 = fx_dynamic(t[i] + h, x[i] + h * (xk1 - xk2 + xk3), y[i] + h * (yk1 - yk2 + yk3), alfa)
        yk4 = fy_dynamic(t[i] + h, x[i] + h * (xk1 - xk2 + xk3), y[i] + h * (yk1 - yk2 + yk3), alfa)
        x[i + 1] = x[i] + h * (xk1 + 3 * xk2 + 3 * xk3 + xk4) / 8#обновляем значения на следующем шаге
        y[i + 1] = y[i] + h * (yk1 + 3 * yk2 + 3 * yk3 + yk4) / 8

    return t, x, y

def exact(t1, t2, h):#вычисляем точное решение
    n = int((t2 - t1) / h + 1)
    t = np.linspace(t1, t2, n)
    x = np.zeros(n)
    y = np.zeros(n)
    for i in range(n):
        x[i] = np.cos(t[i]**2) * np.sqrt(1 + t[i])
        y[i] = np.sin(t[i]**2) * np.sqrt(1 + t[i])
    return x, y

def main_test():
    t1 = 0
    t2 = 2
    x0 = 1 #присутствует определенное количество жертв
    y0 = 0 #хищников пока нет
    n = 25

    
    N = 100
    h = (t2 - t1) / N
    t, x, y = solve_test(t1, t2, h, x0, y0)

    plt.figure()
    plt.plot(t, x, label='x(t)')
    plt.plot(t, y, label='y(t)')
    plt.legend()
    plt.savefig('1.png')
    plt.show()


    x1, y1 = exact(t1, t2, h)
    plt.figure()
    plt.plot(t, x1, label='x(t)')
    plt.plot(t, y1, label='y(t)')
    plt.legend()
    plt.savefig('2.png')
    plt.show()

    
    plt.figure()
    plt.plot(t, x, label='x(t) - RK4')
    plt.plot(t, y, label='y(t) - RK4')
    plt.plot(t, x1, linestyle='--', label='xExact(t)') 
    plt.plot(t, y1, linestyle='--', label='yExact(t)')  
    plt.legend()
    plt.savefig('12overlay.png')
    plt.show()


    N = np.zeros(n)
    e = np.zeros(n)
    h = np.zeros(n)
    e4 = np.zeros(n)

    for i in range(25):#генерация шага n для лучшего отслеживания
        N[i] = i + 1
    for i in range(25, n):
        N[i] = N[i - 1] + 25
        #Вычисление шага h для каждого значения 𝑁
    for i in range(n):
        h[i] = (t2 - t1) / N[i]

    for i in range(n):
        t, x, y = solve_test(t1, t2, h[i], x0, y0)
        x1, y1 = exact(t1, t2, h[i])
        x2 = np.abs(x - x1)
        y2 = np.abs(y - y1)
        maxx = np.max(x2)
        maxy = np.max(y2)
        e[i] = max(maxx, maxy)

    for i in range(len(e)):
        h4 = h[i]**4
        e4[i] = e[i] / h4


    plt.figure()
    plt.plot(N, e, label='e')
    plt.legend()
    plt.savefig('3.png')
    plt.show()

    plt.figure()
    plt.plot(N, e4, label='e/h^4')
    plt.legend()
    plt.savefig('4.png')
    plt.show()

def main_dynamic():
    t1 = 0
    t2 = 100
    x0 = 3
    y0 = 1
    alfa = 0.1

    for i in range(9):
        t, x, y = solve_dynamic(t1, t2, 0.1, x0, y0, alfa)

        plt.figure()
        plt.plot(x, y, label='(X,Y)')
        plt.legend()
        plt.savefig(f'{i + 5}.png')
        plt.show()
        
        alfa += 0.1

    alfa = 0.1

    for i in range(9):
        t, x, y = solve_dynamic(t1, t2, 0.1, x0, y0, alfa)

        plt.figure()
        plt.plot(t, x, label='X(t)')
        plt.plot(t, y, label='Y(t)')
        plt.legend()
        plt.savefig(f'n{i + 5}.png')
        plt.show()
        
        alfa += 0.1

if __name__ == "__main__":
    main_test()
    main_dynamic()
