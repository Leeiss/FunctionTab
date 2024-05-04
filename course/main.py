import numpy as np
import matplotlib.pyplot as plt

# Функция для решения системы ОДУ методом Рунге-Кутты 4-го порядка точности
def runge_kutta(f, t0, y0, h, tf):
    """
    Функция для решения системы ОДУ методом Рунге-Кутты 4-го порядка точности.
    
    Параметры:
        - f: функция, правая часть системы ОДУ, возвращает np.array с производными
        - t0: начальное значение независимой переменной
        - y0: начальное значение зависимых переменных в виде np.array
        - h: шаг интегрирования
        - tf: конечное значение независимой переменной
        
    Возвращает:
        - массив временных значений
        - массив значений зависимых переменных в каждый момент времени
    """
    t = np.arange(t0, tf + h, h)
    y = np.zeros((len(t), len(y0)))
    y[0] = y0

    for i in range(1, len(t)):
        k1 = f(t[i-1], y[i-1])
        k2 = f(t[i-1] + h/3, y[i-1] + (h * k1) / 3)
        k3 = f(t[i-1] + 2*h/3, y[i-1] - (h * k1) / 3 + h * k2)
        k4 = f(t[i-1] + h, y[i-1] + h * k1 - h * k2 + h * k3)
        y[i] = y[i-1] + h * (k1 + 3*k2 + 3*k3 + k4) / 8

    return t, y

# Определение правой части системы ОДУ для тестовой задачи
def f_test(t, y):
    y1, y2 = y
    dy1dt = y1 / (2 + 2*t) - 2*t*y2
    dy2dt = y2 / (2 + 2*t) + 2*t*y1
    return np.array([dy1dt, dy2dt])

# Начальные условия и параметры для тестовой задачи
t0_test = 0
tf_test = 2
y0_test = np.array([np.cos(t0_test**2) * np.sqrt(1 + t0_test), np.sin(t0_test**2) * np.sqrt(1 + t0_test)])
h_test = 0.01

# Решение тестовой задачи
t_test, y_test = runge_kutta(f_test, t0_test, y0_test, h_test, tf_test)

# Точные решения для сравнения
exact_solution = np.array([[np.cos(t**2) * np.sqrt(1 + t) for t in t_test], [np.sin(t**2) * np.sqrt(1 + t) for t in t_test]])

# Построение графика численного решения
plt.figure(figsize=(10, 5))
plt.plot(t_test, y_test[:,0], label='Численное решение для $y_1$')
plt.plot(t_test, y_test[:,1], label='Численное решение для $y_2$')
plt.xlabel('t')
plt.ylabel('y')
plt.title('Численное решение тестовой задачи методом Рунге-Кутты')
plt.legend()
plt.grid(True)
plt.savefig('plot1.png')

# Построение графика точного решения
plt.figure(figsize=(10, 5))
plt.plot(t_test, exact_solution[0], label='Точное решение для $y_1$', linestyle='--')
plt.plot(t_test, exact_solution[1], label='Точное решение для $y_2$', linestyle='--')
plt.xlabel('t')
plt.ylabel('y')
plt.title('Точное решение тестовой задачи')
plt.legend()
plt.grid(True)
plt.savefig('plot2.png')