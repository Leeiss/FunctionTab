## Табулирование трансцендентных функций.

Одна из специальных функций математической физики - интегральный синус, определяется следующим образом

$$
\text{Si}(x) = \int_{0}^{x} \frac{\sin(t)}{t} \, dt
$$


Цель задания - изучить и сравнить различные способы приближенного вычисления этой функции.       
Для этого:

1. Протабулировать Si(x) на отрезке [a,b] с шагом и с точностью є, основываясь на
ряде Тейлора, предварительно вычислив его

$$
\text{Si}(x) = \sum_{n=0}^{\infty} \frac{(-1)^n x^{2n+1}}{(2n+1)(2n+1)!}    
$$


где а = 0, b = 4, h = 0.3, e = 10^(-6), и получить, таким образом, таблицу     

|   x₀   |   x₁   |   x₂   |  ...  |   xₙ   |
|-------|-------|-------|-------|-------|
|   f₀   |   f₁   |   f₂   |  ...  |   fₙ   |




fᵢ=Si(x), xᵢ=a+i*h, i=0,…,n.

______

2. По полученной таблице значений построить интерполяционный полином Ньютона, приближающий Si(x)
        
$$
L_n(x) = f(x_0) + (x - x_0) f(x_0, x_1) + (x - x_0)(x - x_1) f(x_0, x_1, x_2) + \ldots + (x - x_0)(x - x_1)\ldots(x - x_n) f(x_0, x_1, \ldots, x_n)
$$


и вычислить погрешность интерполирования   
   
$$
e_n = \max_{x \in (a, b)} e(x), \quad e(x) = \left| \text{Si}(x) - L_n(x) \right|
$$  


   ______
3. На той же сетке узлов
          
$$
\{x_i\}_{i=0}^{n}
$$

 построить таблицу приближенных значений J₀(x), используя составную квадратурную формулу трапеций

$$
\int_c^d \phi(t) \, dt \approx \sum_{i=1}^{N} \int_{z_{i-1}}^{z_i} \phi(t) \, dt \approx \sum_{i=1}^{N} S_i(\phi)
$$    
  
где   

  $$
S_i(\phi) = h_N\left(\frac{\phi(z_{i-1}) + \phi(z_i)}{2}\right)
$$

zᵢ - точки разбиения отрезка интегрирования на N частей, zᵢ = с + i*h_N, h_N = (c-D)/N
 
Интеграл вычислить с точностью є = 10^(-4). Точность вычисления интеграла определя- ется сравнением результатов при различном числе разбиения отрезка интегрирования.
Именно, точность в считается достигнутой, если

$$
|S_N(\phi) - S^{2N}(\phi)| \leq e
$$

$$
S^N(\phi) = \sum_{i=1}^{N} S_i(\phi)
$$

______

4. Построить таблицу обратной к Si(x) функции F=Si(x) функции F = Si^(-1)
   
|   F₀   |   F₁   |   F₂   |  ...  |   Fₙ   |
|-------|-------|-------|-------|-------|
|   z₀   |   z₁   |   z₂   |  ...  |   zₙ   |

решая уравнения   
Si(z) = F₁, i = 0,...,n, F₁ = f₀ +i*(f_n - f_0)/n   
следующим итерационным методом: решение уравнения g(z) = 0 находится по формуле  
  
$$
z \overset{k+1}{=} z \overset{k}{-} \frac{1}{p^2} \frac{g(z \overset{k}{+} \pi u \overset{k}{})}{g' \overset{k}{(}y\overset{k}{})}
$$

где

$$
u \overset{k}{=} \frac{g(z \overset{k}{})}{g' \overset{k}{(}z \overset{k}{})}
$$


p = &frac12;&radic;5 и k=0,1,...



В качестве начального приближения z к точке z_i, взять х_i

______

Отчет должен содержать:   
+ постановку задачи и исходные данные;
+ Описание методов решения;
+ графики функций Si(x), Si (x), (x);
+ листинги с программой и результатами счета.
   
**Замечание 1**. При вычислении ряда Тейлора   

$$
\sum_{n=0}^{\infty} a_n
$$


учесть, что каждый последующий
член ряда а_(n+1) получается из предыдущего члена а_n, умножением на некоторую величину q_n, т.е.
an+1 = а_n*q_n
Это позволит избежать переполнения при вычислении факториалов, встречающихся в рассматриваемом ряде.   
    
**Замечание 2**. Разделенные разности f(x_k,..., x_(k-1)) в формуле Ньютона удобно вычислять, используя следующий алгоритм.   
Для к= 1,...,n вычисляем f₁ = (fᵢ - fᵢ₋₁) / (xᵢ - xᵢ₋ₖ), где i = n, ..., k.    
В результате в ячейках f₀, f₁,..., fₙ, будут, очевидно, находиться соответственно       
f(x₀), f(x₀, x₁), ..., f(x₀, x₁, ..., xₙ).

