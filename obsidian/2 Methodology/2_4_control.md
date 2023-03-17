Поступательное управление аппарата является основным способом решения задачи трансфера для сложной структуры окружения. Приближение двигательного управления разделяется на импульсное и непрерывное. Первое приближение справедливо для времени работы двигателя $dt \ll \frac{F_{тяги}}{m \partial \upsilon}$, где $\partial \upsilon$ $-$ характерное изменение скорости, второе $-$ для $dt \sim \frac{F_{тяги}}{m \partial \upsilon}$.

В случае импульсного управления, задача трансфера определена если задано правило использования импульса. В данной работе управление основано на следующем правиле: использование одного импульса за перелёт в момент $t_{видимость} = t_{реакция}$ $-$ непрерывное время видимости цели равно времени реакции. Выбор вектора импульса $\boldsymbol{d \upsilon}$ может быть найдено методом пристрелки, методом дифференциальной эволюции, методом спектрального пучка, или совокупностью этих методов.

Непрерывное управление включает в себя методы:
- ПД-регулятор
- линейно-квадратичный регулятор

### 1
Стандартный ***ПД-регулятор*** может быть представлен в виде $$\boldsymbol{u} = - k_r \boldsymbol{\Delta r} - k_v \dot{\boldsymbol{\Delta r}},$$но он не учитывает особенностей динамической системы. ПД-регулятор может быть ляпуновским управлением, если модифицировать стандартную модель ПД-регулятора. В ОСК на аппарат действует ускорение Хилла-Клохесси-Уилтшира и управление:
$${\boldsymbol{a}} = {\boldsymbol{a}}^{ОСК} + \boldsymbol{u}= \begin{bmatrix}  
-2 \omega_{ХКУ} \dot{z} \\  
- \omega_{ХКУ}^2 y\\  
2\omega_{ХКУ} \dot{x} + 3 \omega_{ХКУ}^2z 
\end{bmatrix} + \boldsymbol{u},$$
где $x$, $y$, $z$ - компоненты вектора $\boldsymbol{r}$ аппарата в ОСК. С помощью этого выражается ускорение аппарата в ССК:
$$\boldsymbol{a}^{ССК} = \boldsymbol{S} \boldsymbol{a}^{ОСК} - [\boldsymbol{S}\boldsymbol{\varepsilon} \times \boldsymbol{r}_1'] - [\boldsymbol{S}\boldsymbol{\omega} \times \boldsymbol{S} \boldsymbol{\omega} \times \boldsymbol{r}_1' ] - \boldsymbol{S} \boldsymbol{A}^{ОСК} - 2[\boldsymbol{S}\boldsymbol{\omega} \times \boldsymbol{S}\boldsymbol{\upsilon}],$$
где $\boldsymbol{\varepsilon}=\dot{\boldsymbol{\omega}}$, $\boldsymbol{r}_1' = \boldsymbol{r}_1  - \boldsymbol{r}_с$, $\boldsymbol{A}^{ОСК}$ - ускорение ХКУ для станции. Модифицированный ПД-регулятор будет иметь вид:
$$\boldsymbol{u} = - k_r \boldsymbol{\Delta r} - k_v \dot{\boldsymbol{\Delta r}}-{\boldsymbol{a}}^{ОСК} + {\boldsymbol{A}}^{ОСК}  + \boldsymbol{S}^T([\boldsymbol{S}\boldsymbol{\varepsilon} \times \boldsymbol{r}_1'] + [\boldsymbol{S}\boldsymbol{\omega} \times \boldsymbol{S}\boldsymbol{\omega} \times \boldsymbol{r}_1'  ] + 2[\boldsymbol{S}\boldsymbol{\omega} \times \boldsymbol{S}\boldsymbol{\upsilon}]), $$
где ${\Delta r}$ $-$ расстояние между аппаратом и целью.
Вследствие полёта близ сложных конструкций критически важна ликвидация колебательной моды управляемого движения. Во время работы модифицированного ПД-регулятора движение относительно ССК подчиняется уравнению:
$$\ddot{\boldsymbol{\Delta r}} + k_{\upsilon} \dot{\boldsymbol{\Delta r}} + k_r \boldsymbol{\Delta r} = 0. \hskip 10px (1)$$
Общее решение уравнения $(1)$ есть 
$$x = c_1 e^{\omega_1 t} + c_2 e^{\omega_2 t},$$
$$\omega_{1,2} = -\frac{k_{\upsilon}}{2} \pm \sqrt{\frac{k^2_{\upsilon}}{4} - k_r}.$$
Отсюда следует, что для безколебательного относительного движения необходимо соотношение коэффициентов: $k_{\upsilon} = 2 \sqrt{k_{r}}$.

### 2
Для ***линейно-квадратичного*** регулятора необходимо описать движение в системе координат целевой точки конструкции ЦСК в виде:
$$\dot{\boldsymbol{X}}= \boldsymbol{A} \boldsymbol{X} + \boldsymbol{X}_{нелин},$$
где $\boldsymbol{X} = [\boldsymbol{r}^T_{ССК}, \boldsymbol{\upsilon}^T_{ССК}]^T \in \mathbb{R}^6$ $-$ вектор состояния, $\boldsymbol{A} \boldsymbol{X}$ $-$ линейная правая часть системы уравнений, $\boldsymbol{X}_{нелин}$ $-$ нелинейная. Ускорение аппарата в ССК:
$$\boldsymbol{a}^{ССК} = \begin{bmatrix}  
-2 \omega_{ХКУ} \dot{z} \\  
- \omega_{ХКУ}^2 y\\  
2\omega_{ХКУ} \dot{x} + 3 \omega_{ХКУ}^2z 
\end{bmatrix} - \begin{bmatrix}  
-2 \omega_{ХКУ} \dot{Z} \\  
- \omega_{ХКУ}^2 Y\\  
2\omega_{ХКУ} \dot{X} + 3 \omega_{ХКУ}^2Z 
\end{bmatrix} - [\boldsymbol{\varepsilon} \times \boldsymbol{r}'_1] - [\boldsymbol{\omega} \times \boldsymbol{\omega} \times \boldsymbol{r}'_1] - 2 [\boldsymbol{\omega} \times \boldsymbol{\upsilon}^{ССК}],$$
где $X$, $Y$, $Z$ $-$ компоненты вектора центра масс.
$$2 [\boldsymbol{\omega} \times \boldsymbol{\upsilon}^{ССК}] = \begin{bmatrix}  
2 \omega_y \upsilon_z^{ССК} - 2 \omega_z \upsilon_y^{ССК} \\  
2 \omega_z \upsilon_x^{ССК} - 2 \omega_x \upsilon_z^{ССК} \\  
2 \omega_x \upsilon_y^{ССК} - 2 \omega_y \upsilon_x^{ССК}
\end{bmatrix};$$
$$\begin{bmatrix}  
x\\  
y\\  
z 
\end{bmatrix} = \boldsymbol{r} = \boldsymbol{R} + \boldsymbol{S}^T(\boldsymbol{r}^{ССК} + \boldsymbol{r}_1'), \hskip 10px \begin{bmatrix} 
\dot{x}\\  
\dot{y}\\  
\dot{z} 
\end{bmatrix} = \boldsymbol{\upsilon} = \boldsymbol{V} + [\boldsymbol{\omega} \times \boldsymbol{r}'_1] + \boldsymbol{S}^T \boldsymbol{\upsilon}^{ССК}, $$
Члены $[\boldsymbol{\omega} \times \boldsymbol{\omega} \times \boldsymbol{r}_1]$, $[\boldsymbol{\varepsilon} \times \boldsymbol{r}_1]$ входят в нелинейную состовляющую правой части системы уравнений. 
$$\boldsymbol{a}^{ССК} = \begin{bmatrix}  
-2 \omega_{ХКУ} (\omega_x r'_{1y} - \omega_y r'_{1x} + s_{13} \upsilon^{ССК}_x + s_{23} \upsilon^{ССК}_y + s_{33} \upsilon^{ССК}_z) \\  
-\omega_{ХКУ}^2 (s_{21}x^{ССК} + s_{22}y^{ССК} + s_{23}z^{ССК} + )\\  
2 \omega_{ХКУ} (\omega_y r'_{1z} - \omega_z r'_{1y} + s_{11} \upsilon^{ССК}_x + s_{21} \upsilon^{ССК}_y + s_{31} \upsilon^{ССК}_z)
\end{bmatrix}$$