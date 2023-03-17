

Определим угловые скорости: $\boldsymbol{\omega}$ $-$ вращение ССК относительно ОСК, $\boldsymbol{\Omega}$ $-$ вращение ССК относительно ИСК, и $\boldsymbol{\omega}_{ХКУ}$ $-$ вращение ОСК относительно ИСК. Отметим, что $\boldsymbol{\Omega} = \boldsymbol{\omega} + \boldsymbol{\omega}_{ХКУ}$. Динамическое уравнение Эйлера в ССК описывает вращательное движение тела:
$$\boldsymbol{J} \dot{\boldsymbol{\Omega}} + [\boldsymbol{\Omega} \times \boldsymbol{J} \boldsymbol{\Omega}] = \boldsymbol{M}_{external},$$
где $\boldsymbol{M}_{external}$ $-$ крутящий момент внешних сил, $\boldsymbol{J}$ $-$ тензор инерции станции в ССК, $\boldsymbol{\Omega}$ $-$ угловая скорость станции. Векторы $\boldsymbol{M}_{external}$, $\boldsymbol{\Omega}$ определены в ИСК, подставлены в уравнение в проекции на ССК. В данной работе в момент внешних сил входит только гравитационный момент:
$$\boldsymbol{M}_{grav} = 3 \frac{\mu}{|R|^5} \boldsymbol{R} \times \boldsymbol{J} \boldsymbol{R},$$
где $\boldsymbol{R}$ $-$ вектор между центрами масс Земли и станции в проекции на ССК. Очевидно, можно без потери точности переопределить $\boldsymbol{R}$ как вектор от центра масс Земли до ОСК. В проекции на ОСК вектор $\boldsymbol{R}$ будет иметь только $z$-составляющую. Тогда уравнение угловой скорости запишется в виде:
$$\dot{\boldsymbol{\Omega}} = \boldsymbol{J}^{-1} (3 \frac{\mu}{R^5} [\boldsymbol{S} \begin{bmatrix}0\\0\\R\end{bmatrix} \times \boldsymbol{J} \boldsymbol{S} \begin{bmatrix}0\\0\\R\end{bmatrix}] - [\boldsymbol{\Omega} \times \boldsymbol{J} \boldsymbol{\Omega}]).   $$
В данном уравнении переменные величины: $\boldsymbol{\Omega}$ и $\boldsymbol{S}$. Вращение матриц описывается единичными кватернионами $||\boldsymbol{\Lambda}|| = 1$, или версорами. Для поворота системы координат относительно другой на угол $\theta$ относительно вектора $\boldsymbol{\lambda}=[\lambda_x, \lambda_y, \lambda_z]$ версор задаётся как:
$$\boldsymbol{\Lambda} = [cos \frac{\theta}{2}, \lambda_x sin \frac{\theta}{2}, \lambda_y sin \frac{\theta}{2}, \lambda_z sin \frac{\theta}{2} ].$$
Версор $\boldsymbol{\Lambda}$ связан с угловой скоростью $\boldsymbol{\omega}$ кинематическим уравнением [[quaternions.pdf | [n] ]]:
$$\dot{\boldsymbol{\Lambda}} = \frac{1}{2} \boldsymbol{\omega} \circ \boldsymbol{\Lambda}. $$
Далее, матрица поворота $\boldsymbol{A}$ во вращающуюся систему координат однозначно определяется через компоненты версора $\boldsymbol{\Lambda}_A = [\omega, x, y, z]$:
$$\boldsymbol{A} = \begin{bmatrix}  
1 - 2y^2 - 2z^2 & 2xy + 2z \omega & 2xz - 2y \omega\\  
2xy - 2z \omega & 1 - 2x^2 - 2z^2 & 2yz + 2x \omega\\  
2xz + 2y \omega & 2yz - 2x \omega & 1 - 2x^2 - 2y^2  
\end{bmatrix}.$$
Для описания движения при использовании трёх систем координат необходимо ввести версоры $\boldsymbol{\Lambda}_U$ и $\boldsymbol{\Lambda}_A$, определяющие изменение матриц $\boldsymbol{U}$ и $\boldsymbol{A}$ с угловыми скоростями $\boldsymbol{\omega}_{ХКУ}$ и $\boldsymbol{\Omega}$ соответственно, а вследствие определяющие матрицу $\boldsymbol{S} = \boldsymbol{A}\boldsymbol{U}^T$. 