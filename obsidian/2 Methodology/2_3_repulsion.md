Движение станции и аппарата рассматривается относительно общего центра масс. Для описания орбитального движения в ОСК введены векторы центра масс аппарата и станции $\boldsymbol{r}$ и $\boldsymbol{R}$, поступательные скорости аппарата и станции $\boldsymbol{\upsilon}$ и $\boldsymbol{V}$ соответственно, угловая скорость станции $\boldsymbol{\omega}$. В стадии удерживания аппарата на станции аппарат рассматривается как часть станции с точки зрения механики. 

Параметры положения $\boldsymbol{r}$, $\boldsymbol{R}$ и параметры движения $\boldsymbol{\upsilon}$, $\boldsymbol{V}$, $\boldsymbol{\omega}$ меняются мгновенно при отталкивании от станции или при присоединении к ней.  

Обозначим векторы центра масс конструкции до отталкивания и после как $\boldsymbol{r}_{c}^{p}$ and $\boldsymbol{r}_{c}$ соответственно, в ССК. Тогда $\boldsymbol{r}_c^{rel} = \boldsymbol{r}_{c} - \boldsymbol{r}_{c}^{p}$ есть изменение центра масс конструкции. Аналогично вычисляются тензоры инерции до отталкивания и после $\boldsymbol{J}$ и $\boldsymbol{J}^p$ соответственно, в ССК. В данном случае, при подсчёте $\boldsymbol{r}_{c}$ и $\boldsymbol{J}$ не учитываются аппарат и переносимые балки, в отличии от $\boldsymbol{r}_{c}^{p}$ и $\boldsymbol{J}^p$. Верхний индекс $\{ {\ast}^p \}$ имеет смысл предыдущего состояния (previous), в данном случае - в состоянии удерживания. В соответствии с отталкиванием изменятся и векторы орбитального движения:
$$\boldsymbol{R} = \boldsymbol{R}^{p} + \boldsymbol{S}^T \boldsymbol{r}_c^{rel},$$
$$\boldsymbol{r} = \boldsymbol{R}^p + \boldsymbol{S}^T (\boldsymbol{r}_0 - \boldsymbol{r}_c^p),$$
где $\boldsymbol{R}^p$ - вектор центра масс станции до отталкивания в ОСК, а $\boldsymbol{r}_0$ - радиус-вектор центра масс аппарата до отталкивания в ССК. 

Изменение величин $\boldsymbol{\upsilon}$, $\boldsymbol{V}$, $\boldsymbol{\omega}$ согласовано с законами сохранения импульса и момента импульса:
$$(m_{a}^+ + M^-)\boldsymbol{V}^p = m_{a}^+ \boldsymbol{\upsilon} + M^- \boldsymbol{V}, (1)$$
$$(m_a^+ + M^-)[R^p \times V^p] + \boldsymbol{J}^p \boldsymbol{\omega}^p =  M^-[R \times V] + m_a^+ [\boldsymbol{r} \times \boldsymbol{\upsilon}] + \boldsymbol{J} \boldsymbol{\omega}, (2)$$
где $\boldsymbol{J}^p$ и $\boldsymbol{J}$ - тензоры инерции станции до и после отталкивания, $m_{a}^+$ $-$ масса аппарата с переносимым грузом, $M^-$ $-$ масса станции без учёта переносимого груза. Эти уравнения справедливы для инерциальной системы координат. Первое уравнение $(1)$ справедливо, поскольку оно эквивалентно уравнению в ИСК: 
$$(m_{a}^+ + M^-)(\boldsymbol{V}^p + \boldsymbol{V}_{орб}) = m_{a}^+ (\boldsymbol{\upsilon} + \boldsymbol{V}_{орб}) + M^- (\boldsymbol{V} + \boldsymbol{V}_{орб}), $$
где $\boldsymbol{V}_{орб}$ - поступательная скорость центра ОСК относительно ИСК. Аналогично, справедливо для связной системы координат. Второе уравнение $(2)$ справедливо по теореме Кёнига: 
$$\boldsymbol{L} = \boldsymbol{L}_{цм} + \boldsymbol{L}',$$
где $\boldsymbol{L}$ $-$ момент импульса системы тел, $\boldsymbol{L}_{цм}$ $-$ момент импульса центра масс системы тел, $\boldsymbol{L}'$ $-$ момент импульса тел относительно общего центра масс.

Исходя из постановки задачи и ограничений, скорость отталкивания $\boldsymbol{\upsilon}_0$ должна быть определена в ССК. Справедливы соотношения:
$$\boldsymbol{\upsilon} = \boldsymbol{\upsilon}_0 + [\boldsymbol{\omega}^p \times (\boldsymbol{r}_0 - \boldsymbol{r}_c^p)] + \boldsymbol{V}^p,$$
$$\boldsymbol{V} = \boldsymbol{V}_0 + [\boldsymbol{\omega}^p \times (\boldsymbol{r}_c - \boldsymbol{r}_c^p)] + \boldsymbol{V}^p,$$
$$\boldsymbol{V}_0 = - \boldsymbol{\upsilon}_0 \frac{m_a^+}{M^-}.$$
(?тут точно не надо делать преобразование $v_0$ в $S^T v_0$?)
$$\boldsymbol{\omega} = \boldsymbol{J} ^{-1}((m_a^+ + M^-)[R^p \times V^p] + \boldsymbol{J}^p \boldsymbol{\omega}^p -  M^-[R \times V] - m_a^+ [\boldsymbol{r} \times \boldsymbol{\upsilon}] ) $$
В ***момент захвата аппарата*** за конструкцию величины $\boldsymbol{R}$, $\boldsymbol{V}$, $\boldsymbol{\omega}$ изменяются по тем же законам. Так как после захвата движение аппарата полностью определяется движением станции, исследование величин $\boldsymbol{r}$, $\boldsymbol{\upsilon}$ опускается. Теперь верхний индекс $\{ {\ast}^p \}$ обозначает величины при отсутствии механического взаимодействия аппарата и станции:
$$\boldsymbol{R} = \boldsymbol{R}^{p} + \boldsymbol{S}^T \boldsymbol{r}_c^{rel},$$
$$\boldsymbol{V} = \frac{m_{a}^+ \boldsymbol{\upsilon}^p + M^- \boldsymbol{V}^p}{(m_{a}^+ + M^-)},$$
$$(m_a^+ + M^-)[R \times V] + \boldsymbol{J} \boldsymbol{\omega} =  M^-[R^p \times V^p] + m_a^+ [\boldsymbol{r}^p \times \boldsymbol{\upsilon}^p] + \boldsymbol{J}^p \boldsymbol{\omega}^p, $$
$$  \boldsymbol{\omega} = \boldsymbol{J}^{-1}( M^-[R^p \times V^p] + m_a^+ [\boldsymbol{r}^p \times \boldsymbol{\upsilon}^p] + \boldsymbol{J}^p \boldsymbol{\omega}^p - (m_a^+ + M^-)[R \times V]).$$
