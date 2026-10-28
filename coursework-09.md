# Псевдоевклидово пространство и локальная симметрия сферы

Дифференциальная геометрия и основы тензорного исчисления. Вариант №9.

## Введение

Предмет работы — преобразования, сохраняющие метрику, и связанные с ними инварианты. Вариант №9 предусматривает две задачи: изложить определение и основные свойства псевдоевклидова пространства, связав их со специальной теорией относительности; доказать локальную симметричность обычной сферы $S^2$ и вычислить её гауссову и скалярную кривизны.

В пространстве Минковского условие сохранения интервала определяет преобразования Лоренца. Из них выводятся соотношения для собственного времени, длины и скорости; запись через быстроту раскрывает геометрический смысл и композицию преобразований. Численные модели и экспериментальное сопоставление показывают масштаб полученных следствий.

На сфере метрика задаётся в касательных плоскостях. Её координатный расчёт приводит к тензору кривизны, а явная изометрия обращает геодезические через заданную точку. Для единичной сферы получены $K=1$, $\operatorname{Scal}=2$ и доказана локальная симметричность. Общим языком обеих задач служат метрика и её изометрии; сигнатура и кривизна характеризуют разные свойства пространств.

## 1. Псевдоевклидово пространство

### 1.1. Билинейная форма, сигнатура и изотропные направления

Псевдоевклидово линейное пространство — пара $(V,g)$, где $V$ — вещественное линейное пространство размерности $n$, а $g$ — симметричная невырожденная билинейная форма [1, лекция 12а]. Невырожденность означает $g(u,v)=0$ для всех $v\in V\Rightarrow u=0$; ортогональность определяется равенством $g(u,v)=0$. В аффинной модели форма задана на пространстве переносов, и квадрат интервала равен $s^2(P,Q)=g(\overrightarrow{PQ},\overrightarrow{PQ})$.

В базисе $e=(e_1,\ldots,e_n)$, при $u=ex$, $v=ey$, имеем [1, лекции 5 и 11]
```math
G=(g(e_i,e_j)),\qquad g(u,v)=x^{\mathsf T}Gy,\qquad
G^{\mathsf T}=G,\qquad \det G\ne0.
```
При смене базиса $e'=eC$, где $C$ обратима,
```math
x=Cx',\qquad G'=C^{\mathsf T}GC.
```
По закону инерции Сильвестра существует базис, в котором [1, лекции 11–12]
```math
G=J_{p,q}=\operatorname{diag}(I_p,-I_q),\qquad p+q=n,\qquad
g(u,u)=\sum_{i=1}^{p}(x^i)^2-\sum_{i=p+1}^{n}(x^i)^2.
```
Пара $(p,q)$ не зависит от базиса и называется сигнатурой; такой базис — псевдоортонормированным. При $q=0$ форма положительно определена, при $p,q>0$ принимает значения обоих знаков на квадратах векторов.

Ненулевой вектор $u$ изотропен, если $g(u,u)=0$. При $G=\operatorname{diag}(1,-1)$ вектор $u=(1,1)$ изотропен, но $g(u,(1,0))=1$: изотропность не означает вырожденности формы. Изотропные направления здесь — $y=\pm x$, а $\sqrt{|g(u,u)|}$ не является нормой [1, лекция 12а].

### 1.2. Ортогональные дополнения и изометрии

Для подпространства $W\subset V$ [1, лекция 12а, с.163]
```math
W^\perp=\{v\in V:\ g(v,w)=0\ \forall w\in W\},\qquad
\dim W^\perp=n-\dim W,
```
```math
V=W\oplus W^\perp\quad\Longleftrightarrow\quad
W\cap W^\perp=\{0\}\quad\Longleftrightarrow\quad g|_W\text{ невырожденно}.
```
Для изотропной прямой $W=\mathbb R(1,1)$ в предыдущем примере $W^\perp=W$, поэтому прямой суммы нет.

Линейная изометрия сохраняет форму; её матричное условие имеет вид
```math
g(Lu,Lv)=g(u,v)\ \forall u,v
\quad\Longleftrightarrow\quad L^{\mathsf T}GL=G,
\qquad \det L=\pm1.
```
В псевдоортонормированном базисе такие преобразования образуют группу [1, лекция 12а]
```math
O(p,q)=\{L\in GL(n,\mathbb R):L^{\mathsf T}J_{p,q}L=J_{p,q}\}.
```
Здесь $GL(n,\mathbb R)$ — группа обратимых вещественных матриц. Аффинные изометрии $F(x)=Lx+b$ сохраняют интервалы между точками. Для сигнатуры $(1,3)$ получается пространство Минковского: одно положительное направление связано со временем, три отрицательных — с пространственными координатами. Его линейные изометрии образуют группу Лоренца $O(1,3)$; условие сохранения формы становится основой геометрического описания СТО.

## 2. Геометрия пространства Минковского и СТО

Условие изометрии определяет допустимые переходы между инерциальными системами. Сохраняемый интервал позволяет вывести преобразования координат, а затем сопоставить измерения времени, длины и скорости.

### 2.1. Интервал и световые направления

Событие задаётся координатами $x^0=ct$, $x^1=x$, $x^2=y$, $x^3=z$, где $c$ — скорость света в вакууме. События — точки аффинного пространства; разности координат — векторы разделения с компонентами размерности длины. Форма Минковского равна [3, §§1–2, 6]
```math
\eta=\operatorname{diag}(1,-1,-1,-1),\qquad
\eta(A,B)=A^0B^0-A^1B^1-A^2B^2-A^3B^3.
```
Для двух событий и для бесконечно близких событий соответственно [3, §2]
```math
\Delta s^2=c^2(\Delta t)^2-|\Delta\mathbf r|^2,
\qquad ds^2=c^2dt^2-dx^2-dy^2-dz^2.
```
Ненулевое разделение времениподобно при $\Delta s^2>0$, светоподобно при $\Delta s^2=0$ и пространственноподобно при $\Delta s^2<0$. Нулевой интервал задаёт световой конус; времениподобные направления лежат внутри его полостей. Переходы, сохраняющие интервал, сохраняют эту классификацию.

### 2.2. Преобразование Лоренца и гиперболический поворот

Пусть $K'$ движется относительно $K$ со скоростью $v$ вдоль оси $x$; оси $x,x'$ сонаправлены. Начала совпадают при $t=t'=0$; $y'=y$, $z'=z$. Аффинность перехода принимается в модели инерциальных систем; совпадение начал устраняет перенос [3, §4]. Положим $X=ct$, $\beta=v/c$. Линейная функция $x'$ обращается в нуль на мировой линии начала $K'$, то есть при $x=\beta X$, поэтому
```math
x'=A(x-\beta X),\qquad X'=BX+Dx.
```
Инвариантность интервала даёт
```math
(BX+Dx)^2-A^2(x-\beta X)^2=X^2-x^2
```
Сравнивая коэффициенты, получаем
```math
B^2-A^2\beta^2=1,\qquad D^2-A^2=-1,\qquad BD+A^2\beta=0.
```
Из первого равенства $B\ne0$. Исключая $D=-A^2\beta/B$ из второго, находим $A^2(B^2-A^2\beta^2)=B^2$, откуда $A^2=B^2$; также $A\ne0$. Нормировка $A(0)=B(0)=1$ и непрерывность при $|v|<c$ фиксируют положительную ветвь $A=B>0$, сохраняющую направление времени и сонаправленность осей. Следовательно,
```math
D=-A\beta,\qquad A^2(1-\beta^2)=1,
\qquad A=\gamma=\frac1{\sqrt{1-\beta^2}},\qquad |v|<c.
```
**Преобразование Лоренца.**
```math
x'=\gamma(x-vt),\qquad
t'=\gamma\left(t-\frac{vx}{c^2}\right),\qquad y'=y,\qquad z'=z.
```
Проверка: $X'^2-x'^2=\gamma^2[(X-\beta x)^2-(x-\beta X)^2]=X^2-x^2$. Вместе с $y'=y$, $z'=z$ это равносильно матричному условию $\Lambda^{\mathsf T}\eta\Lambda=\eta$. Обратный переход получается заменой $v$ на $-v$.

На рисунке 1 показано сечение $y=z=0$, $ct\ge0$, при $v=0{,}6c$. Оси $K'$ определяются условиями $x'=0\Leftrightarrow x=\beta ct$ и $t'=0\Leftrightarrow ct=\beta x$; световые прямые остаются границами конуса. Гипербола $c^2t^2-x^2=\mathrm{const}>0$ объединяет события с одинаковым интервалом от начала. В условных единицах $P:(x,ct)=(3/4,5/4)\mapsto(x',ct')=(0,1)$: в обеих системах $s^2=1$.

![Диаграмма Минковского](figures/01-minkowski.svg)

Рисунок 1 — Диаграмма Минковского

Геометрический смысл перехода раскрывает параметр быстроты $\psi$: $\beta=\tanh\psi$, $\cosh\psi=\gamma$, $\sinh\psi=\gamma\beta$. В порядке координат $(X,x)$ получаем
```math
\begin{gathered}
\begin{pmatrix}X'\\x'\end{pmatrix}
=\Lambda(\psi)\begin{pmatrix}X\\x\end{pmatrix},
\qquad
\Lambda(\psi)=
\begin{pmatrix}
\cosh\psi&-\sinh\psi\\
-\sinh\psi&\cosh\psi
\end{pmatrix},\\[1mm]
\Lambda(\psi_2)\Lambda(\psi_1)=\Lambda(\psi_1+\psi_2).
\end{gathered}
```
Матрица описывает переход $K\to K'$, а обратный соответствует $-\psi$ [3, §4, с. 23–24]. На рисунке 2 обычный поворот сохраняет окружность $x_1^2+x_2^2=1$, гиперболический — будущую ветвь $X^2-x^2=1$. Образ $e_0=(1,0)$ равен $(\cosh\psi,-\sinh\psi)$; при горизонтальной оси $x$ и вертикальной $X$ точки с $\psi>0$ лежат слева. Знак согласован со скоростью $-v$ покоящегося в $K$ тела при переходе в $K'$. Отмечены $\psi=0,\ln2,\ln3$; оси фиксированы, масштабы координат равны. Быстрота не измеряется обычным углом между лучами, а форма гиперболы не означает кривизны пространства Минковского.

Формулы сложения гиперболических функций дают последнее равенство в (15): при последовательных сонаправленных бустах складываются быстроты. Если $K_1$ движется относительно $K$ со скоростью $v_1$, а $K_2$ относительно $K_1$ — со скоростью $v_2$ вдоль той же оси, то $\psi=\psi_1+\psi_2$. Из $\tanh(\psi_1+\psi_2)=(\beta_1+\beta_2)/(1+\beta_1\beta_2)$ следует $v=(v_1+v_2)/(1+v_1v_2/c^2)$. Таким образом, композиция изометрий задаёт закон сложения скоростей и сохраняет условие $|v|<c$.

![Орбиты обычного и гиперболического поворотов](figures/07-rotation-orbits.svg)

Рисунок 2 — Орбиты обычного и гиперболического поворотов

### 2.3. Собственное время и измерение длины

На направленной в будущее мировой линии $dt>0$, $\mathbf u=d\mathbf r/dt$, $|\mathbf u|<c$. Из инвариантного соотношения $ds^2=c^2dt^2(1-|\mathbf u|^2/c^2)=c^2d\tau^2$ получаем [3, §3]
```math
\begin{aligned}
d\tau&=dt\sqrt{1-\frac{|\mathbf u|^2}{c^2}},\\
\tau_2-\tau_1&=\int_{t_1}^{t_2}\sqrt{1-\frac{|\mathbf u(t)|^2}{c^2}}\,dt.
\end{aligned}
```
Промежуток $\Delta\tau$ относится к двум событиям на мировой линии одних часов; $\Delta t$ — к координатам тех же событий в выбранной инерциальной системе. Для фиксированной линии собственное время инвариантно; для разных линий между теми же событиями оно может различаться. При постоянной скорости $\Delta t=\gamma_u\Delta\tau$, $\gamma_u=(1-|\mathbf u|^2/c^2)^{-1/2}$. На светоподобной линии $ds^2=0$ собственное время не служит параметром.

Световые часы на рисунке 3 покоятся в $K'$: зеркала разделены расстоянием $h$ поперёк движения. События $A,B,C$ — испускание, отражение и возвращение света; полный цикл длится $\Delta\tau$ в $K'$ и $\Delta t$ в $K$. При движении вдоль $x$ высота сохраняется ($y'=y$). Для половины цикла
```math
\left(\frac{c\Delta t}{2}\right)^2=h^2+\left(\frac{v\Delta t}{2}\right)^2,
\qquad \Delta\tau=\frac{2h}{c}.
```
Отсюда $\Delta t=\gamma\Delta\tau$: между теми же событиями начала и конца цикла свет в $K$ проходит больший путь со скоростью $c$. При $v=0{,}6c$, $\Delta\tau=1$ нс получаем $\Delta t=1{,}25$ нс. Схема дана в плоскости $xy$; положения зеркал справа относятся к разным моментам прохождения света.

![Световые часы](figures/02-light-clock.svg)

Рисунок 3 — Световые часы в двух инерциальных системах

Для измерения длины существенно условие одновременности регистрации концов. Преобразование времени даёт
```math
\Delta t=0\quad\Longrightarrow\quad
\Delta t'=-\gamma\frac{v\Delta x}{c^2}.
```
При $v\Delta x\ne0$ события, одновременные в $K$, не одновременны в $K'$. Пусть стержень покоится в $K'$, направлен вдоль $x'$ и имеет собственную длину $L_0$. В $K$ измерение требует $\Delta t=0$, $\Delta x=L$. Неподвижность концов в $K'$ означает $\Delta x'=L_0$ и для этой пары событий [3, §4], поэтому
```math
L_0=\Delta x'=\gamma(\Delta x-v\Delta t)=\gamma L,
\qquad L=\frac{L_0}{\gamma}.
```
На рисунке 4 пара $A,B$ соответствует одному $t$, а $A,C$ — одному $t'$ (прямая $t'=0$). При $L_0=1$ м, $v=0{,}6c$ имеем $L=0{,}8$ м. Сокращение длины следует из сопоставления разных пар событий при сохранении интервала.

![Измерение длины движущегося стержня](figures/03-length.svg)

Рисунок 4 — Измерение длины движущегося стержня

### 2.4. Скорости и четырёхвекторы

Для скорости частицы $\mathbf u=(u_x,u_y,u_z)$ дифференцирование преобразования Лоренца даёт $dt'=\gamma(1-vu_x/c^2)dt$. Деление пространственных дифференциалов на $dt'$ приводит к формулам [3, §5]
```math
u_x'=\frac{u_x-v}{1-vu_x/c^2},\qquad
u_y'=\frac{u_y}{\gamma(1-vu_x/c^2)},\qquad
u_z'=\frac{u_z}{\gamma(1-vu_x/c^2)}.
```
При $|\mathbf u|<c$, $|v|<c$ знаменатель положителен; из инвариантности интервала следует $|\mathbf u'|<c$. Для $\mathbf u=(v,0,0)$ частица покоится в $K'$, а для луча $u_x=c$ имеем $u_x'=c$. При $u_x=0{,}8c$, $v=0{,}6c$, $u_y=u_z=0$ получаем $u_x'=5c/13\approx0{,}385c$.

Инвариантное описание удобно записать через четырёхвекторы $A'^\mu=\Lambda^\mu{}_{\nu}A^\nu$. Греческие индексы принимают значения $0,1,2,3$, пространственные $i$ — $1,2,3$; по повторяющимся верхнему и нижнему индексам суммируем. Условие изометрии даёт [3, §6]
```math
\eta_{\mu\nu}A'^\mu B'^\nu=A_\mu B^\mu,
\qquad A_\mu=\eta_{\mu\nu}A^\nu.
```
Следовательно, $A_0=A^0$, $A_i=-A^i$. Четырёхскорость $U^\mu=dx^\mu/d\tau=\gamma_u(c,u_x,u_y,u_z)$ имеет инвариантный квадрат $\eta(U,U)=c^2$ [3, §7].

### 2.5. Математические следствия и модели

**Коэффициент Лоренца.** При постоянном модуле скорости отношение лабораторного промежутка времени к собственному определяется функцией $\gamma(\beta)=(1-\beta^2)^{-1/2}$, $0\le\beta<1$:
```math
\gamma(0{,}6)=1{,}25,\qquad
\gamma(0{,}9)\approx2{,}29,\qquad
\gamma(0{,}99)\approx7{,}09.
```
![Коэффициент Лоренца в зависимости от скорости](figures/05-lorentz-factor.svg)

Рисунок 5 — Зависимость коэффициента Лоренца от скорости

График на рисунке 5 показывает переход от $\gamma\approx1+\beta^2/2$ при малых $\beta$ к неограниченному росту при $\beta\to1-0$.

**Модель распада частиц.** Пусть доля нераспавшихся частиц в собственном времени равна $e^{-\tau/\tau_0}$. При постоянной скорости $0<v<c$ из $d\tau=dt/\gamma$ и $s=vt$ следует $\tau=s/(\gamma v)$, поэтому
```math
\frac{N(s)}{N_0}=\exp\left(-\frac{s}{\gamma v\tau_0}\right),
\qquad
\ln\frac{N(s)}{N_0}=-\frac{s}{\gamma v\tau_0}.
```
Здесь $\tau_0$ — среднее, а не предельное время жизни в системе покоя. Для начального числа $N_0$ величина $N(s)$ — ожидаемое число нераспавшихся частиц; $N(s)/N_0$ — вероятность достижения пути $s$ до распада. В полулогарифмических координатах (23) задаёт прямую с наклоном $-1/(\gamma v\tau_0)$.

Рассмотрим модель движения мюонов с собственным средним временем жизни $\tau_0\approx2{,}197$ мкс [4, таблица мюона, с. 2]. Заданы вертикальный путь $h=10$ км и постоянная скорость $v=0{,}998c$; потери энергии и взаимодействия с атмосферой не учитываются. Тогда
```math
\begin{aligned}
\gamma&\approx15{,}82,\qquad
t_{\text{пр}}=\frac hv\approx33{,}42\ \text{мкс},\qquad
\gamma\tau_0\approx34{,}75\ \text{мкс},\\
\frac{N(h)}{N_0}&\approx e^{-33{,}42/34{,}75}\approx0{,}382,\qquad
\left.\frac{N(h)}{N_0}\right|_{\text{без замедления}}
=e^{-h/(v\tau_0)}\approx2{,}47\cdot10^{-7}.
\end{aligned}
```
В системе мюона $h'=h/\gamma\approx0{,}632$ км, а пролёт занимает $h'/v\approx2{,}11$ мкс собственного времени. На рисунке 6 сопоставлены расчётные зависимости с учётом СТО и без замедления времени; вертикальная шкала логарифмическая.

![Доля нераспавшихся мюонов на пути к поверхности](figures/06-muon-survival.svg)

Рисунок 6 — Доля нераспавшихся мюонов в двух моделях

В выбранной модели учёт собственного времени увеличивает долю нераспавшихся частиц примерно на шесть порядков. График представляет теоретическое сопоставление, а не измерение атмосферного потока.

**Экспериментальное сопоставление.** В опыте CERN для положительных мюонов при $\gamma=29{,}33$ измерено среднее время жизни $(64{,}419\pm0{,}058)$ мкс [5, аннотация]; расчёт $\gamma\tau_0$ даёт около $64{,}44$ мкс. Постоянство модуля скорости на круговой орбите позволяет применить тот же интеграл собственного времени.

**Распространение высокоэнергичных нейтронов.** Для свободного нейтрона $\tau_n\approx878{,}4$ с и $m_nc^2\approx939{,}6$ МэВ [4, таблица нейтрона, с. 1, 3]. При полной энергии $E=10^{18}$ эВ имеем $\gamma=E/(m_nc^2)$ [3, §9]. Для условного расстояния $D=8$ кпк, пренебрегая гравитацией и взаимодействиями в пути, получаем
```math
\begin{aligned}
\gamma&\approx1{,}06\cdot10^9,\qquad
\bar t_n=\gamma\tau_n\approx3{,}0\cdot10^4\ \text{лет},\\
L_n&\approx c\bar t_n\approx9{,}1\ \text{кпк},\qquad
t_D\approx D/c\approx2{,}6\cdot10^4\ \text{лет},\\
\Delta\tau&=t_D/\gamma\approx13\ \text{мин},\qquad
P=\exp(-D/L_n)\approx0{,}41.
\end{aligned}
```
Величины $\bar t_n$ и $L_n$ — средние время жизни в системе Земли и путь до распада; 1 кпк — около 3260 световых лет. Тысячи лет координатного времени соответствуют минутам собственного. Вероятность $P$ относится к модели: в поиске Pierre Auger 2012 года значимого избытка от точечных нейтронных источников не обнаружено [7].

**Релятивистские поправки GPS.** Вклад движения по СТО составляет около $-7{,}1$ мкс/сут относительно наземных часов, гравитационный вклад ОТО — около $+45{,}7$ мкс/сут; сумма равна $+38{,}6$ мкс/сут [6, с. 204–205, табл. 1]. Первый связан с рассмотренным собственным временем, второй требует гравитационной модели и не выводится из плоской метрики Минковского.

## 3. Метрика, связность и кривизна

На сфере метрика задаётся в каждой касательной плоскости; дифференциалы изометрий сохраняют эти формы. Для вычисления кривизны необходима согласованная с метрикой связность; для доказательства локальной симметрии — обращение геодезических через заданную точку.

### 3.1. Метрика и связность

Регулярная взаимно однозначная параметризация $r:U\subset\mathbb R^2\to M$ задаёт координатный участок и базис $\partial_i=r_{x^i}$ касательной плоскости. Далее $i,j,k,l,m=1,2$; по повторяющимся верхнему и нижнему индексам суммируем.

Риманова метрика — гладкая положительно определённая симметричная форма на касательных пространствах. Для поверхности в евклидовом пространстве [2, лекция 11]
```math
g_{ij}=\langle r_{x^i},r_{x^j}\rangle,\qquad
ds^2=g_{ij}\,dx^i dx^j,\qquad
(g^{ij})=G^{-1},\qquad g^{ik}g_{kj}=\delta^i{}_j.
```
Здесь $G=(g_{ij})$, $\delta^i{}_j$ — символ Кронекера. Связность Леви-Чивиты позволяет дифференцировать касательные поля; она единственна при условиях отсутствия кручения и согласованности с метрикой [2, лекция 11]:
```math
\nabla_XY-\nabla_YX=[X,Y],\qquad \nabla g=0,\qquad
\nabla_{\partial_i}\partial_j=\Gamma^k{}_{ij}\partial_k.
```
Ковариантная производная учитывает изменение координатного базиса, поэтому $\nabla g=0$ совместимо с переменными компонентами $g_{ij}$. Она линейна по $X$ над гладкими функциями, по $Y$ над $\mathbb R$ и удовлетворяет правилу $\nabla_X(fY)=X(f)Y+f\nabla_XY$; $[X,Y]$ — коммутатор полей. Для координатных полей $[\partial_i,\partial_j]=0$, откуда $\Gamma^k{}_{ij}=\Gamma^k{}_{ji}$ и [2, с.174–175]
```math
\Gamma^k{}_{ij}=\frac12g^{k\ell}
\left(\partial_i g_{j\ell}+\partial_jg_{i\ell}-\partial_\ell g_{ij}\right).
```
Эта формула определяет связность по коэффициентам метрики и будет применена к сферическим координатам.

### 3.2. Кривизна и геодезическая симметрия

Фиксируем знак кривизны и порядок индексов [2, лекция 2]:
```math
R(X,Y)Z=\nabla_X\nabla_YZ-\nabla_Y\nabla_XZ-\nabla_{[X,Y]}Z,
\qquad R(\partial_i,\partial_j)\partial_k=R^l{}_{kij}\partial_l,
```
В $R^l{}_{kij}$ индексы $i,j$ задают направления дифференцирования, $k$ — аргумент $\partial_k$, $l$ — компоненту результата. В координатном базисе
```math
R^l{}_{kij}=\partial_i\Gamma^l{}_{jk}-\partial_j\Gamma^l{}_{ik}
+\Gamma^m{}_{jk}\Gamma^l{}_{im}-\Gamma^m{}_{ik}\Gamma^l{}_{jm}.
```
Тензор $R$ имеет тип $(1,3)$ и кососимметричен по $i,j$. Его свёртки и гауссова кривизна для независимых $X,Y\in T_pM$ задаются формулами [2, лекции 15–17]
```math
\operatorname{Ric}_{kj}=\sum_iR^i{}_{kij},\qquad
\operatorname{Scal}=\sum_{k,j}g^{kj}\operatorname{Ric}_{kj},\qquad
K=\frac{g(R(X,Y)Y,X)}{g(X,X)g(Y,Y)-g(X,Y)^2}.
```
На поверхности $T_pM$ двумерно, поэтому $K$ не зависит от выбора базиса. В инерциальных координатах Минковского $\Gamma=0$ и $R=0$: постоянная неопределённая форма задаёт плоскую геометрию. Для сферы (30)–(31) позволят вычислить кривизну из её метрики.

Геодезическая определяется уравнением $\nabla_{\dot\gamma}\dot\gamma=0$. Экспонента связывает начальное направление с точкой геодезической:
```math
\exp_p(v)=\gamma_v(1),\qquad \gamma_v(0)=p,\qquad
\dot\gamma_v(0)=v,\qquad (d\exp_p)_0=\operatorname{id}.
```
На достаточно малом шаре $B\subset T_pM$ экспонента — диффеоморфизм на нормальную окрестность [2, лекция 1, с.19–20]. Геодезическая симметрия обращает начальное направление:
```math
s_p(\exp_pv)=\exp_p(-v),\qquad v\in B,\qquad
s_p(p)=p,\qquad (ds_p)_p=-\operatorname{id}.
```
Риманово пространство локально симметрично, если все $s_p$ — локальные изометрии. Для сферы это свойство будет доказано непосредственно; равенство $\nabla R=0$ отдельно проверяется для полного тензора кривизны.

## 4. Геометрия и локальная симметрия сферы

### 4.1. Сфера и её метрика

Вычислим кривизну сферы постоянного радиуса $a>0$ с метрикой, индуцированной евклидовым скалярным произведением; для сферы из задания положим $a=1$.

Поверхность $S_a^2=\{x\in\mathbb R^3:\langle x,x\rangle=a^2\}$ регулярна: для $F(x)=\langle x,x\rangle-a^2$ имеем $\operatorname{grad}F(p)=2p\ne0$. Поэтому
```math
T_pS_a^2=\ker dF_p=\{v:\langle p,v\rangle=0\}=p^\perp,
\qquad g_p(v,w)=\langle v,w\rangle.
```
При $0\ne v\in T_pS_a^2$ получаем $g_p(v,v)=|v|^2>0$, то есть метрика риманова на всей сфере. В координатах $x^1=\theta$, $x^2=\varphi$ возьмём
```math
r(\theta,\varphi)=a(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta),
\qquad 0<\theta<\pi,\quad \alpha<\varphi<\alpha+2\pi.
```
Карта исключает полюса и меридианный разрез $\varphi=\alpha$; в указанной области она взаимно однозначна. Дифференцирование даёт
```math
\begin{aligned}
r_\theta&=a(\cos\theta\cos\varphi,\cos\theta\sin\varphi,-\sin\theta),\\
r_\varphi&=a(-\sin\theta\sin\varphi,\sin\theta\cos\varphi,0).
\end{aligned}
```
Скалярные произведения координатных векторов равны
```math
\begin{aligned}
g_{\theta\theta}
&=a^2\bigl(\cos^2\theta(\cos^2\varphi+\sin^2\varphi)+\sin^2\theta\bigr)=a^2,\\
g_{\theta\varphi}
&=a^2\cos\theta\sin\theta(-\cos\varphi\sin\varphi+\sin\varphi\cos\varphi)=0,\\
g_{\varphi\varphi}
&=a^2\sin^2\theta(\sin^2\varphi+\cos^2\varphi)=a^2\sin^2\theta.
\end{aligned}
```
Следовательно,
```math
G=\operatorname{diag}(a^2,a^2\sin^2\theta),\qquad
ds^2=a^2(d\theta^2+\sin^2\theta\,d\varphi^2),
```
```math
\det G=a^4\sin^2\theta>0,\qquad
G^{-1}=\operatorname{diag}\left(a^{-2},(a^2\sin^2\theta)^{-1}\right).
```
При $a>0$, $0<\theta<\pi$ базис невырожден; деления на $a$ и $\sin\theta$ допустимы.

### 4.2. Связность, тензор кривизны и его свёртки

Единственная не тождественно нулевая производная метрики имеет вид
```math
\partial_\theta g_{\varphi\varphi}=2a^2\sin\theta\cos\theta
```
Подстановка в (28) даёт
```math
\begin{aligned}
\Gamma^\theta{}_{\varphi\varphi}
&=-\frac12g^{\theta\theta}\partial_\theta g_{\varphi\varphi}
=-\frac1{2a^2}\,2a^2\sin\theta\cos\theta
=-\sin\theta\cos\theta,\\
\Gamma^\varphi{}_{\theta\varphi}
&=\frac12g^{\varphi\varphi}\partial_\theta g_{\varphi\varphi}
=\frac{2a^2\sin\theta\cos\theta}{2a^2\sin^2\theta}
=\cot\theta.
\end{aligned}
```
По симметрии нижних индексов $\Gamma^\varphi{}_{\varphi\theta}=\cot\theta$; остальные символы равны нулю. Применяя (30), получаем
```math
\begin{aligned}
R^\theta{}_{\varphi\theta\varphi}
&=\partial_\theta\Gamma^\theta{}_{\varphi\varphi}
-\Gamma^\varphi{}_{\theta\varphi}\Gamma^\theta{}_{\varphi\varphi}\\
&=\partial_\theta(-\sin\theta\cos\theta)
-\cot\theta(-\sin\theta\cos\theta)\\
&=\sin^2\theta-\cos^2\theta+\cos^2\theta=\sin^2\theta,\\[2mm]
R^\varphi{}_{\theta\varphi\theta}
&=-\partial_\theta\Gamma^\varphi{}_{\theta\varphi}
-(\Gamma^\varphi{}_{\theta\varphi})^2\\
&=-\partial_\theta(\cot\theta)-\cot^2\theta
=\frac1{\sin^2\theta}-\frac{\cos^2\theta}{\sin^2\theta}=1.
\end{aligned}
```
Оставшиеся компоненты при $(i,j)=(\theta,\varphi)$:
```math
R^\theta{}_{\theta\theta\varphi}=0,\qquad
R^\varphi{}_{\varphi\theta\varphi}=-\partial_\varphi(\cot\theta)=0.
```
При $i=j$ компоненты нулевые, перестановка $i,j$ меняет знак. Тем самым определены все компоненты, в частности
```math
R(\partial_\theta,\partial_\varphi)\partial_\theta=-\partial_\varphi,
\qquad R(\partial_\theta,\partial_\varphi)\partial_\varphi
=\sin^2\theta\,\partial_\theta.
```
Для $X=\partial_\theta$, $Y=\partial_\varphi$ формула (31) даёт
```math
K=\frac{g(\sin^2\theta\,\partial_\theta,\partial_\theta)}
{g_{\theta\theta}g_{\varphi\varphi}-g_{\theta\varphi}^2}
=\frac{a^2\sin^2\theta}{a^4\sin^2\theta}
=\frac1{a^2}.
```
Свёртка $\operatorname{Ric}_{kj}=\sum_iR^i{}_{kij}$ приводит к
```math
\begin{aligned}
\operatorname{Ric}_{\theta\theta}&=R^\varphi{}_{\theta\varphi\theta}=1,\\
\operatorname{Ric}_{\varphi\varphi}&=R^\theta{}_{\varphi\theta\varphi}=\sin^2\theta,\\
\operatorname{Ric}_{\theta\varphi}&=\operatorname{Ric}_{\varphi\theta}=0.
\end{aligned}
```
Следовательно,
```math
(\operatorname{Ric}_{ij})=\operatorname{diag}(1,\sin^2\theta),
\qquad \operatorname{Ric}=\frac1{a^2}g.
```
Сворачивая с обратной метрикой, получаем скалярную кривизну:
```math
\operatorname{Scal}
=g^{\theta\theta}\operatorname{Ric}_{\theta\theta}
+g^{\varphi\varphi}\operatorname{Ric}_{\varphi\varphi}
=\frac1{a^2}+\frac{\sin^2\theta}{a^2\sin^2\theta}
=\frac2{a^2}=2K.
```

### 4.3. Тензор кривизны на всей сфере

При $\theta=0,\pi$ вырождается карта, а не метрика на $p^\perp$. Для каждого $p$ можно выбрать полярную ось и разрез, не содержащие $p$. Карта $\widetilde r=Qr$ при ортогональном $Q$ имеет ту же метрику. Такие карты покрывают сферу и распространяют полученные тензорные результаты на всю $S_a^2$.

Проверим по вычисленным компонентам равенство
```math
R(X,Y)Z=\frac1{a^2}\bigl(g(Y,Z)X-g(X,Z)Y\bigr).
```
При $X=\partial_\theta$, $Y=\partial_\varphi$ правая часть на двух базисных $Z$ равна
```math
\begin{aligned}
Z=\partial_\theta:&\qquad
\frac1{a^2}(g_{\varphi\theta}\partial_\theta-g_{\theta\theta}\partial_\varphi)
=-\partial_\varphi,\\
Z=\partial_\varphi:&\qquad
\frac1{a^2}(g_{\varphi\varphi}\partial_\theta-g_{\theta\varphi}\partial_\varphi)
=\sin^2\theta\,\partial_\theta.
\end{aligned}
```
Это совпадает с (44). Обе части кососимметричны по $X,Y$, при равных базисных $X,Y$ равны нулю; по тензорности равенство верно для любых $X,Y,Z$. Проверка применима в каждой из поворотных карт. В координатах
```math
R^l{}_{kij}=\frac1{a^2}\left(\delta^l{}_i g_{jk}-\delta^l{}_j g_{ik}\right).
```
Проверим ковариантное постоянство этого тензора. Радиус $a$ постоянен и $\nabla g=0$; полная производная имеет вид
```math
\begin{aligned}
(\nabla_WR)(X,Y)Z
&=\nabla_W(R(X,Y)Z)-R(\nabla_WX,Y)Z\\
&\quad-R(X,\nabla_WY)Z-R(X,Y)(\nabla_WZ).
\end{aligned}
```
В (52) учитывается изменение результата $R(X,Y)Z$ и трёх аргументов. После подстановки (49) правило Лейбница сокращает слагаемые с производными $X,Y,Z$; остаются только производные метрики:
```math
(\nabla_WR)(X,Y)Z=\frac1{a^2}\bigl((\nabla_Wg)(Y,Z)X-(\nabla_Wg)(X,Z)Y\bigr)=0.
```
Таким образом, $\nabla R=0$ на всей сфере. Это свойство полного тензора: отдельные координатные компоненты могут быть непостоянны.

### 4.4. Локальная геодезическая симметрия

Ортогональное преобразование $\mathbb R^3$, сохраняющее сферу, сохраняет и её индуцированную метрику [2, с. 381]. Для фиксированного $p\in S_a^2$ рассмотрим [2, с. 87, пример-задача 2]
```math
s_p(q)=2\frac{\langle p,q\rangle}{a^2}p-q.
```
Разложение $q=\lambda p+v$, $v\perp p$, даёт $s_p(q)=\lambda p-v$. Это ортогональный поворот на $\pi$ вокруг $Op$, сохраняющий сферу и метрику. При этом $s_p(p)=p$ и $(ds_p)_p=-\operatorname{id}$ на $p^\perp$. Изометрия переводит геодезическую с начальными данными $(p,v)$ в геодезическую с данными $(p,-v)$. По единственности получаем $s_p(\exp_pv)=\exp_p(-v)$ при достаточно малом $v$. Следовательно, геодезическая симметрия является изометрией.

![Геодезическая симметрия сферы](figures/04-sphere-symmetry.svg)

Рисунок 7 — Геодезическая симметрия сферы

В сечении большой окружностью (рисунок 7) дуги $pq$ и $ps_p(q)$ имеют равные длины и противоположные начальные направления в $p$.

**Итог для единичной сферы ($a=1$).**
```math
K=1,\qquad \operatorname{Scal}=2,\qquad \nabla R=0.
```
Стандартная сфера $S^2$ с обычной метрикой локально симметрична.

## Заключение

Сохранение метрики определяет изометрии и связанные с ними инварианты. В пространстве Минковского из инвариантности интервала выведены преобразования Лоренца, соотношения для собственного времени, длины и скорости. Быстрота выражает композицию сонаправленных бустов, а расчётные модели показывают количественные следствия преобразования времени.

Для сферы $S_a^2$ координатным методом получены $K=a^{-2}$, $\operatorname{Ric}=a^{-2}g$, $\operatorname{Scal}=2a^{-2}$ и проверено $\nabla R=0$. Явная изометрия $s_p$ обращает геодезические через $p$ и доказывает локальную симметричность; при $a=1$ имеем $K=1$, $\operatorname{Scal}=2$.

Общий принцип сохранения метрики реализуется в разных геометриях: пространство Минковского имеет сигнатуру $(1,3)$ и $R=0$, сфера — положительно определённую метрику и $K>0$. Буст Лоренца и геодезическая симметрия сферы являются разными изометриями; сигнатура не определяет кривизну.

## Использованная литература

1. Постников М. М. Лекции по геометрии. Семестр II. Линейная алгебра. 3-е изд., испр. СПб.: Лань, 2009. 400 с. [Каталог изданий](https://www.publ.lib.ru/ARCHIVES/P/POSTNIKOV_Mihail_Mihaylovich_%28matematik%29/_Postnikov_M.M..html) (дата обращения: 24.10.2026).
2. Постников М. М. Лекции по геометрии. Семестр V. Риманова геометрия. М.: Факториал, 1998. 496 с. [Электронный текст](https://djvu.online/file/ogoEPpiGy7UEj) (дата обращения: 24.10.2026).
3. Ландау Л. Д., Лифшиц Е. М. Теоретическая физика. Том II. Теория поля. 6-е изд., испр. и доп. М.: Наука, 1973. [Скан издания](https://radfiz.org.ua/files/k2/s3/TeopMex/Landau%2CLifshic/Landay_II.pdf) (дата обращения: 24.10.2026).

4. Particle Data Group. Review of Particle Physics // Phys. Rev. D. 2024. Vol. 110. 030001; обновление 2025 года. Таблицы: [мюон, с. 2](https://pdg.lbl.gov/2025/listings/rpp2025-list-muon.pdf); [нейтрон, с. 1, 3](https://pdg.lbl.gov/2025/listings/rpp2025-list-n.pdf) (дата обращения: 25.10.2026).
5. Bailey J., Borer K., Combley F. et al. Measurements of relativistic time dilatation for positive and negative muons in a circular orbit // Nature. 1977. Vol. 268. P. 301–305. DOI: [10.1038/268301a0](https://doi.org/10.1038/268301a0) (дата обращения: 25.10.2026).
6. Francis S., Ramsey B., Stein S. et al. Timekeeping and Time Dissemination in a Distributed Space-Based Clock Ensemble // Proc. 34th Annual PTTI Meeting. Reston, 2002. P. 201–214. [Текст NIST](https://tf.nist.gov/general/pdf/1869.pdf) (дата обращения: 25.10.2026).

7. Pierre Auger Collaboration. A Search for Point Sources of EeV Neutrons // Astrophysical Journal. 2012. Vol. 760. Art. 148. [arXiv:1211.4901](https://arxiv.org/abs/1211.4901) (дата обращения: 25.10.2026).

Для дополнительного чтения

O’Neill B. Semi-Riemannian Geometry with Applications to Relativity. Academic Press, 1983. [Страница издателя](https://shop.elsevier.com/books/semi-riemannian-geometry-with-applications-to-relativity/oneill/978-0-12-526740-3) (дата обращения: 24.10.2026).

## Приложение А. Краткая математическая часть

### А.1. Преобразование Лоренца и его следствия

Стандартная конфигурация $K\to K'$; скорость $v$ вдоль $x$, сонаправленные оси; начала совпадают при $t=t'=0$.

```math
\begin{aligned}
&X=ct,\qquad \beta=v/c,\qquad |\beta|<1,\qquad y'=y,\quad z'=z,\\
&x'=A(x-\beta X),\qquad X'=BX+Dx,\\
&(BX+Dx)^2-A^2(x-\beta X)^2=X^2-x^2,\\
&B^2-A^2\beta^2=1,\qquad D^2-A^2=-1,\qquad BD+A^2\beta=0,\\
&B\ne0,\qquad D=-\frac{A^2\beta}{B},\qquad
A^2(B^2-A^2\beta^2)=B^2\ \Longrightarrow\ A^2=B^2,\\
&A(0)=B(0)=1,\quad A,B\ne0,\quad A,B\text{ непрерывны}
\ \Longrightarrow\ A=B>0,\\
&D=-A\beta,\qquad A^2(1-\beta^2)=1,\qquad
A=B=\gamma=\frac1{\sqrt{1-\beta^2}},\\
&x'=\gamma(x-vt),\qquad
t'=\gamma\left(t-\frac{vx}{c^2}\right),\qquad y'=y,\quad z'=z,\\
&X'^2-x'^2=\gamma^2\bigl((X-\beta x)^2-(x-\beta X)^2\bigr)\\
&\hspace{25mm}=\gamma^2(1-\beta^2)(X^2-x^2)=X^2-x^2,\\
&\Delta s'^2=c^2\Delta t'^2-\Delta x'^2-\Delta y'^2-\Delta z'^2
=\Delta s^2.
\end{aligned}
```

Собственное время; световые часы; стержень, покоящийся в $K'$.

```math
\begin{aligned}
&dt>0,\quad |\mathbf u|<c,\qquad
c^2d\tau^2=ds^2=c^2dt^2\left(1-\frac{|\mathbf u|^2}{c^2}\right),\\
&d\tau=dt\sqrt{1-\frac{|\mathbf u|^2}{c^2}},\qquad
\tau_2-\tau_1=\int_{t_1}^{t_2}\sqrt{1-\frac{|\mathbf u(t)|^2}{c^2}}\,dt,\\
&\mathbf u=\mathrm{const}\ \Longrightarrow\
\Delta t=\gamma_u\Delta\tau,\qquad
\gamma_u=\left(1-|\mathbf u|^2/c^2\right)^{-1/2},\\
&\left(\frac{c\Delta t}{2}\right)^2
=h^2+\left(\frac{v\Delta t}{2}\right)^2,\qquad
\Delta\tau=\frac{2h}{c},\\
&\Delta t=\frac{2h}{c\sqrt{1-v^2/c^2}}=\gamma\Delta\tau
\qquad (\text{полный цикл часов}),\\
&\Delta t=0\ \Longrightarrow\
\Delta t'=-\gamma\frac{v\Delta x}{c^2},\\
&L_0=\Delta x'=\gamma(\Delta x-v\Delta t)
=\gamma L\ \Longrightarrow\ L=\frac{L_0}{\gamma}
\qquad (\Delta t=0).
\end{aligned}
```

Скорости: $\mathbf u=d\mathbf r/dt$, $|\mathbf u|\le c$, $|v|<c$.

```math
\begin{aligned}
&dx'=\gamma(dx-v\,dt),\qquad
dt'=\gamma\left(dt-\frac{v}{c^2}dx\right),\qquad dy'=dy,\quad dz'=dz,\\
&dt'=\gamma\left(1-\frac{vu_x}{c^2}\right)dt>0
\qquad (dt>0,\ |\mathbf u|\le c),\\
&u_x'=\frac{dx'}{dt'}=\frac{u_x-v}{1-vu_x/c^2},\\
&u_y'=\frac{u_y}{\gamma(1-vu_x/c^2)},\qquad
u_z'=\frac{u_z}{\gamma(1-vu_x/c^2)},\\
&|\mathbf u|<c:\quad c^2dt'^2-|d\mathbf r'|^2
=c^2dt^2-|d\mathbf r|^2>0\ \Longrightarrow\ |\mathbf u'|<c,\\
&\mathbf u=(v,0,0)\ \Longrightarrow\ \mathbf u'=0,\\
&\text{луч вдоль }x:\quad (u_x,u_y,u_z)=(c,0,0)
\ \Longrightarrow\ u_x'=\frac{c-v}{1-v/c}=c.
\end{aligned}
```

### А.2. Сфера: метрика, кривизна и локальная симметрия

Параметризация и метрика; $a>0$, $0<\theta<\pi$, $\alpha<\varphi<\alpha+2\pi$.

```math
\begin{aligned}
&S_a^2=\{p:\langle p,p\rangle=a^2\},\qquad
T_pS_a^2=p^\perp,\qquad g_p(v,w)=\langle v,w\rangle,\\
&r(\theta,\varphi)
=a(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta),\\
&r_\theta=a(\cos\theta\cos\varphi,\cos\theta\sin\varphi,-\sin\theta),\\
&r_\varphi=a(-\sin\theta\sin\varphi,\sin\theta\cos\varphi,0),\\
&g_{\theta\theta}=\langle r_\theta,r_\theta\rangle
=a^2(\cos^2\theta+\sin^2\theta)=a^2,\\
&g_{\theta\varphi}
=a^2\cos\theta\sin\theta(-\cos\varphi\sin\varphi+\sin\varphi\cos\varphi)=0,\\
&g_{\varphi\varphi}
=a^2\sin^2\theta(\sin^2\varphi+\cos^2\varphi)=a^2\sin^2\theta,\\
&G=\operatorname{diag}(a^2,a^2\sin^2\theta),\qquad
G^{-1}=\operatorname{diag}\left(a^{-2},(a^2\sin^2\theta)^{-1}\right),\\
&ds^2=a^2(d\theta^2+\sin^2\theta\,d\varphi^2),\qquad
\det G=a^4\sin^2\theta>0.
\end{aligned}
```

Связность и кривизна; $R(\partial_i,\partial_j)\partial_k=R^l{}_{kij}\partial_l$.

```math
\begin{aligned}
&\partial_\theta g_{\varphi\varphi}=2a^2\sin\theta\cos\theta,
\qquad \text{остальные }\partial_k g_{ij}=0,\\
&\Gamma^\theta{}_{\varphi\varphi}
=-\frac1{2a^2}\,2a^2\sin\theta\cos\theta
=-\sin\theta\cos\theta,\\
&\Gamma^\varphi{}_{\theta\varphi}
=\Gamma^\varphi{}_{\varphi\theta}
=\frac{2a^2\sin\theta\cos\theta}{2a^2\sin^2\theta}
=\cot\theta,\qquad \text{остальные }\Gamma^k{}_{ij}=0,\\
&R^\theta{}_{\varphi\theta\varphi}
=\partial_\theta\Gamma^\theta{}_{\varphi\varphi}
-\Gamma^\varphi{}_{\theta\varphi}\Gamma^\theta{}_{\varphi\varphi}\\
&\hspace{17mm}=\partial_\theta(-\sin\theta\cos\theta)
-\cot\theta(-\sin\theta\cos\theta)\\
&\hspace{17mm}=\sin^2\theta-\cos^2\theta+\cos^2\theta=\sin^2\theta,\\
&R^\varphi{}_{\theta\varphi\theta}
=-\partial_\theta\Gamma^\varphi{}_{\theta\varphi}
-(\Gamma^\varphi{}_{\theta\varphi})^2\\
&\hspace{17mm}=-\partial_\theta(\cot\theta)-\cot^2\theta
=\frac1{\sin^2\theta}-\frac{\cos^2\theta}{\sin^2\theta}=1,\\
&R^\theta{}_{\theta\theta\varphi}=0,\qquad
R^\varphi{}_{\varphi\theta\varphi}=-\partial_\varphi(\cot\theta)=0,\\
&R^l{}_{kji}=-R^l{}_{kij},\qquad R^l{}_{kii}=0.
\end{aligned}
```

Гауссова кривизна, Риччи и скалярная кривизна.

```math
\begin{aligned}
&K=\frac{g(R(\partial_\theta,\partial_\varphi)\partial_\varphi,
\partial_\theta)}{g_{\theta\theta}g_{\varphi\varphi}-g_{\theta\varphi}^2}
=\frac{a^2\sin^2\theta}{a^4\sin^2\theta}=\frac1{a^2},\\
&\operatorname{Ric}_{kj}=\sum_iR^i{}_{kij},\\
&\operatorname{Ric}_{\theta\theta}=R^\varphi{}_{\theta\varphi\theta}=1,
\qquad \operatorname{Ric}_{\varphi\varphi}
=R^\theta{}_{\varphi\theta\varphi}=\sin^2\theta,\\
&\operatorname{Ric}_{\theta\varphi}
=\operatorname{Ric}_{\varphi\theta}=0,\qquad
(\operatorname{Ric}_{ij})=\operatorname{diag}(1,\sin^2\theta)
=\frac1{a^2}G,\\
&\operatorname{Scal}=g^{\theta\theta}\operatorname{Ric}_{\theta\theta}
+g^{\varphi\varphi}\operatorname{Ric}_{\varphi\varphi}
=\frac1{a^2}+\frac{\sin^2\theta}{a^2\sin^2\theta}
=\frac2{a^2}=2K.
\end{aligned}
```

Тензорная формула; $Q$ — поворот $\mathbb R^3$.

```math
\begin{aligned}
&\frac1{a^2}(g_{\varphi\theta}\partial_\theta
-g_{\theta\theta}\partial_\varphi)
=-\partial_\varphi
=R(\partial_\theta,\partial_\varphi)\partial_\theta,\\
&\frac1{a^2}(g_{\varphi\varphi}\partial_\theta
-g_{\theta\varphi}\partial_\varphi)
=\sin^2\theta\,\partial_\theta
=R(\partial_\theta,\partial_\varphi)\partial_\varphi,\\
&R(X,Y)Z=\frac1{a^2}\bigl(g(Y,Z)X-g(X,Z)Y\bigr),\\
&R^l{}_{kij}=\frac1{a^2}
\left(\delta^l{}_i g_{jk}-\delta^l{}_j g_{ik}\right),\\
&U=(0,\pi)\times(\alpha,\alpha+2\pi),\qquad
S_a^2=\bigcup_Q Qr(U),\\
&\widetilde r=Qr,\qquad
\langle\widetilde r_i,\widetilde r_j\rangle
=\langle Qr_i,Qr_j\rangle=\langle r_i,r_j\rangle,\\
&K=a^{-2},\quad \operatorname{Ric}=a^{-2}g,\quad
\operatorname{Scal}=2a^{-2}\quad\text{на }S_a^2,\\
&R(X,Y)Z=a^{-2}(g(Y,Z)X-g(X,Z)Y)\quad\text{на }S_a^2.
\end{aligned}
```

Полная ковариантная производная; $a=\mathrm{const}$, $\nabla g=0$.

```math
\begin{aligned}
(\nabla_WR)(X,Y)Z
&=\nabla_W(R(X,Y)Z)-R(\nabla_WX,Y)Z\\
&\quad-R(X,\nabla_WY)Z-R(X,Y)(\nabla_WZ),\\
a^2(\nabla_WR)(X,Y)Z
&=W(g(Y,Z))X+g(Y,Z)\nabla_WX\\
&\quad-W(g(X,Z))Y-g(X,Z)\nabla_WY\\
&\quad-g(Y,Z)\nabla_WX+g(\nabla_WX,Z)Y\\
&\quad-g(\nabla_WY,Z)X+g(X,Z)\nabla_WY\\
&\quad-g(Y,\nabla_WZ)X+g(X,\nabla_WZ)Y\\
&=\bigl(W(g(Y,Z))-g(\nabla_WY,Z)-g(Y,\nabla_WZ)\bigr)X\\
&\quad-\bigl(W(g(X,Z))-g(\nabla_WX,Z)-g(X,\nabla_WZ)\bigr)Y\\
&=(\nabla_Wg)(Y,Z)X-(\nabla_Wg)(X,Z)Y=0,\\
\nabla R&=0\ \Longrightarrow\ S_a^2\text{ локально симметрична},\\
 a=1:&\qquad K=1,\qquad \operatorname{Scal}=2,\qquad \nabla R=0.
\end{aligned}
```
