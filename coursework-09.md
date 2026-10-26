# Псевдоевклидово пространство и локальная симметрия сферы

Дифференциальная геометрия и основы тензорного исчисления. Вариант №9.

## Введение

Билинейная форма задаёт скалярные произведения и квадраты интервалов, а изометрии сохраняют эту форму. В положительно определённом случае она определяет также обычные длины и углы. В псевдоевклидовом пространстве основой служит невырожденная билинейная форма; в пространстве Минковского её сохраняют преобразования Лоренца. Поэтому связь геометрии со специальной теорией относительности можно проследить через свойства формы и её изометрий.

Цель работы — раскрыть эту связь и применить метрический и тензорный аппарат к задаче варианта №9: доказать локальную симметричность обычной сферы $S^2$, вычислить её гауссову и скалярную кривизны. Сначала рассматриваются сигнатура и изотропные направления, затем выводится преобразование Лоренца и исследуются его геометрический смысл и следствия. Иллюстративные расчёты собственного времени дополняют математическое изложение.

Кривизна сферы вычисляется координатным методом: от метрики к связности, тензору кривизны и его свёрткам. Локальная симметричность доказывается построением явной изометрии. Обе части объединяет понятие преобразования, сохраняющего метрику; различие между плоским пространством Минковского и сферой выражается их тензорами кривизны.

## 1. Псевдоевклидово пространство

### 1.1. Билинейная форма и координаты

Псевдоевклидово линейное пространство — пара $(V,g)$, где $V$ — вещественное линейное пространство размерности $n$, а $g$ — симметричная невырожденная билинейная форма [1, лекция 12а]. Невырожденность означает $g(u,v)=0$ для всех $v\in V\Rightarrow u=0$; ортогональность — $g(u,v)=0$. В аффинной модели $g$ задана на пространстве переносов, и квадрат интервала равен $s^2(P,Q)=g(\overrightarrow{PQ},\overrightarrow{PQ})$.

В базисе $e=(e_1,\ldots,e_n)$, при $u=ex$, $v=ey$, получаем [1, лекции 5 и 11]

```math
G=(g(e_i,e_j)),\qquad g(u,v)=x^{\mathsf T}Gy,\qquad
G^{\mathsf T}=G,\qquad \det G\ne0.
```

Для нового базиса $e'=eC$ с обратимой матрицей $C$ координаты и матрица формы связаны равенствами

```math
x=Cx',\qquad G'=C^{\mathsf T}GC.
```

### 1.2. Канонический вид и сигнатура

По закону инерции Сильвестра в подходящем базисе [1, лекции 11–12]

```math
G=J_{p,q}=\operatorname{diag}(I_p,-I_q),\qquad p+q=n,\qquad
g(u,u)=\sum_{i=1}^{p}(x^i)^2-\sum_{i=p+1}^{n}(x^i)^2.
```

Числа $p,q$ не зависят от базиса; пару $(p,q)$ называем сигнатурой, а такой базис — псевдоортонормированным. При $q=0$ форма положительно определённая; при $p,q>0$ — неопределённая, то есть принимает на квадратах векторов значения обоих знаков.

### 1.3. Изотропные векторы и ортогональные дополнения

Ненулевой $u$ изотропен, если $g(u,u)=0$. Например, при $G=\operatorname{diag}(1,-1)$ для $u=(1,1)$ имеем $g(u,u)=0$, но $g(u,(1,0))=1$. Поэтому изотропность вектора не означает вырожденности формы; $\sqrt{|g(u,u)|}$ не является нормой. Изотропные направления этого примера — $y=\pm x$ [1, лекция 12а].

Для подпространства $W\subset V$ [1, лекция 12а, с.163]

```math
W^\perp=\{v\in V:\ g(v,w)=0\ \forall w\in W\},\qquad
\dim W^\perp=n-\dim W,
```

```math
V=W\oplus W^\perp\quad\Longleftrightarrow\quad
W\cap W^\perp=\{0\}\quad\Longleftrightarrow\quad g|_W\text{ невырожденно}.
```

Для изотропной прямой $W=\mathbb R(1,1)$ в приведённом примере $W^\perp=W$, и прямой суммы нет.

### 1.4. Изометрии и пространство Минковского

Линейная изометрия сохраняет форму; в координатах

```math
g(Lu,Lv)=g(u,v)\ \forall u,v
\quad\Longleftrightarrow\quad L^{\mathsf T}GL=G,
\qquad \det L=\pm1.
```

В псевдоортонормированном базисе изометрии образуют группу [1, лекция 12а]

```math
O(p,q)=\{L\in GL(n,\mathbb R):L^{\mathsf T}J_{p,q}L=J_{p,q}\}.
```

Здесь $GL(n,\mathbb R)$ — группа обратимых вещественных матриц. Аффинные изометрии имеют вид $F(x)=Lx+b$ и сохраняют интервалы. Частный случай — пространство Минковского сигнатуры $(1,3)$: положительное направление связано со временем, три отрицательных — с пространственными координатами. Его линейные изометрии образуют группу Лоренца $O(1,3)$; эта геометрия используется в СТО.

## 2. Геометрия пространства Минковского и СТО

В СТО пространство событий снабжается формой Минковского, а переходы между инерциальными системами сохраняют её: $\Lambda^{\mathsf T}\eta\Lambda=\eta$. Из этого условия получим преобразование Лоренца и его физические следствия.

### 2.1. События и пространство Минковского

Событие задаётся временем и местом. Введём $x^0=ct$, $x^1=x$, $x^2=y$, $x^3=z$, где $c$ — скорость света в вакууме. События — точки аффинного пространства; разности координат — векторы разделения, все компоненты которых имеют размерность длины. Форма Минковского равна [3, §§1–2, 6]

```math
\eta=\operatorname{diag}(1,-1,-1,-1),\qquad
\eta(A,B)=A^0B^0-A^1B^1-A^2B^2-A^3B^3.
```

### 2.2. Интервал и световой конус

Между двумя событиями и для бесконечно близких событий соответственно [3, §2]

```math
\Delta s^2=c^2(\Delta t)^2-|\Delta\mathbf r|^2,
\qquad ds^2=c^2dt^2-dx^2-dy^2-dz^2.
```

Интервал инвариантен. Ненулевое разделение времениподобно при $\Delta s^2>0$, светоподобно при $\Delta s^2=0$ и пространственноподобно при $\Delta s^2<0$. Нулевой интервал задаёт световой конус; времениподобные направления лежат внутри его полостей.

### 2.3. Вывод преобразования Лоренца

Пусть $K'$ движется относительно $K$ со скоростью $v$ вдоль сонаправленных осей $x,x'$, начала совпадают при $t=t'=0$, а $y'=y$, $z'=z$. Переход аффинен; совпадение начал устраняет перенос [3, §4]. Положим $X=ct$, $\beta=v/c$. Поскольку начало $K'$ удовлетворяет $x=\beta X$, $x'=0$, ищем

```math
x'=A(x-\beta X),\qquad X'=BX+Dx.
```

Из инвариантности интервала

```math
(BX+Dx)^2-A^2(x-\beta X)^2=X^2-x^2
```

сравнением коэффициентов получаем

```math
B^2-A^2\beta^2=1,\qquad D^2-A^2=-1,\qquad BD+A^2\beta=0.
```

Первое равенство даёт $B\ne0$. Подставляя $D=-A^2\beta/B$ во второе, получаем $A^2(B^2-A^2\beta^2)=B^2$, откуда $A^2=B^2$; второе равенство также даёт $A\ne0$. Выбираем непрерывную ветвь с $A(0)=B(0)=1$, сохраняющую направление времени и сонаправленность осей. На ней $A=B>0$, поэтому

```math
D=-A\beta,\qquad A^2(1-\beta^2)=1,
\qquad A=\gamma=\frac1{\sqrt{1-\beta^2}},\qquad |v|<c.
```

**Преобразование Лоренца.**

```math
x'=\gamma(x-vt),\qquad
t'=\gamma\left(t-\frac{vx}{c^2}\right),\qquad y'=y,\qquad z'=z.
```

**Проверка инвариантности.** $X'^2-x'^2=\gamma^2[(X-\beta x)^2-(x-\beta X)^2]=X^2-x^2$. С учётом $y'=y$, $z'=z$ полный интервал сохраняется. Обратное преобразование получается заменой $v$ на $-v$.

На рисунке 1 дано сечение $y=z=0$, $ct\ge0$, при $v=0{,}6c$. Оси $K'$ задаются $x'=0\Leftrightarrow x=\beta ct$ и $t'=0\Leftrightarrow ct=\beta x$; световые прямые остаются границами конуса. Гипербола $c^2t^2-x^2=\mathrm{const}>0$ объединяет события с одинаковым интервалом от начала. В условных единицах $P:(x,ct)=(3/4,5/4)\mapsto(x',ct')=(0,1)$, и $s^2=1$ в обеих системах.

![Диаграмма Минковского](figures/01-minkowski.svg)

Рисунок 1 — Диаграмма Минковского

### 2.4. Гиперболический поворот и композиция

Полученное преобразование Лоренца удобно рассмотреть как линейную изометрию псевдоевклидовой плоскости. Положим $X=ct$ и введём параметр быстроты $\psi$ равенством $\beta=v/c=\tanh\psi$. Тогда $\cosh\psi=\gamma$, $\sinh\psi=\gamma\beta$, и уже выведенные координатные формулы получают матричную запись:

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

Здесь сохранено направление перехода от $K$ к $K'$. У Ландау–Лифшица гиперболический поворот сначала записан в обратном направлении; обращение матрицы соответствует замене $\psi$ на $-\psi$ [3, §4, с. 23–24].

Смысл параметра виден по орбитам изометрий на рисунке 2. При обычном повороте единичный вектор проходит окружность $x_1^2+x_2^2=1$, а угол $\alpha$ задаёт его положение. Действие $\Lambda(\psi)$ на времениподобный вектор $e_0=(1,0)$ даёт $X=\cosh\psi$, $x=-\sinh\psi$; поэтому орбита лежит на будущей ветви $X^2-x^2=1$. При положительном $\psi$ выбранные точки идут влево вверх. Подписи $0$, $\ln2$, $\ln3$ показывают возрастание быстроты вдоль этой орбиты. Оси на обоих графиках фиксированы, масштабы координат равны; гипербола является орбитой в плоском пространстве Минковского. Параметр $\psi$ характеризует гиперболический поворот и не измеряется обычным углом между лучами рисунка.

Последнее матричное равенство следует из формул сложения гиперболических функций: при последовательных сонаправленных бустах складываются быстроты. Пусть $K_1$ движется относительно $K$ со скоростью $v_1$, а $K_2$ относительно $K_1$ — со скоростью $v_2$ вдоль той же оси. Тогда $\psi=\psi_1+\psi_2$ и $\tanh(\psi_1+\psi_2)=(\beta_1+\beta_2)/(1+\beta_1\beta_2)$. Следовательно, скорость $K_2$ относительно $K$ равна $v=(v_1+v_2)/(1+v_1v_2/c^2)$. Быстрота связывает геометрическую композицию изометрий с физическим законом сложения скоростей; условие $|v|<c$ выполняется при любой конечной сумме быстрот.

![Орбиты обычного и гиперболического поворотов](figures/07-rotation-orbits.svg)

Рисунок 2 — Орбиты обычного и гиперболического поворотов

### 2.5. Собственное время

На направленной в будущее мировой линии частицы $dt>0$, $\mathbf u=d\mathbf r/dt$, $|\mathbf u|<c$. Из $ds^2=c^2dt^2(1-|\mathbf u|^2/c^2)=c^2d\tau^2$ следует [3, §3]

```math
\begin{aligned}
d\tau&=dt\sqrt{1-\frac{|\mathbf u|^2}{c^2}},\\
\tau_2-\tau_1&=\int_{t_1}^{t_2}\sqrt{1-\frac{|\mathbf u(t)|^2}{c^2}}\,dt.
\end{aligned}
```

Собственное время $\tau$ отсчитывают часы, движущиеся с частицей; оно зависит от мировой линии. При постоянной скорости $\Delta t=\gamma_u\Delta\tau$, $\gamma_u=(1-|\mathbf u|^2/c^2)^{-1/2}$. Для светоподобной линии $ds^2=0$ собственное время не служит параметром.

Световые часы на рисунке 3 покоятся в $K'$: два зеркала разделены расстоянием $h$ поперёк движения. События $A,B,C$ — испускание, отражение и возвращение света; полный цикл длится $\Delta\tau$ в $K'$ и $\Delta t$ в $K$. В $K$ часы движутся вдоль $x$ со скоростью $v$, высота сохраняется ($y'=y$). Для половины цикла

```math
\left(\frac{c\Delta t}{2}\right)^2=h^2+\left(\frac{v\Delta t}{2}\right)^2,
\qquad \Delta\tau=\frac{2h}{c}.
```

Отсюда $\Delta t=\gamma\Delta\tau$: между теми же событиями начала и конца цикла свет в $K$ проходит больший путь со скоростью $c$. При $v=0{,}6c$ имеем $\gamma=1{,}25$; циклу $\Delta\tau=1$ нс соответствует $\Delta t=1{,}25$ нс.

![Световые часы](figures/02-light-clock.svg)

Рисунок 3 — Световые часы в двух инерциальных системах

Схема дана в пространственной плоскости $xy$; положения зеркал справа относятся к разным моментам прохождения света.

### 2.6. Одновременность и измерение длины

Преобразование времени даёт

```math
\Delta t=0\quad\Longrightarrow\quad
\Delta t'=-\gamma\frac{v\Delta x}{c^2}.
```

При $v\Delta x\ne0$ одновременные в $K$ события не одновременны в $K'$. Стержень покоится в $K'$, направлен вдоль $x'$ и имеет собственную длину $L_0$. Его длину $L$ в $K$ измеряют при одном времени $t$: $\Delta t=0$, $\Delta x=L$. Поскольку концы неподвижны в $K'$, разность их координат равна $L_0$ и для этой пары событий [3, §4]:

```math
L_0=\Delta x'=\gamma(\Delta x-v\Delta t)=\gamma L,
\qquad L=\frac{L_0}{\gamma}.
```

На рисунке 4 пара $A,B$ задаёт измерение при одном $t$, а $A,C$ — при одном $t'$ (прямая $t'=0$). При $L_0=1$ м, $v=0{,}6c$ получаем $L=0{,}8$ м. Речь идёт об одновременной регистрации концов в $K$; фотография может показывать свет, испущенный ими в разные моменты.

![Измерение длины движущегося стержня](figures/03-length.svg)

Рисунок 4 — Измерение длины движущегося стержня

### 2.7. Преобразование скоростей

Для скорости частицы $\mathbf u=(u_x,u_y,u_z)$ дифференцирование преобразования Лоренца даёт $dt'=\gamma(1-vu_x/c^2)dt$. Делим пространственные дифференциалы на $dt'$ [3, §5]:

```math
u_x'=\frac{u_x-v}{1-vu_x/c^2},\qquad
u_y'=\frac{u_y}{\gamma(1-vu_x/c^2)},\qquad
u_z'=\frac{u_z}{\gamma(1-vu_x/c^2)}.
```

При $|\mathbf u|<c$, $|v|<c$ знаменатель положителен; инвариантность интервала даёт $|\mathbf u'|<c$. Для $\mathbf u=(v,0,0)$ частица покоится в $K'$, а для луча $u_x=c$ получаем $u_x'=c$. Например, при $u_x=0{,}8c$, $v=0{,}6c$, $u_y=u_z=0$ имеем $u_x'=(0{,}8-0{,}6)c/(1-0{,}48)=5c/13\approx0{,}385c$.

### 2.8. Четырёхвекторы и метрический тензор

Греческие индексы принимают значения $0,1,2,3$, пространственные $i$ — $1,2,3$; по повторяющимся верхнему и нижнему индексам суммируем. Четырёхвектор преобразуется как $A'^\mu=\Lambda^\mu{}_{\nu}A^\nu$. Из $\Lambda^{\mathsf T}\eta\Lambda=\eta$ следует [3, §6]

```math
\eta_{\mu\nu}A'^\mu B'^\nu=A_\mu B^\mu,
\qquad A_\mu=\eta_{\mu\nu}A^\nu.
```

Поэтому $A_0=A^0$, $A_i=-A^i$. Четырёхскорость $U^\mu=dx^\mu/d\tau=\gamma_u(c,u_x,u_y,u_z)$ имеет инвариантный квадрат $\eta(U,U)=c^2$ [3, §7].

### 2.9. Математические следствия и примеры

**Коэффициент Лоренца.** При постоянном модуле скорости отношение лабораторного промежутка времени к собственному задаётся функцией $\gamma(\beta)=(1-\beta^2)^{-1/2}$, $0\le\beta<1$. Например,
```math
\gamma(0{,}6)=1{,}25,\qquad
\gamma(0{,}9)\approx2{,}29,\qquad
\gamma(0{,}99)\approx7{,}09.
```
![Коэффициент Лоренца в зависимости от скорости](figures/05-lorentz-factor.svg)

Рисунок 5 — Зависимость коэффициента Лоренца от скорости

**Вывод.** При малых $\beta$ имеем $\gamma\approx1+\beta^2/2$, а при $\beta\to1-0$ коэффициент неограниченно возрастает.

**Экспоненциальная модель собственного времени.** Если доля нераспавшихся частиц в собственном времени равна $e^{-\tau/\tau_0}$, то при постоянном модуле скорости $0<v<c$ путь $s$ соответствует $\tau=s/(\gamma v)$. Поэтому
```math
\frac{N(s)}{N_0}=\exp\left(-\frac{s}{\gamma v\tau_0}\right),
\qquad
\ln\frac{N(s)}{N_0}=-\frac{s}{\gamma v\tau_0}.
```
Здесь $N_0$ — начальное, $N(s)$ — ожидаемое число нераспавшихся частиц. В полулогарифмических координатах получаются прямые; коэффициент $\gamma$ изменяет их наклон.

Для наглядного расчёта возьмём мюоны с собственным средним временем жизни $\tau_0\approx2{,}197$ мкс [4, таблица мюона, с. 2]. Пусть они движутся вертикально вниз с высоты $h=10$ км со скоростью $v=0{,}998c$; высота и скорость выбраны для модели, потери энергии и взаимодействия с атмосферой не учитываются. Тогда
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
В системе мюона путь сокращён до $h'=h/\gamma\approx0{,}632$ км, и пролёт занимает $h'/v\approx2{,}11$ мкс собственного времени. На рисунке 6 обе зависимости рассчитаны по экспоненциальной модели; вертикальная шкала логарифмическая.

![Доля нераспавшихся мюонов на пути к поверхности](figures/06-muon-survival.svg)

Рисунок 6 — Доля нераспавшихся мюонов в двух моделях

**Вывод.** В выбранной модели учёт собственного времени увеличивает долю нераспавшихся частиц примерно на шесть порядков. Это теоретическое сопоставление, а не измерение атмосферного потока.

**Сопоставление с экспериментом.** В лабораторном опыте CERN для положительных мюонов при $\gamma=29{,}33$ измерено среднее время жизни $(64{,}419\pm0{,}058)$ мкс [5, аннотация]. Формула $\gamma\tau_0$ даёт около $64{,}44$ мкс. Здесь модуль скорости постоянен, поэтому интеграл собственного времени приводит к тому же множителю $\gamma$.

**Один закон на космическом масштабе.** Для свободного нейтрона $\tau_n\approx878{,}4$ с, то есть около 15 минут, и $m_nc^2\approx939{,}6$ МэВ [4, таблица нейтрона, с. 1, 3]. Полная энергия $E=10^{18}$ эВ задаёт $\gamma=E/(m_nc^2)$ [3, §9]. Пренебрежём гравитацией и взаимодействиями в пути. При условном расстоянии $D=8$ кпк до источника получаем
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
Здесь $\bar t_n$ и $L_n$ — средние время жизни в системе Земли и путь до распада; 1 кпк — около 3260 световых лет. Тысячи лет координатного времени совместимы с минутами собственного. Доля $P$ описывает расчётное выживание, а не регистрацию источника: в поиске Pierre Auger 2012 года значимого избытка от точечных нейтронных источников не обнаружено [7].

**Граница применимости.** Практические поправки GPS включают вклад движения по СТО около $-7{,}1$ мкс/сут относительно наземных часов и гравитационный вклад ОТО около $+45{,}7$ мкс/сут; сумма составляет $+38{,}6$ мкс/сут [6, с. 204–205, табл. 1]. Первый связан с рассмотренным собственным временем, второй требует гравитационной модели. Это различие существенно при применении плоской метрики Минковского.

## 3. Метрика, связность и кривизна

Преобразования Лоренца сохраняют постоянную форму на пространстве событий. Для поверхности метрика задаётся на каждой касательной плоскости, а для дифференцирования полей нужна связность. Введём необходимые понятия и применим их к обычной римановой сфере.

### 3.1. Метрика и связность

Регулярная параметризация $r:U\subset\mathbb R^2\to M$, взаимно однозначная на координатном участке, задаёт базис $\partial_i=r_{x^i}$ касательной плоскости. Далее $i,j,k,l,m=1,2$; по повторяющимся верхнему и нижнему индексам суммируем.

Риманова метрика $g$ — гладкая положительно определённая симметричная форма на касательных пространствах. Для поверхности в евклидовом пространстве [2, лекция 11]

```math
g_{ij}=\langle r_{x^i},r_{x^j}\rangle,\qquad
ds^2=g_{ij}\,dx^i dx^j,\qquad
(g^{ij})=G^{-1},\qquad g^{ik}g_{kj}=\delta^i{}_j.
```

Здесь $G=(g_{ij})$, $\delta^i{}_j$ — символ Кронекера. Связность Леви-Чивиты — единственная связность без кручения, согласованная с метрикой [2, лекция 11]:

```math
\nabla_XY-\nabla_YX=[X,Y],\qquad \nabla g=0,\qquad
\nabla_{\partial_i}\partial_j=\Gamma^k{}_{ij}\partial_k.
```

Ковариантная производная линейна по $X$ над гладкими функциями, по $Y$ над $\mathbb R$ и удовлетворяет правилу $\nabla_X(fY)=X(f)Y+f\nabla_XY$; $[X,Y]$ — коммутатор полей. Для координатных полей $[\partial_i,\partial_j]=0$, поэтому $\Gamma^k{}_{ij}=\Gamma^k{}_{ji}$, причём [2, с.174–175]

```math
\Gamma^k{}_{ij}=\frac12g^{k\ell}
\left(\partial_i g_{j\ell}+\partial_jg_{i\ell}-\partial_\ell g_{ij}\right).
```

### 3.2. Кривизна и её свёртки

Фиксируем знак кривизны и порядок индексов [2, лекция 2]:

```math
R(X,Y)Z=\nabla_X\nabla_YZ-\nabla_Y\nabla_XZ-\nabla_{[X,Y]}Z,
\qquad R(\partial_i,\partial_j)\partial_k=R^l{}_{kij}\partial_l,
```

```math
R^l{}_{kij}=\partial_i\Gamma^l{}_{jk}-\partial_j\Gamma^l{}_{ik}
+\Gamma^m{}_{jk}\Gamma^l{}_{im}-\Gamma^m{}_{ik}\Gamma^l{}_{jm}.
```

Тензор $R$ имеет тип $(1,3)$ и кососимметричен по $i,j$. Свёртки Риччи и скалярной кривизны [2, лекции 15–17], а также гауссова кривизна для независимых $X,Y\in T_pM$ равны

```math
\operatorname{Ric}_{kj}=\sum_iR^i{}_{kij},\qquad
\operatorname{Scal}=\sum_{k,j}g^{kj}\operatorname{Ric}_{kj},\qquad
K=\frac{g(R(X,Y)Y,X)}{g(X,X)g(Y,Y)-g(X,Y)^2}.
```

На поверхности плоскость $T_pM$ единственна, и $K$ не зависит от выбора её базиса. В инерциальных координатах метрика Минковского постоянна, поэтому $\Gamma=0$ и $R=0$: неопределённая сигнатура сама по себе не означает кривизну.

### 3.3. Локальная геодезическая симметрия

Геодезическая определяется уравнением $\nabla_{\dot\gamma}\dot\gamma=0$. Экспонента задаётся начальными условиями

```math
\exp_p(v)=\gamma_v(1),\qquad \gamma_v(0)=p,\qquad
\dot\gamma_v(0)=v,\qquad (d\exp_p)_0=\operatorname{id}.
```

На достаточно малом шаре $B\subset T_pM$ она является диффеоморфизмом на нормальную окрестность [2, лекция 1, с.19–20]. Геодезическая симметрия имеет вид

```math
s_p(\exp_pv)=\exp_p(-v),\qquad v\in B,\qquad
s_p(p)=p,\qquad (ds_p)_p=-\operatorname{id}.
```

Риманово пространство локально симметрично, если все $s_p$ — локальные изометрии. Для сферы такую изометрию построим непосредственно.

## 4. Геометрия и локальная симметрия сферы

### 4.1. Сфера, параметризация и метрика

Рассмотрим сферу постоянного радиуса $a>0$ с метрикой, индуцированной евклидовым скалярным произведением. Для сферы из задания затем положим $a=1$.

Сфера $S_a^2=\{x\in\mathbb R^3:\langle x,x\rangle=a^2\}$ — регулярная поверхность уровня $F(x)=\langle x,x\rangle-a^2$: при $p\in S_a^2$ имеем $\operatorname{grad}F(p)=2p\ne0$. Поэтому

```math
T_pS_a^2=\ker dF_p=\{v:\langle p,v\rangle=0\}=p^\perp,
\qquad g_p(v,w)=\langle v,w\rangle.
```

Для $0\ne v\in T_pS_a^2$ имеем $g_p(v,v)=|v|^2>0$, поэтому метрика риманова на всей сфере.

Введём $x^1=\theta$, $x^2=\varphi$:

```math
r(\theta,\varphi)=a(\sin\theta\cos\varphi,\sin\theta\sin\varphi,\cos\theta),
\qquad 0<\theta<\pi,\quad \alpha<\varphi<\alpha+2\pi.
```

Карта исключает полюса и меридианный разрез $\varphi=\alpha$; на указанной области она взаимно однозначна. Её производные:

```math
\begin{aligned}
r_\theta&=a(\cos\theta\cos\varphi,\cos\theta\sin\varphi,-\sin\theta),\\
r_\varphi&=a(-\sin\theta\sin\varphi,\sin\theta\cos\varphi,0).
\end{aligned}
```

Вычислим их скалярные произведения:

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

При $a>0$ и $0<\theta<\pi$ векторы $r_\theta,r_\varphi$ независимы; деления на $a$ и $\sin\theta$ допустимы.

### 4.2. Символы Кристоффеля

Единственная не тождественно нулевая производная метрики:

```math
\partial_\theta g_{\varphi\varphi}=2a^2\sin\theta\cos\theta
```

Подстановка в формулу Кристоффеля даёт

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

По симметрии нижних индексов $\Gamma^\varphi{}_{\varphi\theta}=\cot\theta$; остальные символы равны нулю.

### 4.3. Компоненты тензора кривизны

При соглашении $R(\partial_i,\partial_j)\partial_k=R^l{}_{kij}\partial_l$ получаем

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

При $i=j$ компоненты равны нулю; перестановка $i,j$ меняет знак. Поэтому все компоненты определены и

```math
R(\partial_\theta,\partial_\varphi)\partial_\theta=-\partial_\varphi,
\qquad R(\partial_\theta,\partial_\varphi)\partial_\varphi
=\sin^2\theta\,\partial_\theta.
```

### 4.4. Кривизна и её свёртки

Для $X=\partial_\theta$, $Y=\partial_\varphi$

```math
K=\frac{g(\sin^2\theta\,\partial_\theta,\partial_\theta)}
{g_{\theta\theta}g_{\varphi\varphi}-g_{\theta\varphi}^2}
=\frac{a^2\sin^2\theta}{a^4\sin^2\theta}
=\frac1{a^2}.
```

Свёртка $\operatorname{Ric}_{kj}=\sum_iR^i{}_{kij}$ даёт

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

Сворачивая с обратной метрикой, получаем

```math
\operatorname{Scal}
=g^{\theta\theta}\operatorname{Ric}_{\theta\theta}
+g^{\varphi\varphi}\operatorname{Ric}_{\varphi\varphi}
=\frac1{a^2}+\frac{\sin^2\theta}{a^2\sin^2\theta}
=\frac2{a^2}=2K.
```

### 4.5. Переход ко всей сфере

При $\theta=0,\pi$ вырождается карта, а не метрика на $p^\perp$. Для каждого $p$ выберем полярную ось и разрез, не содержащие $p$. Поворотная карта $\widetilde r=Qr$ имеет ту же метрику, поскольку $Q$ сохраняет скалярные произведения. Такие карты покрывают сферу; по тензорности результаты справедливы на всей $S_a^2$.

### 4.6. Тензорная формула кривизны сферы

Проверим по вычисленным компонентам равенство

```math
R(X,Y)Z=\frac1{a^2}\bigl(g(Y,Z)X-g(X,Z)Y\bigr).
```

При $X=\partial_\theta$, $Y=\partial_\varphi$ правая часть для двух базисных $Z$ равна

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

Это совпадает с §4.3. При перестановке $X,Y$ обе части меняют знак, при равных базисных $X,Y$ равны нулю; по тензорности равенство верно для любых $X,Y,Z$. Поворотные карты распространяют его на всю сферу. В координатах

```math
R^l{}_{kij}=\frac1{a^2}\left(\delta^l{}_i g_{jk}-\delta^l{}_j g_{ik}\right).
```

### 4.7. Ковариантное постоянство кривизны

Радиус $a$ постоянен, $\nabla g=0$. Для полного тензора

```math
\begin{aligned}
(\nabla_WR)(X,Y)Z
&=\nabla_W(R(X,Y)Z)-R(\nabla_WX,Y)Z\\
&\quad-R(X,\nabla_WY)Z-R(X,Y)(\nabla_WZ).
\end{aligned}
```

Подставляем тензорную формулу и применяем правило Лейбница; слагаемые с производными $X,Y,Z$ сокращаются:

```math
(\nabla_WR)(X,Y)Z=\frac1{a^2}\bigl((\nabla_Wg)(Y,Z)X-(\nabla_Wg)(X,Z)Y\bigr)=0.
```

Таким образом, $\nabla R=0$ на всей сфере. Это проверка ковариантного постоянства вычисленного тензора: отдельные координатные компоненты $R$ могут быть непостоянны.

### 4.8. Локальная симметрия сферы

Ортогональные преобразования $\mathbb R^3$, сохраняющие сферу, сохраняют и её индуцированную метрику [2, с. 381]. Для фиксированного $p\in S_a^2$ возьмём стандартную симметрию [2, с. 87, пример-задача 2]:

```math
s_p(q)=2\frac{\langle p,q\rangle}{a^2}p-q.
```

Разложение $q=\lambda p+v$, $v\perp p$, даёт $s_p(q)=\lambda p-v$. Поэтому это ортогональный поворот на $\pi$ вокруг $Op$, сохраняющий сферу и её метрику. При этом $s_p(p)=p$ и $(ds_p)_p=-\operatorname{id}$ на $p^\perp$. Изометрия сохраняет геодезические, поэтому по единственности геодезической с заданными начальными данными $s_p(\exp_pv)=\exp_p(-v)$ при малом $v$. Следовательно, геодезическая симметрия является изометрией.

![Геодезическая симметрия сферы](figures/04-sphere-symmetry.svg)

Рисунок 7 — Геодезическая симметрия сферы

В сечении большой окружностью (рисунок 7) дуги $pq$ и $ps_p(q)$ имеют равные длины и противоположные начальные направления в $p$.

**Итог для единичной сферы ($a=1$).**

```math
K=1,\qquad \operatorname{Scal}=2,\qquad \nabla R=0.
```

Стандартная сфера $S^2$ с обычной метрикой локально симметрична.

## Заключение

Псевдоевклидова форма приводит к пространству Минковского и группе Лоренца. Из сохранения интервала получено преобразование между инерциальными системами; его запись через гиперболический угол объясняет композицию сонаправленных бустов. Собственное время, сокращение длины и преобразование скоростей рассмотрены как математические следствия этой геометрии. Графики коэффициента Лоренца и экспоненциальной модели показывают масштаб соответствующих изменений.

Для сферы радиуса $a$ координатным методом вычислены связность, тензор кривизны и его свёртки: $K=a^{-2}$, $\operatorname{Ric}=a^{-2}g$, $\operatorname{Scal}=2a^{-2}$. Проверено $\nabla R=0$, а явная изометрия $s_p$ обращает геодезические через $p$ и доказывает локальную симметричность. Для единичной сферы ответ равен $K=1$, $\operatorname{Scal}=2$.

Таким образом, метрика и её изометрии дают общий язык обеих частей работы. Пространство Минковского с постоянной метрикой имеет нулевой тензор кривизны, тогда как обычная риманова сфера имеет положительную постоянную кривизну.

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
