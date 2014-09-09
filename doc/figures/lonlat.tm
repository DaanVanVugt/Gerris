<TeXmacs|1.0.6.11>

<style|article>

<\body>
  <doc-data|<doc-title|Derivation of metric coefficients for
  a<next-line>longitude-latitude coordinate
  system>|<doc-author-data|<author-name|Stéphane
  Popinet>>|<doc-date|<date>>|>

  A system of conservation laws of the form

  <\equation>
    \<partial\><rsub|t><big|int><rsub|\<Omega\>>\<b-U\>*dV+<big|int><rsub|\<partial\>\<Omega\>>\<b-F\>(\<b-U\>)\<cdot\>\<b-n\>*dS=0<label|conservation>
  </equation>

  is discretised in Gerris as

  <\equation*>
    h<rsup|2>*c<rsub|i,j>*\<partial\><rsub|t>\<b-U\>+<big|sum><rsub|faces>s<rsub|f>*h*\<b-F\><rsub|f>(\<b-U\>)=0
  </equation*>

  or

  <\equation>
    h*c<rsub|i,j>*\<partial\><rsub|t>\<b-U\>+<big|sum><rsub|faces>s<rsub|f>\<b-F\><rsub|f>(\<b-U\>)*=0<label|gerris>
  </equation>

  with <math|h> the dimensional cell size. The figure below illustrates the
  lengths of the faces of a surface element in a longitude-latitude
  coordinate system with <math|\<lambda\>> the longitude and <math|\<theta\>>
  the latitude.

  <big-figure|<postscript|lonlat.fig|0.3par|||||>|Elementary lengths for a
  surface element.>

  This leads to the discrete form for (<reference|conservation>)

  <\equation*>
    r<rsup|2>*cos \<theta\>*d\<theta\>*d\<lambda\>*\<partial\><rsub|t>\<b-U\>+<big|sum><rsub|f<rsub|\<theta\>>>r*d\<theta\>*\<b-F\><rsub|f<rsub|\<theta\>>>(\<b-U\>)+<big|sum><rsub|f<rsub|\<lambda\>>>r*cos
    \<theta\>*d\<lambda\>*\<b-F\><rsub|f<rsub|\<lambda\>>>(\<b-U\>)=0,
  </equation*>

  or

  <\equation>
    r*cos \<theta\>*d\<theta\>*d\<lambda\>*\<partial\><rsub|t>\<b-U\>+<big|sum><rsub|f<rsub|\<theta\>>>d\<theta\>*\<b-F\><rsub|f<rsub|\<theta\>>>(\<b-U\>)+<big|sum><rsub|f<rsub|\<lambda\>>>cos
    \<theta\>*d\<lambda\>*\<b-F\><rsub|f<rsub|\<lambda\>>>(\<b-U\>)=0.<label|lonlat>
  </equation>

  Equating the terms with (<reference|gerris>) gives

  <\with|color|blue>
    <\eqnarray*>
      <tformat|<table|<row|<cell|c<rsub|i,j>>|<cell|\<equiv\>>|<cell|<frac|r|h>*cos
      \<theta\> *d\<theta\>*d\<lambda\>,>>|<row|<cell|s<rsub|f<rsub|\<theta\>>>>|<cell|\<equiv\>>|<cell|d\<theta\>,>>|<row|<cell|s<rsub|f<rsub|\<lambda\>>>>|<cell|\<equiv\>>|<cell|cos
      \<theta\>* d\<lambda\>.>>>>
    </eqnarray*>
  </with>

  <subsection|Application to the Saint-Venant equations>

  For the Saint-Venant equations <math|\<b-U\>> and <math|\<b-F\>(\<b-U\>)>
  are given by

  <\equation*>
    \<b-U\>\<equiv\><matrix|<tformat|<table|<row|<cell|\<phi\>>>|<row|<cell|\<phi\>*u>>|<row|<cell|\<phi\>*v>>>>>,
    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \<b-F\><rsub|x>(\<b-U\>)\<equiv\><matrix|<tformat|<table|<row|<cell|\<phi\>*u>>|<row|<cell|\<phi\>*u<rsup|2>+g*<frac|\<phi\><rsup|2>|2>>>|<row|<cell|\<phi\>*u*v>>>>>,
    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \<b-F\><rsub|y>(\<b-U\>)\<equiv\><matrix|<tformat|<table|<row|<cell|\<phi\>*v>>|<row|<cell|\<phi\>*u*v>>|<row|<cell|\<phi\>*v<rsup|2>+g*<frac|\<phi\><rsup|2>|2>>>>>>.
  </equation*>

  Using (<reference|lonlat>), the equations in longitude-latitude coordinates
  can be written

  <\equation>
    r*cos \<theta\>*\<partial\><rsub|t>\<b-U\>+<frac|1|d\<lambda\>><big|sum><rsub|f<rsub|x>>*\<b-F\><rsub|x>(\<b-U\>)+<frac|1|d\<theta\>><big|sum><rsub|f<rsub|y>>cos
    \<theta\>*\<b-F\><rsub|y>(\<b-U\>)=0.<label|lonlat>
  </equation>

  The differential form can be recovered by taking the limit as
  <math|d\<lambda\>> and <math|d\<theta\>> tend to zero

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|r*cos
    \<theta\>*\<partial\><rsub|t>\<phi\>+\<partial\><rsub|\<lambda\>>(\<phi\>*u)+\<partial\><rsub|\<theta\>>(\<phi\>*v*cos
    \<theta\>)=0,>|<cell|>>|<row|<cell|>|<cell|r*cos
    \<theta\>*\<partial\><rsub|t>(\<phi\>*u)+\<partial\><rsub|\<lambda\>><left|(>\<phi\>*u<rsup|2>+g*<frac|\<phi\><rsup|2>|2><right|)>+\<partial\><rsub|\<theta\>><left|(>\<phi\>*u*v*cos
    \<theta\><right|)>=0,>|<cell|>>|<row|<cell|>|<cell|r*cos
    \<theta\>*\<partial\><rsub|t>(\<phi\>*v)+\<partial\><rsub|\<lambda\>>(\<phi\>*u*v)+\<partial\><rsub|\<theta\>><left|[><left|(>\<phi\>*v<rsup|2>+g*<frac|\<phi\><rsup|2>|2><right|)>*cos
    \<theta\><right|]>=0.>|<cell|>>>>
  </eqnarray*>

  This can be rewritten in advective form as

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t>\<phi\>+<frac|1|r*cos
    \<theta\>*>*<left|[>\<partial\><rsub|\<lambda\>>(\<phi\>*u)+\<partial\><rsub|\<theta\>>(\<phi\>*v*cos
    \<theta\>)<right|]>=0,>|<cell|>>|<row|<cell|>|<cell|\<partial\><rsub|t>u+<frac|1|r*cos
    \<theta\>*>*(u*\<partial\><rsub|\<lambda\>>u+v*cos
    \<theta\>*\<partial\><rsub|\<theta\>>
    u)<with|color|blue|<with|color|green|-<frac|u*v|r>*tan
    \<theta\>>>+<frac|g|r*cos \<theta\>*>*\<partial\><rsub|\<lambda\>>\<phi\>=0,>|<cell|>>|<row|<cell|>|<cell|\<partial\><rsub|t>v+<frac|1|r*cos
    \<theta\>*>*(u*\<partial\><rsub|\<lambda\>>v+v*cos \<theta\>
    \<partial\><rsub|\<theta\>> v)<with|color|green|+<frac|u<rsup|2>|r>*tan
    \<theta\>>+<frac|g|r>*\<partial\><rsub|\<theta\>>\<phi\><with|color|red|-<frac|g*\<phi\>|2*r>*tan
    \<theta\>>=0,>|<cell|>>>>
  </eqnarray*>

  where additional terms (compared to the standard form) are indicated in red
  and missing terms in green.

  <subsection|Derivation using the ``manifold approach'' of Rossmanith et
  al., 2004>

  System of ``balance laws'' on general manifolds (equation (95))

  <\equation>
    \<partial\><rsub|t>\<b-q\>+<frac|1|<sqrt|h>>*\<partial\><rsub|x<rsup|k>>\<b-f\><rsup|k>=<with|math-font-series|bold|\<psi\>><rsub|c>,<label|manifold>
  </equation>

  with <math|h> the determinant of the metric tensor and
  <with|mode|math|<with|math-font-series|bold|\<psi\>><rsub|c>> the geometric
  source terms. For the shallow-water equations

  <\equation*>
    \<b-q\>(\<b-x\>,t)\<equiv\><matrix|<tformat|<table|<row|<cell|\<phi\>>>|<row|<cell|\<phi\>*u<rsup|1>>>|<row|<cell|\<phi\>*u<rsup|2>>>>>>,<hspace|10pt>\<b-f\><rsup|k>(\<b-q\>,\<b-x\>)\<equiv\><sqrt|h>*<matrix|<tformat|<table|<row|<cell|\<phi\>*u<rsup|k>>>|<row|<cell|T<rsup|k1>>>|<row|<cell|T<rsup|k2>>>>>>,<hspace|10pt><with|math-font-series|bold|\<psi\>><rsub|c>(\<b-q\>,\<b-x\>)\<equiv\><matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|-\<Gamma\><rsub|m
    n><rsup|1>*T<rsup|n m>>>|<row|<cell|-\<Gamma\><rsub|m n><rsup|2>*T<rsup|n
    m>>>>>>,
  </equation*>

  and

  <\equation*>
    T<rsup|n m>\<equiv\>\<phi\>*u<rsup|n>*u<rsup|m>+<frac|1|2>*g*\<phi\><rsup|2>*h<rsup|n
    m>.
  </equation*>

  For a spherical coordinate system with <math|x<rsup|1>=\<lambda\>> the
  longitude and <math|x<rsup|2>=\<theta\>> the latitude, the metric is

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-h\>>|<cell|\<equiv\>>|<cell|<matrix|<tformat|<table|<row|<cell|r<rsup|2>*cos<rsup|2>
    \<theta\>>|<cell|0>>|<row|<cell|0>|<cell|r<rsup|2>>>>>>,>>|<row|<cell|<sqrt|h>>|<cell|\<equiv\>>|<cell|r<rsup|2>*cos
    \<theta\>,>>>>
  </eqnarray*>

  so that

  <\equation*>
    h<rsup|11>=<frac|1|h<rsub|11>>=<frac|1|r<rsup|2>*cos<rsup|2>
    \<theta\>>,<hspace|10pt>h<rsup|22>=<frac|1|h<rsub|22>>=<frac|1|r<rsup|2>>,<hspace|10pt>h<rsup|12>=h<rsup|21>=0.
  </equation*>

  The corresponding Christofell symbols are given by (for a general
  orthogonal metric)

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<Gamma\><rsub|m
    n><rsup|k>>|<cell|\<equiv\>>|<cell|<frac|1|2>*h<rsup|\<alpha\>k>(\<partial\><rsub|x<rsup|n>>h<rsub|\<alpha\>m>+\<partial\><rsub|x<rsup|m>>h<rsub|\<alpha\>n>-\<partial\><rsub|x<rsup|\<alpha\>>>h<rsub|m
    n>),>>|<row|<cell|\<Gamma\><rsup|1>>|<cell|\<equiv\>>|<cell|<frac|1|2*h<rsub|11>>*<matrix|<tformat|<table|<row|<cell|\<partial\><rsub|x<rsup|1>>h<rsub|11>>|<cell|\<partial\><rsub|x<rsup|2>>h<rsub|11>>>|<row|<cell|\<partial\><rsub|x<rsup|2>>h<rsub|11>>|<cell|-\<partial\><rsub|x<rsup|1>>h<rsub|2
    2>>>>>>=<matrix|<tformat|<table|<row|<cell|0>|<cell|-tan
    \<theta\>*>>|<row|<cell|-tan \<theta\>*>|<cell|0>>>>>,>>|<row|<cell|\<Gamma\><rsup|2>>|<cell|\<equiv\>>|<cell|<frac|1|2*h<rsub|22>>*<matrix|<tformat|<table|<row|<cell|-\<partial\><rsub|x<rsup|2>>h<rsub|11>>|<cell|\<partial\><rsub|x<rsup|1>>h<rsub|22>>>|<row|<cell|\<partial\><rsub|x<rsup|1>>h<rsub|22>>|<cell|\<partial\><rsub|x<rsup|2>>h<rsub|2
    2>>>>>>=<matrix|<tformat|<table|<row|<cell|*sin \<theta\>*cos
    \<theta\>>|<cell|0>>|<row|<cell|0>|<cell|0>>>>>.>>>>
  </eqnarray*>

  The geometric source term is then

  <\equation*>
    <hspace|10pt><with|math-font-series|bold|\<psi\>><rsub|c>(\<b-q\>,\<b-x\>)\<equiv\><matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|tan
    \<theta\>*(T<rsup|1 2>+*T<rsup|2 1>)>>|<row|<cell|-sin \<theta\>*cos
    \<theta\>*T<rsup|11>>>>>>=<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|2*tan
    \<theta\>*\<phi\>*u<rsup|1>*u<rsup|2>>>|<row|<cell|-sin \<theta\>*cos
    \<theta\>*\<phi\>*u<rsup|1>*u<rsup|1>-<frac|g*\<phi\><rsup|2>|2*r<rsup|2>*>*tan
    \<theta\>>>>>>
  </equation*>

  Equation (<reference|manifold>) can be developed as

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>\<phi\>+<frac|1|cos
    \<theta\>>*<left|[>\<partial\><rsub|\<lambda\>>(*cos
    \<theta\>*\<phi\>*u<rsup|1>)+\<partial\><rsub|\<theta\>>(cos
    \<theta\>*\<phi\>*u<rsup|2>)<right|]>>|<cell|=>|<cell|0,>>|<row|<cell|\<partial\><rsub|t>(\<phi\>*u<rsup|1>)+<frac|1|cos
    \<theta\>>*<left|[>\<partial\><rsub|\<lambda\>>(cos
    \<theta\>*(\<phi\>*u<rsup|1>*u<rsup|1>+<frac|g*\<phi\><rsup|2>|2*r<rsup|2>*cos<rsup|2>\<theta\>>))+\<partial\><rsub|\<theta\>>(cos
    \<theta\>*\<phi\>*u<rsup|1>*u<rsup|2>)<right|]>>|<cell|=>|<cell|-\<Gamma\><rsub|m
    n><rsup|1>*T<rsup|n m>,>>|<row|<cell|\<partial\><rsub|t>(\<phi\>*u<rsup|2>)+<frac|1|cos
    \<theta\>>*<left|[>\<partial\><rsub|\<lambda\>>(cos
    \<theta\>**\<phi\>*u<rsup|1>*u<rsup|2>)+\<partial\><rsub|\<theta\>>(cos
    \<theta\>*(\<phi\>*u<rsup|2>*u<rsup|2>+<frac|g*\<phi\><rsup|2>|2*r<rsup|2>*>))<right|]>>|<cell|=>|<cell|-\<Gamma\><rsub|m
    n><rsup|2>*T<rsup|n m>.>>>>
  </eqnarray*>

  Using

  <\equation*>
    u\<equiv\>u<rsub|1>\<equiv\>r*cos \<theta\>
    u<rsup|1>,<hspace|10pt>v\<equiv\>u<rsub|2>\<equiv\>r*u<rsup|2>,
  </equation*>

  (where exponents and indices indicate the covariant and contravariant
  vector coordinates respectively) then gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|t>\<phi\>+<frac|1|r*cos
    \<theta\>>*<left|[>\<partial\><rsub|\<lambda\>>(\<phi\>*u)+\<partial\><rsub|\<theta\>>(cos
    \<theta\>*\<phi\>*v)<right|]>>|<cell|=>|<cell|0,>>|<row|<cell|\<partial\><rsub|t>u+<frac|1|r*cos
    \<theta\>*>*(u*\<partial\><rsub|\<lambda\>>u+cos
    \<theta\>*v*\<partial\><rsub|\<theta\>>u)+<frac|g|r*cos
    \<theta\>>*\<partial\><rsub|\<lambda\>>\<phi\>-<frac|u*v|r>*tan
    \<theta\>>|<cell|=>|<cell|*0,>>|<row|<cell|\<partial\><rsub|t>v+<frac|1|r*cos
    \<theta\>>*(u*\<partial\><rsub|\<lambda\>>v+cos
    \<theta\>*v**\<partial\><rsub|\<theta\>>*v)+<frac|g|r>*\<partial\><rsub|\<theta\>>\<phi\>>|<cell|=>|<cell|-<frac|u*u|r*>*tan
    \<theta\>,>>>>
  </eqnarray*>

  which is the standard advective form.

  The geometric source term can be re-expressed in the covariant coordinate
  system as

  <\equation*>
    <hspace|10pt><with|math-font-series|bold|<wide|\<psi\>|^>>\<equiv\><matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|2*tan
    \<theta\>*\<phi\>*<frac|u<rsub|1>*u<rsub|2>|r*>>>|<row|<cell|-tan
    \<theta\>*\<phi\>*<frac|u<rsub|1>*u<rsub|1>|r*>-<frac|g*\<phi\><rsup|2>|2*r*>*tan
    \<theta\>>>>>>
  </equation*>

  <subsection|The Saint-Venant equations in general orthogonal coordinates>

  \;

  <\equation*>
    h<rsub|\<lambda\>>*h<rsub|\<theta\>>*d\<lambda\>*d\<theta\>*\<partial\><rsub|t>\<b-U\>+<big|sum><rsub|f<rsub|\<lambda\>>>h<rsub|\<theta\>>*d\<theta\>*\<b-F\><rsub|\<lambda\>>(\<b-U\>)+<big|sum><rsub|f<rsub|\<theta\>>>h<rsub|\<lambda\>>*d\<lambda\>*\<b-F\><rsub|\<theta\>>(\<b-U\>)=0
  </equation*>

  can be rewritten

  <\equation*>
    h<rsub|\<lambda\>>*h<rsub|\<theta\>>*\<partial\><rsub|t>\<b-U\>+<frac|1|d\<lambda\>>*<big|sum><rsub|f<rsub|\<lambda\>>>h<rsub|\<theta\>>*\<b-F\><rsub|\<lambda\>>(\<b-U\>)+<frac|1|d\<theta\>>*<big|sum><rsub|f<rsub|\<theta\>>>h<rsub|\<lambda\>>*\<b-F\><rsub|\<theta\>>(\<b-U\>)=0.
  </equation*>

  The differential form can be recovered by making <math|d\<lambda\>> and
  <math|d\<theta\>> tend to zero

  <\with|mode|math>
    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|h<rsub|\<lambda\>>*h<rsub|\<theta\>>*\<partial\><rsub|t>\<phi\>+\<partial\><rsub|\<lambda\>>(h<rsub|\<theta\>>*\<phi\>*u)+\<partial\><rsub|\<theta\>>(h<rsub|\<lambda\>>*\<phi\>*v)=0,>|<cell|>>|<row|<cell|>|<cell|h<rsub|\<lambda\>>*h<rsub|\<theta\>>*\<partial\><rsub|t>(\<phi\>*u)+\<partial\><rsub|\<lambda\>><left|[>h<rsub|\<theta\>>*<left|(>\<phi\>*u<rsup|2>+g*<frac|\<phi\><rsup|2>|2><right|)><right|]>+\<partial\><rsub|\<theta\>><left|(>h<rsub|\<lambda\>>*\<phi\>*u*v<right|)>=0,>|<cell|>>|<row|<cell|>|<cell|h<rsub|\<lambda\>>*h<rsub|\<theta\>>*\<partial\><rsub|t>(\<phi\>*v)+\<partial\><rsub|\<lambda\>>(h<rsub|\<theta\>>*\<phi\>*u*v)+\<partial\><rsub|\<theta\>><left|[>h<rsub|\<lambda\>>*<left|(>\<phi\>*v<rsup|2>+g*<frac|\<phi\><rsup|2>|2><right|)><right|]>=0.>|<cell|>>>>
    </eqnarray*>
  </with>

  which can be expanded as

  <\with|mode|math>
    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t>\<phi\>+<frac|u*|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>\<phi\>*+<frac|v*|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>\<phi\>+<frac|\<phi\>*|h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*<left|[>\<partial\><rsub|\<lambda\>>(h<rsub|\<theta\>>*u)+\<partial\><rsub|\<theta\>>(h<rsub|\<lambda\>>*v)<right|]>=0,>|<cell|>>|<row|<cell|>|<cell|\<partial\><rsub|t>u*+<frac|u*|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>u+<frac|v*|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>u*+<frac|g|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>\<phi\>+g*<frac|\<phi\>|2*h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*\<partial\><rsub|\<lambda\>>h<rsub|\<theta\>>=0,>|<cell|>>|<row|<cell|>|<cell|\<partial\><rsub|t>v+<frac|u|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>v+<frac|v|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>v+<frac|g|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>\<phi\>+g*<frac|\<phi\>|2*h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>h<rsub|\<lambda\>>=0.>|<cell|>>>>
    </eqnarray*>
  </with>

  Introducing the notation

  <\equation*>
    d<rsub|t>=\<partial\><rsub|t>+<frac|u*|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>+<frac|v*|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>,
  </equation*>

  this can be rewritten

  <\with|mode|math>
    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|d<rsub|t>\<phi\>+<frac|\<phi\>*|h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*<left|[>\<partial\><rsub|\<lambda\>>(h<rsub|\<theta\>>*u)+\<partial\><rsub|\<theta\>>(h<rsub|\<lambda\>>*v)<right|]>=0,>|<cell|>>|<row|<cell|>|<cell|d<rsub|t>u*<with|color|blue|<with|color|green|<with|color|green|-f<rsub|G>*v>>>+<frac|g|h<rsub|\<lambda\>>>*\<partial\><rsub|\<lambda\>>\<phi\><with|color|red|+<frac|g*\<phi\>|2*h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*\<partial\><rsub|\<lambda\>>h<rsub|\<theta\>>>=0,>|<cell|<eq-number><label|general-u>>>|<row|<cell|>|<cell|d<rsub|t>v<with|color|green|<with|color|blue|<with|color|green|+f<rsub|G>*u>>>+<frac|g|h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>\<phi\><with|color|red|+*<frac|g*\<phi\>|2*h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*\<partial\><rsub|\<theta\>>h<rsub|\<lambda\>>>=0.>|<cell|>>>>
    </eqnarray*>
  </with>

  where additional and missing terms (respective to equations 43, 44 and 45
  of Williamson et al, 1991) are indicated in red and green respectively and

  <\equation*>
    f<rsub|G>\<equiv\><frac|v*\<partial\><rsub|\<lambda\>>h<rsub|\<theta\>>-u*\<partial\><rsub|\<theta\>>h<rsub|\<lambda\>>|h<rsub|\<lambda\>>*h<rsub|\<theta\>>>*.*
  </equation*>

  <subsubsection|Application to spherical coordinates>

  For spherical coordinates

  <\equation*>
    h<rsub|\<lambda\>>=r*cos \<theta\>,<hspace|10pt>h<rsub|\<theta\>>=r,<hspace|10pt>\<partial\><rsub|\<lambda\>>h<rsub|\<theta\>>=0,<hspace|10pt>\<partial\><rsub|\<theta\>>h<rsub|\<lambda\>>=-r*sin
    \<theta\>,
  </equation*>

  which gives

  <\with|mode|math>
    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|d<rsub|t>\<phi\>+<frac|\<phi\>*|r*cos
      \<theta\>>*<left|[>\<partial\><rsub|\<lambda\>>u+\<partial\><rsub|\<theta\>>(v*cos
      \<theta\>)<right|]>=0,>|<cell|>>|<row|<cell|>|<cell|d<rsub|t>u*<with|color|blue|<with|color|green|<with|color|green|-<frac|u*v|r>*tan
      \<theta\>>>>+<frac|g|r*cos \<theta\>>*\<partial\><rsub|\<lambda\>>\<phi\>=0,>|<cell|>>|<row|<cell|>|<cell|d<rsub|t>v<with|color|green|<with|color|blue|<with|color|green|<with|color|green|+<frac|u<rsup|2>|r>*tan
      \<theta\>>>>>+<frac|g|r>*\<partial\><rsub|\<theta\>>\<phi\><with|color|red|-<frac|g*\<phi\>|2*r>*tan
      \<theta\>>=0.>|<cell|>>>>
    </eqnarray*>
  </with>

  <subsubsection|Application to polar coordinates>

  For polar coordinates

  <\equation*>
    h<rsub|\<lambda\>>=1,<hspace|10pt>h<rsub|\<theta\>>=\<lambda\>,<hspace|10pt>\<partial\><rsub|\<lambda\>>h<rsub|\<theta\>>=1,<hspace|10pt>\<partial\><rsub|\<theta\>>h<rsub|\<lambda\>>=0,
  </equation*>

  which gives

  <\with|mode|math>
    <\eqnarray*>
      <tformat|<table|<row|<cell|>|<cell|d<rsub|t>\<phi\>+<frac|\<phi\>*|\<lambda\>>*<left|[>\<partial\><rsub|\<lambda\>>(\<lambda\>*u)+\<partial\><rsub|\<theta\>>v<right|]>=0,>|<cell|>>|<row|<cell|>|<cell|d<rsub|t>u*<with|color|blue|<with|color|green|<with|color|green|-<frac|v<rsup|2>|\<lambda\>>>>>+g*\<partial\><rsub|\<lambda\>>\<phi\><with|color|red|+g*<frac|\<phi\>|2*\<lambda\>>>=0,>|<cell|>>|<row|<cell|>|<cell|d<rsub|t>v<with|color|green|<with|color|blue|<with|color|green|<with|color|green|+<frac|u*v|\<lambda\>>>>>>+<frac|g|\<lambda\>>*\<partial\><rsub|\<theta\>>\<phi\>=0.>|<cell|>>>>
    </eqnarray*>
  </with>

  <subsection|``Well-balanced'' scheme in general orthogonal coordinates>

  From Audusse et al, 2004, the one-dimensional scheme in Cartesian
  coordinates can be written

  <\equation>
    \<Delta\>x<rsub|i>*d<rsub|t>U<rsub|i>+\<cal-F\><rsub|l>-\<cal-F\><rsub|r>=S<rsub|ci>,<label|audusse>
  </equation>

  with left and right numerical fluxes

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|\<cal-F\>(U<rsub|i+1/2->,U<rsub|i+1/2+>)+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<phi\><rsup|2><rsub|i,r>-\<phi\><rsup|2><rsub|i+1/2->>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|\<cal-F\>(U<rsub|i-1/2->,U<rsub|i-1/2+>)+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<phi\><rsup|2><rsub|i,l>-\<phi\><rsup|2><rsub|i-1/2+>>>>>>,>>>>
  </eqnarray*>

  and centered source term (necessary to correct second-order imbalances)

  <\equation*>
    S<rsub|ci>=<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|(\<phi\><rsub|i,l>+\<phi\><rsub|i,r>)*(z<rsub|i,l>-z<rsub|i,r>)>>>>>.
  </equation*>

  The lake-at-rest steady state is defined by <math|u<rsub|i>=0> and
  <math|\<phi\><rsub|i,l>+z<rsub|i,l>=\<phi\><rsub|i,r>+z<rsub|i,r>=\<phi\><rsub|i>+z<rsub|i>=H>
  for all <math|i>. In this case

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<cal-F\><rsub|l>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\<phi\><rsup|2><rsub|i+1/2->>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<phi\><rsup|2><rsub|i,r>-\<phi\><rsup|2><rsub|i+1/2->>>>>>,>>|<row|<cell|\<cal-F\><rsub|r>>|<cell|=>|<cell|<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<phi\><rsup|2><rsub|i-1/2+>>>>>>+<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|*\<phi\><rsup|2><rsub|i,l>-\<phi\><rsup|2><rsub|i-1/2+>>>>>>.>>>>
  </eqnarray*>

  The centered source term becomes

  <\equation*>
    S<rsub|ci>=<frac|g|2>**<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|(\<phi\><rsub|i,l>+\<phi\><rsub|i,r>)*(H-\<phi\><rsub|i,l>-H+\<phi\><rsub|i,r>)>>>>>=<frac|g|2>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|\<phi\><rsup|2><rsub|i,r>-\<phi\><rsup|2><rsub|i,l>>>>>>,
  </equation*>

  and restores hydrostatic balance.

  For general orthogonal coordinates (<reference|audusse>) can be rewritten

  <\equation*>
    \<Delta\>x<rsub|i>*c<rsub|i>*d<rsub|t>U<rsub|i>+s<rsub|l>*\<cal-F\><rsub|l>-s<rsub|r>*\<cal-F\><rsub|r>=S<rsub|ci>+S<rsub|g>.
  </equation*>

  For the lake-at-rest steady state, this gives for the velocity component
  <math|u>

  <\equation*>
    \<Delta\>x<rsub|i>*c<rsub|i>*d<rsub|t>u+<frac|g|2>*(s<rsub|l>*\<phi\><rsup|2><rsub|i,r>-s<rsub|r>*\<phi\><rsup|2><rsub|i,l>)=S<rsub|ci>+S<rsub|g>,
  </equation*>

  with the geometric source term (red term in (<reference|general-u>))
  discretised as

  <\equation*>
    S<rsub|g>\<equiv\><frac|g|4>*(\<phi\><rsup|2><rsub|i,r>+\<phi\><rsup|2><rsub|i,l>)*(s<rsub|l>-s<rsub|r>).
  </equation*>

  The definition of <math|S<rsub|ci>> thus needs to be modified to maintain
  balance, a simple choice is

  <\equation*>
    S<rsub|ci>=<frac|g|4>*<matrix|<tformat|<table|<row|<cell|0>>|<row|<cell|(s<rsub|l>+s<rsub|r>)*(\<phi\><rsub|i,l>+\<phi\><rsub|i,r>)*(z<rsub|i,l>-z<rsub|i,r>)>>>>>,
  </equation*>

  which for the lake-at-rest steady-state gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|2|g>*\<Delta\>x<rsub|i>*c<rsub|i>*d<rsub|t>u>|<cell|=>|<cell|<frac|1|2>*(s<rsub|l>+s<rsub|r>)*(\<phi\><rsup|2><rsub|i,r>-\<phi\><rsup|2><rsub|i,l>)*+<frac|\<phi\><rsup|2><rsub|i,r>+\<phi\><rsup|2><rsub|i,l>|2>*(s<rsub|l>-s<rsub|r>)-s<rsub|l>*\<phi\><rsup|2><rsub|i,r>+s<rsub|r>*\<phi\><rsup|2><rsub|i,l>>>|<row|<cell|>|<cell|=>|<cell|0.>>>>
  </eqnarray*>

  This scheme also reduces to the Cartesian scheme for
  <math|s<rsub|l>=s<rsub|r>=c<rsub|i>=1> and to the first-order scheme for
  <math|z<rsub|i,l>=z<rsub|i,r>=z<rsub|i>> and
  <math|<rsub|>\<phi\><rsub|i,r>=\<phi\><rsub|i,l>=\<phi\><rsub|i>>.
</body>

<\references>
  <\collection>
    <associate|audusse|<tuple|7|4>>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-3|<tuple|2|2>>
    <associate|auto-4|<tuple|3|3>>
    <associate|auto-5|<tuple|3.1|4>>
    <associate|auto-6|<tuple|3.2|4>>
    <associate|auto-7|<tuple|4|4>>
    <associate|conservation|<tuple|1|1>>
    <associate|general-u|<tuple|6|4>>
    <associate|gerris|<tuple|2|1>>
    <associate|lonlat|<tuple|4|2>>
    <associate|manifold|<tuple|5|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal|Elementary lengths for a surface
      element.|<pageref|auto-1>>
    </associate>
    <\associate|toc>
      <with|par-left|<quote|1.5fn>|1<space|2spc>Application to the
      Saint-Venant equations <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|2<space|2spc>Derivation using the
      ``manifold approach'' of Rossmanith et al., 2004
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1.5fn>|3<space|2spc>The Saint-Venant equations in
      general orthogonal coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|3fn>|3.1<space|2spc>Application to spherical
      coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|3fn>|3.2<space|2spc>Application to polar
      coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1.5fn>|4<space|2spc>``Well-balanced'' scheme in
      general orthogonal coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>
    </associate>
  </collection>
</auxiliary>