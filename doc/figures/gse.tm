<TeXmacs|1.0.6.11>

<style|<tuple|article|maxima>>

<\body>
  <section|GSE alleviation using spatial filtering (or is it diffusion?)>

  Following Tolman, 2002, Appendix A, the filtered value is defined as

  <\equation>
    F<rsub|avg>(\<b-x\>)=<frac|1|3>*F(\<b-x\>)+<frac|1|6>*<big|sum><rsup|4><rsub|n=1>F(\<b-x\>+\<b-r\><rsub|n>),<label|avg>
  </equation>

  with

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<b-r\><rsub|1>>|<cell|=>|<cell|\<b-s\>+\<b-n\>,>>|<row|<cell|\<b-r\><rsub|2>>|<cell|=>|<cell|-\<b-s\>+\<b-n\>,>>|<row|<cell|\<b-r\><rsub|3>>|<cell|=>|<cell|-\<b-s\>-\<b-n\>,>>|<row|<cell|\<b-r\><rsub|4>>|<cell|=>|<cell|\<b-s\>-\<b-n\>>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|math-font-series|bold|s>>|<cell|=>|<cell|\<alpha\>*<matrix|<tformat|<table|<row|<cell|cos
    \<theta\>>>|<row|<cell|sin \<theta\>>>>>>,>>|<row|<cell|\<b-n\>>|<cell|=>|<cell|\<beta\>*<matrix|<tformat|<table|<row|<cell|-sin
    \<theta\>>>|<row|<cell|cos \<theta\>>>>>>.>>>>
  </eqnarray*>

  To third-order the Taylor expansion of <math|F> can be written

  <\equation>
    F(\<b-x\>+\<b-r\>)=F(\<b-x\>)+\<b-r\>*\<b-F\><rprime|'>(\<b-x\>)+<frac|1|2>*\<b-r\><rsup|T>*\<b-F\><rprime|''>(\<b-x\>)*\<b-r\>,<label|taylor>
  </equation>

  with the gradient <math|\<b-F\><rprime|'>> and Hessian
  <math|\<b-F\><rprime|''>> given by

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|mode|text|<math|\<b-F\><rprime|'>>>>|<cell|\<equiv\>>|<cell|<matrix|<tformat|<table|<row|<cell|<frac|\<partial\>F|\<partial\>x>>>|<row|<cell|<frac|\<partial\>F|\<partial\>y>>>>>>,>>|<row|<cell|<with|mode|text|<math|\<b-F\><rprime|''>>>>|<cell|\<equiv\>>|<cell|<matrix|<tformat|<table|<row|<cell|<frac|\<partial\><rsup|2>F|\<partial\>x<rsup|2>>>|<cell|<frac|\<partial\><rsup|2>F|\<partial\>x\<partial\>y>>>|<row|<cell|<frac|\<partial\><rsup|2>F|\<partial\>x\<partial\>y>>|<cell|<frac|\<partial\><rsup|2>F|\<partial\>y<rsup|2>>>>>>>.>>>>
  </eqnarray*>

  Using (<reference|taylor>) in (<reference|avg>) gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|F<rsub|avg>(\<b-x\>)>|<cell|=>|<cell|<frac|1|3>*F(\<b-x\>)+<frac|1|6>*<big|sum><rsup|4><rsub|n=1>F(\<b-x\>)+\<b-r\><rsub|n>*\<b-F\><rprime|'>(\<b-x\>)+<frac|1|2>*\<b-r\><rsub|n><rsup|T>*\<b-F\><rprime|''>(\<b-x\>)*\<b-r\><rsub|n>>>|<row|<cell|>|<cell|=>|<cell|F(\<b-x\>)+<frac|1|6>*\<b-F\><rprime|'>(\<b-x\>)*<big|sum><rsup|4><rsub|n=1>\<b-r\><rsub|n>+<frac|1|12>*<big|sum><rsup|4><rsub|n=1>\<b-r\><rsub|n><rsup|T>*\<b-F\><rprime|''>(\<b-x\>)*\<b-r\><rsub|n>>>|<row|<cell|>|<cell|=>|<cell|F(\<b-x\>)+<frac|1|12>*<big|sum><rsup|4><rsub|n=1>\<b-r\><rsub|n><rsup|T>*\<b-F\><rprime|''>(\<b-x\>)*\<b-r\><rsub|n>>>>>
  </eqnarray*>

  I am lazy so I use Maxima to obtain the simplified form

  <with|prog-language|maxima|prog-session|default|<\session>
    <\input|<with|color|red|(<with|math-font-family|rm|%i>1)
    <with|color|black|>>>
      s:matrix([alpha*cos(theta)],[alpha*sin(theta)]);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o1>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|\<alpha\>*cos
      <left|(>\<vartheta\><right|)>>>|<row|<cell|\<alpha\>*sin
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>2)
    <with|color|black|>>>
      n:matrix([-beta*sin(theta)],[beta*cos(theta)]);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o2>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|-\<beta\>*sin
      <left|(>\<vartheta\><right|)>>>|<row|<cell|\<beta\>*cos
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>3)
    <with|color|black|>>>
      F:matrix( [Fxx,Fxy], [Fxy,Fyy] );
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o3>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|<with|math-font-family|rm|Fxx>>|<cell|<with|math-font-family|rm|Fxy>>>|<row|<cell|<with|math-font-family|rm|Fxy>>|<cell|<with|math-font-family|rm|Fyy>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>4)
    <with|color|black|>>>
      r1: s+n;
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o4>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|\<alpha\>*cos
      <left|(>\<vartheta\><right|)>-\<beta\>*sin
      <left|(>\<vartheta\><right|)>>>|<row|<cell|\<alpha\>*sin
      <left|(>\<vartheta\><right|)>+\<beta\>*cos
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>5)
    <with|color|black|>>>
      r2: -s+n;
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o5>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|-\<beta\>*sin
      <left|(>\<vartheta\><right|)>-\<alpha\>*cos
      <left|(>\<vartheta\><right|)>>>|<row|<cell|\<beta\>*cos
      <left|(>\<vartheta\><right|)>-\<alpha\>*sin
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>6)
    <with|color|black|>>>
      r3: -s-n;
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o6>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|\<beta\>*sin
      <left|(>\<vartheta\><right|)>-\<alpha\>*cos
      <left|(>\<vartheta\><right|)>>>|<row|<cell|-\<alpha\>*sin
      <left|(>\<vartheta\><right|)>-\<beta\>*cos
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>7)
    <with|color|black|>>>
      r4: s-n;
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o10>)
      <with|color|black|>><left|(><tabular*|<tformat|<table|<row|<cell|\<beta\>*sin
      <left|(>\<vartheta\><right|)>+\<alpha\>*cos
      <left|(>\<vartheta\><right|)>>>|<row|<cell|\<alpha\>*sin
      <left|(>\<vartheta\><right|)>-\<beta\>*cos
      <left|(>\<vartheta\><right|)>>>>>><right|)>>
    </output>

    <\input|<with|color|red|(<with|math-font-family|rm|%i>11)
    <with|color|black|>>>
      ratsimp(transpose(r1).F.r1+transpose(r2).F.r2+transpose(r3).F.r3+transpose(r4).F.r4);
    </input>

    <\output>
      <with|mode|math|math-display|true|<with|mode|text|font-family|tt|color|red|(<with|math-font-family|rm|%o12>)
      <with|color|black|>><left|(>4*\<alpha\><rsup|2>*<with|math-font-family|rm|Fyy>+4*\<beta\><rsup|2>*<with|math-font-family|rm|Fxx><right|)>*sin
      <left|(>\<vartheta\><right|)><rsup|2>+<left|(>8*\<alpha\><rsup|2>-8*\<beta\><rsup|2><right|)>*<with|math-font-family|rm|Fxy>*cos
      <left|(>\<vartheta\><right|)>*sin <left|(>\<vartheta\><right|)>+<left|(>4*\<beta\><rsup|2>*<with|math-font-family|rm|Fyy>+4*\<alpha\><rsup|2>*<with|math-font-family|rm|Fxx><right|)>*cos
      <left|(>\<vartheta\><right|)><rsup|2>>
    </output>
  </session>>

  We get

  <\eqnarray*>
    <tformat|<table|<row|<cell|F<rsub|avg>(\<b-x\>)>|<cell|=>|<cell|F(\<b-x\>)+>>|<row|<cell|>|<cell|>|<cell|<frac|1|3>**(\<alpha\><rsup|2>*cos<rsup|2>
    \<theta\>+\<beta\><rsup|2>*sin<rsup|2> \<theta\>)*F<rsub|x
    x>*+>>|<row|<cell|>|<cell|>|<cell|<frac|1|3>**(\<alpha\><rsup|2>*sin<rsup|2>
    \<theta\>+\<beta\><rsup|2>*cos<rsup|2> \<theta\>)*F<rsub|y
    y>+>>|<row|<cell|>|<cell|>|<cell|<frac|2|3>
    (\<alpha\><rsup|2>-\<beta\><rsup|2>)*cos \<theta\>*sin \<theta\>*F<rsub|x
    y>,>>>>
  </eqnarray*>

  which can be rewritten

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|F<rsub|avg>(\<b-x\>)-F(\<b-x\>)|\<Delta\>t>>|<cell|=>|<cell|D<rsub|x
    x>*F<rsub|x x>+D<rsub|y y>*F<rsub|y y>+2*D<rsub|x y>*F<rsub|x y>,>>>>
  </eqnarray*>

  with

  <\eqnarray*>
    <tformat|<table|<row|<cell|D<rsub|x*x>>|<cell|\<equiv\>>|<cell|D<rsub|s*s>*cos<rsup|2>
    \<theta\>+D<rsub|n*n>*sin<rsup|2> \<theta\>,>>|<row|<cell|D<rsub|y*y>>|<cell|\<equiv\>>|<cell|D<rsub|s*s>*sin<rsup|2>
    \<theta\>+D<rsub|n*n>*cos<rsup|2> \<theta\>,>>|<row|<cell|D<rsub|x*y>>|<cell|\<equiv\>>|<cell|(D<rsub|s*s>-D<rsub|n*n>)*cos
    \<theta\>*sin \<theta\>,>>>>
  </eqnarray*>

  and

  <\eqnarray*>
    <tformat|<table|<row|<cell|D<rsub|s*s>>|<cell|\<equiv\>>|<cell|<frac|\<alpha\><rsup|2>|3*\<Delta\>t>,>>|<row|<cell|D<rsub|n*n>>|<cell|\<equiv\>>|<cell|<frac|\<beta\><rsup|2>|3*\<Delta\>t>.>>>>
  </eqnarray*>

  This means that to third-order accuracy the spatial-filtering scheme of
  Tolman is formally equivalent to a first-order time discretisation of the
  diffusion equation of Booij and Holthuijsen, 1987. Furthermore Tolman takes
  <math|\<alpha\>> and <math|\<beta\>> as

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<alpha\>>|<cell|\<equiv\>>|<cell|\<alpha\><rsub|s>*\<Delta\>c<rsub|g>*\<Delta\>t,>>|<row|<cell|\<beta\>>|<cell|\<equiv\>>|<cell|\<alpha\><rsub|n>*c<rsub|g>*\<Delta\>\<theta\>*\<Delta\>t,>>>>
  </eqnarray*>

  with <math|\<Delta\>c<rsub|g>\<equiv\>><with|mode|math|(\<gamma\>-\<gamma\><rsup|-1>)*c<rsub|g>/2>,
  whereas Booij and Holthuijsen take

  <\eqnarray*>
    <tformat|<table|<row|<cell|<wide|D|^><rsub|s*s>>|<cell|\<equiv\>>|<cell|<frac|1|12>*(\<Delta\>c<rsub|g>)<rsup|2>*T<rsub|s>,>>|<row|<cell|<wide|D|^><rsub|n*n>>|<cell|\<equiv\>>|<cell|<frac|1|12>*(c<rsub|g>*\<Delta\>\<theta\>)<rsup|2>*T<rsub|s>.>>>>
  </eqnarray*>

  The two schemes are actually equivalent when <math|D<rsub|s*s>>,
  <math|D<rsub|n*n>> are identical to <math|<wide|D|^><rsub|s*s>>,
  <math|<wide|D|^><rsub|n*n>> which gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|T<rsub|s>>|<cell|=>|<cell|4*\<alpha\><rsup|2><rsub|s>*\<Delta\>t,>>|<row|<cell|\<alpha\><rsub|n>>|<cell|=>|<cell|\<alpha\><rsub|s>.>>>>
  </eqnarray*>

  This is obviously a much smaller value for <math|T<rsub|s>> than the 4 days
  used in Tolman 2002 and I am not sure I understand where this discrepancy
  comes from... The stability criterion for the explicit diffusion scheme
  then becomes

  <\equation*>
    \<Delta\>t\<less\><frac|(\<Delta\>x)<rsup|2>|(c<rsub|g>*\<Delta\>\<theta\>)<rsup|2>*4*\<alpha\><rsup|2><rsub|s>*\<Delta\>t>
  </equation*>

  which simplifies as

  <\equation*>
    \<Delta\>t\<less\><frac|\<Delta\>x|2*\<alpha\><rsub|s>*c<rsub|g>*\<Delta\>\<theta\>>
  </equation*>

  which is the standard CFL stability criterion with a CFL number of
  <math|(2*\<alpha\><rsub|s>*\<Delta\>\<theta\>)<rsup|-1>>. For 24 wave
  directions this number is larger than one for
  <math|\<alpha\><rsub|s>\<lesssim\>2>, which means that this scheme does not
  restrict the timestep if wave advection is resolved using an explicit
  scheme (as is usually the case).
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1>>
    <associate|avg|<tuple|1|1>>
    <associate|taylor|<tuple|2|1>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>GSE
      alleviation using spatial filtering (or is it diffusion?)>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>