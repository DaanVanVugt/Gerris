<TeXmacs|1.0.7.3>

<style|generic>

<\body>
  <section|Vertical diffusion>

  We consider the one-dimensional diffusion equation

  <\equation*>
    \<partial\><rsub|t>u=\<partial\><rsub|z>\<mu\>*\<partial\><rsub|z>u
  </equation*>

  with boundary conditions

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<partial\><rsub|z>u<mid|\|><rsub|z<rsub|t>>=<wide|u|\<dot\>><rsub|t>>|<cell|>|<cell|>>|<row|<cell|u(z<rsub|b>)=u<rsub|b>+\<lambda\><rsub|b>*\<partial\><rsub|z>u<mid|\|><rsub|z<rsub|b>>>|<cell|>|<cell|>>>>
  </eqnarray*>

  This is discretised implicitly as

  <\equation*>
    u<rsup|n+1>-u<rsup|n>=\<Delta\>t*\<partial\><rsub|z>\<mu\>*\<partial\><rsub|z>u<rsup|n+1>
  </equation*>

  Using the finite volume approximation for each level <math|l>

  <\equation*>
    \<partial\><rsub|z>\<mu\>*\<partial\><rsub|z>u\<simeq\><frac|2|\<delta\><rsub|l>>*<left|[>\<mu\><rsub|l+1/2>*<frac|u<rsub|l+1>-u<rsub|l>|\<delta\><rsub|l>+\<delta\><rsub|l+1>>-\<mu\><rsub|l-1/2>*<frac|u<rsub|l>-u<rsub|l-1>|\<delta\><rsub|l>+\<delta\><rsub|l-1>><right|]>
  </equation*>

  we get

  <\equation*>
    u<rsub|l><rsup|n+1>*<left|[>1+a<rsub|l+1/2>+a<rsub|l-1/2><right|]>-u<rsup|n+1><rsub|l+1>*a<rsub|l+1/2>-u<rsup|n+1><rsub|l-1>*a<rsub|l-1/2>=u<rsub|l><rsup|n>
  </equation*>

  with

  <\equation*>
    a<rsub|l+1/2>\<equiv\><frac|2*\<Delta\>t|\<delta\><rsub|l>>*<frac|\<mu\><rsub|l+1/2>|\<delta\><rsub|l>+\<delta\><rsub|l+1>>
  </equation*>

  for <math|0\<less\>l\<less\>N-1> where <math|N> is the number of levels.
  This tridiagonal system is closed by the boundary conditions which can be
  discretised as

  <\eqnarray*>
    <tformat|<table|<row|<cell|<frac|u<rsub|N>-u<rsub|N-1>|\<delta\><rsub|N-1>>\<simeq\><wide|u|\<dot\>><rsub|t>>|<cell|>|<cell|>>|<row|<cell|<frac|u<rsub|-1>+u<rsub|0>|2>\<simeq\>u<rsub|b>+\<lambda\><rsub|b>*<frac|u<rsub|0>-u<rsub|-1>|\<delta\><rsub|0>>>|<cell|>|<cell|>>>>
  </eqnarray*>

  from which we get the (ghost) boundary values <math|u<rsub|N>> and
  <math|u<rsub|-1>>

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<rsub|N>*>|<cell|=>|<cell|u<rsub|N-1>+<wide|u|\<dot\>><rsub|t>*\<delta\><rsub|N-1>>>|<row|<cell|u<rsub|-1>*>|<cell|=>|<cell|<frac|2*\<delta\><rsub|0>|2*\<lambda\><rsub|b>+\<delta\><rsub|0>>*u<rsub|b>+<left|(><frac|2*\<lambda\><rsub|b>-\<delta\><rsub|0>|2*\<lambda\><rsub|b>+\<delta\><rsub|0>><right|)>*u<rsub|0><eq-number><label|ub>>>>>
  </eqnarray*>

  This gives the additional equations

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|u<rsub|N-1><rsup|n+1>*<left|[>1+a<rsub|N-3/2><right|]>-u<rsup|n+1><rsub|N-2>*a<rsub|N-3/2>=u<rsub|N-1><rsup|n>+<wide|u|\<dot\>><rsub|t>*\<delta\><rsub|N-1>*a<rsub|N-1/2>>|<cell|>>|<row|<cell|>|<cell|u<rsub|0><rsup|n+1>*<left|[>1+a<rsub|1/2>+a<rsub|-1/2>-<left|(><frac|2*\<lambda\><rsub|b>-\<delta\><rsub|0>|2*\<lambda\><rsub|b>+\<delta\><rsub|0>><right|)>*a<rsub|-1/2><right|]>-u<rsup|n+1><rsub|1>*a<rsub|1/2>=u<rsub|0><rsup|n>+<frac|2*\<delta\><rsub|0>|2*\<lambda\><rsub|b>+\<delta\><rsub|0>>*u<rsub|b>*a<rsub|-1/2>>|<cell|>>>>
  </eqnarray*>

  <subsection|Bottom friction>

  According to Audusse, Bristeau and Decoene (2008) equation (20), the
  friction terms for the bottom layer (<math|><math|l=0>) can be discretised
  as

  <\with|mode|math>
    <\equation*>
      \<partial\><rsub|t>u<rsub|0>=-<frac|k|\<delta\><rsub|0>>*u<rsub|0>+<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>><right|]>,
    </equation*>
  </with>

  This gives the additional equation

  <\equation*>
    u<rsup|n+1><rsub|0>*(1+a<rsub|1/2>+<frac|k|\<delta\><rsub|0>>*\<Delta\>t)-u<rsup|n+1><rsub|1>*a<rsub|1/2>=u<rsup|n><rsub|0>
  </equation*>

  <subsection|Link between Navier condition and bottom friction>

  According to Audusse, Bristeau and Decoene (2008) equation (20), the
  friction terms for the bottom layer (<math|><math|l=0>) can be discretised
  as

  <\equation*>
    \<partial\><rsub|t>u<rsub|0>=-<frac|k|\<delta\><rsub|0>>*u<rsub|0>+<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>><right|]>,
  </equation*>

  with <math|k> the bottom friction coefficient. The boundary conditions for
  our system must then verify

  <\equation*>
    -<frac|k|\<delta\><rsub|0>>*u<rsub|0>+<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>><right|]>=<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>>-\<mu\><rsub|-1/2>*<frac|u<rsub|0>-u<rsub|-1>|\<delta\><rsub|0>+\<delta\><rsub|0>><right|]>
  </equation*>

  which gives

  <\equation*>
    u<rsub|-1>=<left|(>1-<frac|k*\<delta\><rsub|0>|\<mu\><rsub|-1/2>><right|)>*u<rsub|0>
  </equation*>

  Identifying with (<reference|ub>) we get

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<rsub|b>>|<cell|=>|<cell|0>>|<row|<cell|\<lambda\><rsub|b>>|<cell|=>|<cell|<frac|\<mu\><rsub|-1/2>|k*>*-<frac|\<delta\><rsub|0>|2>>>>>
  </eqnarray*>

  A more sensible discretisation seems to be

  <\with|mode|math>
    <\equation*>
      -<frac|k|\<delta\><rsub|0>>*(u<rsub|0>+u<rsub|-1>)/2+<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>><right|]>=<frac|2|\<delta\><rsub|0>>*<left|[>\<mu\><rsub|1/2>*<frac|u<rsub|1>-u<rsub|0>|\<delta\><rsub|0>+\<delta\><rsub|1>>-\<mu\><rsub|-1/2>*<frac|u<rsub|0>-u<rsub|-1>|\<delta\><rsub|0>+\<delta\><rsub|0>><right|]>
    </equation*>

    which gives

    <\equation*>
      u<rsub|-1>=<left|(><frac|2*\<mu\><rsub|-1/2>-k*\<delta\><rsub|0>|2*\<mu\><rsub|-1/2>+k*\<delta\><rsub|0>><right|)>*u<rsub|0>
    </equation*>

    Identifying with (<reference|ub>) we get

    <\eqnarray*>
      <tformat|<table|<row|<cell|u<rsub|b>>|<cell|=>|<cell|0>>|<row|<cell|\<lambda\><rsub|b>>|<cell|=>|<cell|<frac|\<mu\><rsub|-1/2>|k>*>>>>
    </eqnarray*>
  </with>

  which behaves as we expect i.e. <math|\<lambda\><rsub|b>\<rightarrow\>0>
  when <math|k\<rightarrow\>\<infty\>> (in contrast with the previous
  discretisation).
</body>

<\initial>
  <\collection>
    <associate|sfactor|5>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|1.1|?>>
    <associate|auto-3|<tuple|1.2|?>>
    <associate|ub|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Vertical
      diffusion> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1.5fn>|Bottom friction
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|Link between Navier condition and bottom
      friction <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>
    </associate>
  </collection>
</auxiliary>