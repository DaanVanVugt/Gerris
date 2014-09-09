<TeXmacs|1.0.6.11>

<style|article>

<\body>
  The incompressible Navier--Stokes equations in cylindrical coordinates are

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t>v<rsub|r>+<frac|1|r>*\<partial\><rsub|r>(r*v<rsub|r><rsup|2>)+\<partial\><rsub|z>(v<rsub|r*>*v<rsub|z>)=-\<partial\><rsub|r>\<phi\>+<frac|1|r>*\<partial\><rsub|r>(r*S<rsub|rr>)+\<partial\><rsub|z>S<rsub|zr>-<frac|S<rsub|\<theta\>\<theta\>>|r>,>|<cell|<eq-number><label|momr>>>|<row|<cell|>|<cell|\<partial\><rsub|t>v<rsub|z>+<frac|1|r>*\<partial\><rsub|r>(r*v<rsub|r>*v<rsub|z>)+\<partial\><rsub|z>(v<rsup|2><rsub|z>)=-\<partial\><rsub|z>\<phi\>+<frac|1|r>*\<partial\><rsub|r>(r*S<rsub|zr>)+\<partial\><rsub|z>S<rsub|zz>,>|<cell|<eq-number><label|momz>>>|<row|<cell|>|<cell|<frac|1|r>*\<partial\><rsub|r>(r*v<rsub|r>)+\<partial\><rsub|z>v<rsub|z>=0,>|<cell|<eq-number><label|continuity>>>>>
  </eqnarray*>

  with <math|\<phi\>=p/\<rho\>> and the stress tensor

  <\eqnarray*>
    <tformat|<table|<row|<cell|S<rsub|rr>>|<cell|=>|<cell|2*\<nu\>*\<partial\><rsub|r>v<rsub|r>,>>|<row|<cell|S<rsub|\<theta\>\<theta\>>>|<cell|=>|<cell|2*\<nu\>*<frac|v<rsub|r>|r>,>>|<row|<cell|S<rsub|zz>>|<cell|=>|<cell|2*\<nu\>*\<partial\><rsub|z>v<rsub|z>,>>|<row|<cell|S<rsub|zr>>|<cell|=>|<cell|\<nu\>*<left|(>\<partial\><rsub|r>v<rsub|z>+\<partial\><rsub|z>v<rsub|r><right|)>.>>>>
  </eqnarray*>

  Considering a control volume <math|\<Omega\>> with boundary
  <math|\<partial\>\<Omega\>>, the integral equations can then be written

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t><big|int><rsub|\<Omega\>>v<rsub|r>*r*dr*dz+<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|r>*v<rsub|z>*dr+<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|r>*v<rsub|r>*dz=>|<cell|>>|<row|<cell|>|<cell|-<big|int><rsub|\<partial\>\<Omega\>>\<phi\>*r*dz+<big|int><rsub|\<Omega\>>\<phi\>*dr*dz+<big|int><rsub|\<partial\>\<Omega\>>r*S<rsub|zr>*dr+<big|int><rsub|\<partial\>\<Omega\>>r*S<rsub|rr>*dz-<big|int><rsub|\<Omega\>>S<rsub|\<theta\>\<theta\>>*dr*dz,>|<cell|>>|<row|<cell|>|<cell|\<partial\><rsub|t><big|int><rsub|\<Omega\>>v<rsub|z>*r*dr*dz+<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|z>*v<rsub|z>*dr+<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|z>*v<rsub|r>*dz=>|<cell|>>|<row|<cell|>|<cell|-<big|int><rsub|\<partial\>\<Omega\>>\<phi\>*r*dr+<big|int><rsub|\<partial\>\<Omega\>>r*S<rsub|zz>*dr+<big|int><rsub|\<partial\>\<Omega\>>r*S<rsub|zr>*dz,>|<cell|>>|<row|<cell|>|<cell|<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|z>*dr+<big|int><rsub|\<partial\>\<Omega\>>r*v<rsub|r>*dz=0.>|<cell|>>|<row|<cell|>|<cell|>|<cell|>>>>
  </eqnarray*>

  Posing

  <\eqnarray*>
    <tformat|<table|<row|<cell|<with|mode|text|<math|a<rsup|i,j>>>>|<cell|\<equiv\>>|<cell|<big|int><rsub|\<Omega\>>dr*dz,>>|<row|<cell|c<rsup|i,j>>|<cell|\<equiv\>>|<cell|<big|int><rsub|\<Omega\>>r*dr*dz\<approx\>r<rsup|j>*a<rsup|i,j>,>>|<row|<cell|s<rsub|z><rsup|i,j>>|<cell|\<equiv\>>|<cell|<big|int><rsub|\<partial\>\<Omega\><rsup|j>>r*dr\<approx\>r<rsup|j>*<big|int><rsub|\<partial\>\<Omega\><rsup|j>>dr,>>|<row|<cell|s<rsup|i,j-1/2><rsub|r>>|<cell|\<equiv\>>|<cell|<big|int><rsub|\<partial\>\<Omega\><rsup|j-1/2>>r*dz<rsub|>=r<rsup|j-1/2>*<big|int><rsub|\<partial\>\<Omega\><rsup|j-1/2>>dz,>>>>
  </eqnarray*>

  we get the discrete form of the equations as

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t>(c*v<rsub|r>)<rsup|i,j>+(s<rsub|z>*v<rsub|r>*v<rsub|z>)<rsup|i+1/2,j>-(s<rsub|z>*v<rsub|r>*v<rsub|z>)<rsup|i-1/2,j>+(s<rsub|r>*v<rsub|r><rsup|2>)<rsup|i,j+1/2>-(s<rsub|r>*v<rsub|r><rsup|2>)<rsup|i,j-1/2>=>|<cell|>>|<row|<cell|>|<cell|(s<rsub|r>\<phi\>)<rsup|i,j-1/2>-(s<rsub|r>\<phi\>)<rsup|i,j+1/2>+<with|color|blue|<with|math-font-series|bold|(a*\<phi\>)<rsup|i,j>>>+>|<cell|>>|<row|<cell|>|<cell|(s<rsub|z>*S<rsub|zr>)<rsup|i+1/2,j>-(s<rsub|z>*S<rsub|zr>)<rsup|i-1/2,j>+(s<rsub|r>*S<rsub|rr>)<rsup|i,j+1/2>-(s<rsub|r>*S<rsub|rr>)<rsup|i,j-1/2>-<with|color|blue|<with|math-font-series|bold|(a*S<rsub|\<theta\>\<theta\>>)<rsup|i,j>>>,>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|\<partial\><rsub|t>(c*v<rsub|z>)<rsup|i,j>+(s<rsub|z>*v<rsup|2><rsub|z>)<rsup|i+1/2,j>-(s<rsub|z>*v<rsup|2><rsub|z>)<rsup|i-1/2,j>+(s<rsub|r>*v<rsub|r><rsup|>*v<rsub|z>)<rsup|i,j+1/2>-(s<rsub|r>*v<rsub|r><rsup|>*v<rsub|z>)<rsup|i,j-1/2>=>|<cell|>>|<row|<cell|>|<cell|(s<rsub|z>\<phi\>)<rsup|i-1/2,j>-(s<rsub|z>\<phi\>)<rsup|i+1/2,j>+>|<cell|>>|<row|<cell|>|<cell|(s<rsub|z>*S<rsub|zz>)<rsup|i+1/2,j>-(s<rsub|z>*S<rsub|zz>)<rsup|i-1/2,j>+(s<rsub|r>*S<rsub|zr>)<rsup|i,j+1/2>-(s<rsub|r>*S<rsub|zr>)<rsup|i,j-1/2>>|<cell|>>>>
  </eqnarray*>

  <\equation*>
    (s<rsub|z>*v<rsub|z>)<rsup|i+1/2,j>-(s<rsub|z>*v<rsub|z>)<rsup|i-1/2,j>+(s<rsub|r>*v<rsub|r>)<rsup|i,j+1/2>-(s<rsub|r>*v<rsub|r>)<rsup|i,j-1/2>=0,
  </equation*>

  where only the bold terms differ from the discretisation of the N--S
  equations in Cartesian coordinates in \ two dimensions.

  The projection method relies on a decomposition of the velocity as

  <\eqnarray*>
    <tformat|<table|<row|<cell|v<rsup|i,j><rsub|r>>|<cell|=>|<cell|(v<rsub|r><rsup|\<star\>>)<rsup|i,j>+<frac|\<Delta\>t|c<rsup|i,j>><left|[>(s<rsub|r>\<phi\>)<rsup|i,j-1/2>-(s<rsub|r>\<phi\>)<rsup|i,j+1/2>+<with|color|blue|<with|math-font-series|bold|(a*\<phi\>)<rsup|i,j>>><right|]>,>>|<row|<cell|v<rsup|i,j><rsub|z>>|<cell|=>|<cell|(v<rsub|z><rsup|\<star\>>)<rsup|i,j>+<frac|\<Delta\>t|c<rsup|i,j>><left|[>(s<rsub|z>\<phi\>)<rsup|i-1/2,j>-(s<rsub|z>\<phi\>)<rsup|i+1/2,j><right|]>,>>>>
  </eqnarray*>

  wich leads to the Poisson-like equation

  <\equation*>
    <frac|1|\<Delta\>t>*<left|[>s<rsup|i+1/2,j><rsub|z>**(v<rsub|z><rsup|\<star\>>)<rsup|i+1/2,j>-s<rsup|i-1/2,j><rsub|z>*(v<rsub|z><rsup|\<star\>>)<rsup|i-1/2,j>+s<rsup|i,j+1/2><rsub|r>*(v<rsub|r><rsup|\<star\>>)<rsup|i,j+1/2>-s<rsup|i,j-1/2><rsub|r>*(v<rsub|r><rsup|\<star\>>)<rsup|i,j-1/2><right|]>+\<phi\><rsup|i,j>*<left|(>s<rsup|i,j><rsub|z>*<left|[>(s<rsub|z>/c)<rsup|i+1/2,j>+(s<rsub|z>/c)<rsup|i-1/2,j><right|]>+s<rsup|i,j><rsub|r>*<left|[>(s<rsub|r>/c)<rsup|i,j+1/2>+(s<rsub|r>/c)<rsup|i,j-1/2><right|]><right|)>+s<rsup|i,j+1/2><rsub|r>*(a*\<phi\>)<rsup|i,j+1/2>-s<rsup|i,j-1/2><rsub|r>*(a*\<phi\>)<rsup|i,j-1/2>-(s<rsub|z>/c)<rsup|i+1/2,j>*(s<rsub|z>\<phi\>)<rsup|i+1,j>-(s<rsub|z>/c)<rsup|i-1/2,j>*(s<rsub|z>\<phi\>)<rsup|i-1,j>-(s<rsub|r>/c)<rsup|i,j+1/2>*(s<rsub|r>\<phi\>)<rsup|i,j+1>-(s<rsub|r>/c)<rsup|i,j-1/2>*(s<rsub|r>\<phi\>)<rsup|i,j-1>=0.
  </equation*>

  This looks complicated though... What about just splitting the velocity as

  <\eqnarray*>
    <tformat|<table|<row|<cell|v<rsup|i,j><rsub|r>>|<cell|=>|<cell|(v<rsub|r><rsup|\<star\>>)<rsup|i,j>+\<Delta\>t*<left|(>\<phi\><rsup|i,j-1/2>-\<phi\><rsup|i,j+1/2><right|)>,>>|<row|<cell|v<rsup|i,j><rsub|z>>|<cell|=>|<cell|(v<rsub|z><rsup|\<star\>>)<rsup|i,j>+\<Delta\>t*<left|(>\<phi\><rsup|i-1/2,j>-\<phi\><rsup|i+1/2,j><right|)>,>>>>
  </eqnarray*>

  which is still consistent with equations (<reference|momr>) and
  (<reference|momz>). This gives

  <\equation*>
    <frac|1|\<Delta\>t>*<left|[>s<rsup|i+1/2,j><rsub|z>**(v<rsub|z><rsup|\<star\>>)<rsup|i+1/2,j>-s<rsup|i-1/2,j><rsub|z>*(v<rsub|z><rsup|\<star\>>)<rsup|i-1/2,j>+s<rsup|i,j+1/2><rsub|r>*(v<rsub|r><rsup|\<star\>>)<rsup|i,j+1/2>-s<rsup|i,j-1/2><rsub|r>*(v<rsub|r><rsup|\<star\>>)<rsup|i,j-1/2><right|]>+\<phi\><rsup|i,j>*(s<rsup|i+1/2,j><rsub|z>+s<rsup|i-1/2,j><rsub|z>+s<rsup|i,j+1/2><rsub|r>+s<rsup|i,j-1/2><rsub|r>)-s<rsub|z><rsup|i+1/2,j>*\<phi\><rsup|i+1,j>-s<rsub|z><rsup|i-1/2,j>*\<phi\><rsup|i-1,j>-s<rsub|r><rsup|i,j+1/2>*\<phi\><rsup|i,j+1>-s<rsub|r><rsup|i,j-1/2>*\<phi\><rsup|i,j-1>=0,
  </equation*>

  together with

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c<rsup|i,j>*<frac|(v<rsub|r><rsup|\<star\>>)<rsup|i,j>-v<rsup|i,j><rsub|r>|\<Delta\>t>+(s<rsub|z>*v<rsub|r>*v<rsub|z>)<rsup|i+1/2,j>-(s<rsub|z>*v<rsub|r>*v<rsub|z>)<rsup|i-1/2,j>+(s<rsub|r>*v<rsub|r><rsup|2>)<rsup|i,j+1/2>-(s<rsub|r>*v<rsub|r><rsup|2>)<rsup|i,j-1/2>=>|<cell|>>|<row|<cell|>|<cell|(s<rsub|z>*S<rsub|zr>)<rsup|i+1/2,j>-(s<rsub|z>*S<rsub|zr>)<rsup|i-1/2,j>+(s<rsub|r>*S<rsub|rr>)<rsup|i,j+1/2>-(s<rsub|r>*S<rsub|rr>)<rsup|i,j-1/2>-<with|color|blue|<with|math-font-series|bold|(a*S<rsub|\<theta\>\<theta\>>)<rsup|i,j>>>,>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c<rsup|i,j>*<frac|(v<rsub|z><rsup|\<star\>>)<rsup|i,j>-v<rsup|i,j><rsub|z>|\<Delta\>t>+(s<rsub|z>*v<rsup|2><rsub|z>)<rsup|i+1/2,j>-(s<rsub|z>*v<rsup|2><rsub|z>)<rsup|i-1/2,j>+(s<rsub|r>*v<rsub|r><rsup|>*v<rsub|z>)<rsup|i,j+1/2>-(s<rsub|r>*v<rsub|r><rsup|>*v<rsub|z>)<rsup|i,j-1/2>=>|<cell|>>|<row|<cell|>|<cell|(s<rsub|z>*S<rsub|zz>)<rsup|i+1/2,j>-(s<rsub|z>*S<rsub|zz>)<rsup|i-1/2,j>+(s<rsub|r>*S<rsub|zr>)<rsup|i,j+1/2>-(s<rsub|r>*S<rsub|zr>)<rsup|i,j-1/2>.>|<cell|>>>>
  </eqnarray*>

  <subsection|Advection term>

  <\equation*>
    u<rsup|n+1/2><rsub|d>=u<rsup|n>+<frac|h|2>*\<partial\><rsub|d>u<rsup|n>+<frac|\<Delta\>t|2>*\<partial\><rsub|t>u<rsup|n>+\<cal-O\>(h<rsup|2>,\<Delta\>t<rsup|2>),
  </equation*>

  using equations (<reference|momr>), (<reference|momz>) and
  (<reference|continuity>) then gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|v<rsup|n+1/2><rsub|r>=v<rsup|n><rsub|r>+<frac|1|2>*(h*-v<rsub|r>*\<Delta\>t)*\<partial\><rsub|r>v<rsub|r><rsup|n>+<frac|\<Delta\>t|2>*(-v<rsup|n><rsub|z>*\<partial\><rsub|z>v<rsup|n><rsub|r*>+src<rsup|n>),>|<cell|>|<cell|>>|<row|<cell|v<rsup|n+1/2><rsub|z>=v<rsub|z><rsup|n>+<frac|1|2>*(h-v<rsub|z>*\<Delta\>t)*\<partial\><rsub|z>v<rsub|z><rsup|n>+<frac|\<Delta\>t|2>*<left|(>-v<rsup|n><rsub|r>*\<partial\><rsub|r>v<rsup|n><rsub|z>+src<rsup|n><right|)>,>|<cell|>|<cell|>>>>
  </eqnarray*>

  which is the same as in the two-dimensional case.

  <subsection|Implicit diffusion>

  <subsubsection|Original implicit diffusion>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c*<frac|v<rsup|\<star\>>-v<rsup|n>|\<Delta\>t>-\<beta\>*\<nabla\>s*S<rsup|\<star\>><rsup|*>=src<rsup|n+1/2>+(1-\<beta\>)*\<nabla\>s*S<rsup|n>>|<cell|>>>>
  </eqnarray*>

  which can be rewritten

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c*v<rsup|\<star\>>-\<nabla\><wide|s|^>*S<rsup|\<star\>>=c*v<rsup|n>+\<Delta\>t*src<rsup|n+1/2>+<frac|1-\<beta\>|\<beta\>>*\<nabla\><wide|s|^>*S<rsup|n>,>|<cell|>>>>
  </eqnarray*>

  with

  <\equation*>
    <wide|s|^>\<equiv\>\<beta\>*\<Delta\>t*s.
  </equation*>

  This finally gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|v<rsup|\<star\>>-<frac|1|c>\<nabla\><wide|s|^>*S<rsup|\<star\>>=v<rsup|n>+\<Delta\>t*<wide|src|^><rsup|n+1/2>+<frac|<wide|\<beta\>|^>|c>*\<nabla\><wide|s|^>*S<rsup|n>,>|<cell|>>>>
  </eqnarray*>

  with <math|><with|mode|math|<wide|\<beta\>|^>\<equiv\>(1-\<beta\>)/\<beta\>>.

  <subsubsection|Axisymmetric version>

  To account for the axisymmetric term this needs to be rewritten

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c*<frac|v<rsup|\<star\>>-v<rsup|n>|\<Delta\>t>-\<beta\>*\<nabla\>s*S<rsup|\<star\>>+\<beta\>*a*2*\<nu\>*<frac|v<rsup|\<star\>>|r>=src<rsup|n+1/2>+(1-\<beta\>)*\<nabla\>s*S<rsup|n>-(1-\<beta\>)*a*2*\<nu\>*<frac|v<rsup|n>|r>>|<cell|>>>>
  </eqnarray*>

  which can be rewritten

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|c*v<rsup|\<star\>>*-\<nabla\><wide|s|^>*S<rsup|\<star\>>+2*\<beta\>*\<Delta\>t*a*\<nu\>*<frac|v<rsup|\<star\>>|r>=c*v<rsup|n>+\<Delta\>t*src<rsup|n+1/2>+<frac|1-\<beta\>|\<beta\>>*\<nabla\><wide|s|^>*S<rsup|n>-2*(1-\<beta\>)*\<Delta\>t*a*\<nu\>*<frac|v<rsup|n>|r>>|<cell|>>>>
  </eqnarray*>

  this finally gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|v<rsup|\<star\>>*-<frac|1|c>*\<nabla\><wide|s|^>*S<rsup|\<star\>>+2*\<beta\>*\<Delta\>t*\<nu\>*<frac|v<rsup|\<star\>>|r<rsup|2>>=v<rsup|n>+\<Delta\>t*<wide|src|^><rsup|n+1/2>+<frac|1-\<beta\>|\<beta\>*c>*\<nabla\><wide|s|^>*S<rsup|n>-2*(1-\<beta\>)*\<Delta\>t*\<nu\>*<frac|v<rsup|n>|r<rsup|2>>>|<cell|>>>>
  </eqnarray*>

  Posing <math|d\<equiv\>2*\<beta\>*\<Delta\>t*\<nu\>/r<rsup|2>> then gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|v<rsup|\<star\>>*(1+d)-<frac|1|c>*\<nabla\><wide|s|^>*S<rsup|\<star\>>=v<rsup|n>*<left|(>1-<wide|\<beta\>*|^>*d<right|)>+\<Delta\>t*<wide|src|^><rsup|n+1/2>+<frac|<wide|\<beta\>*|^>*|c>*\<nabla\><wide|s|^>*S<rsup|n>.>|<cell|>>>>
  </eqnarray*>
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|2.1|?>>
    <associate|auto-4|<tuple|2.2|?>>
    <associate|auto-5|<tuple|3|?>>
    <associate|continuity|<tuple|3|?>>
    <associate|momr|<tuple|1|?>>
    <associate|momz|<tuple|2|?>>
    <associate|mon|<tuple|?|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <with|par-left|<quote|3fn>|1<space|2spc>Advection term
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>>
    </associate>
  </collection>
</auxiliary>