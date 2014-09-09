<TeXmacs|1.0.7.3>

<style|article>

<\body>
  The viscous terms require the expression of

  <\equation*>
    <with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)
  </equation*>

  with

  <\equation*>
    2*<with|math-font-series|bold|D>=<with|math-font-series|bold|\<nabla\>u>+<with|math-font-series|bold|\<nabla\>><rsup|T><with|math-font-series|bold|u>
  </equation*>

  in general orthogonal coordinates. This can be expressed in the general
  case (Germain and Muller, chapter XIII.5)

  <\equation*>
    <with|math-font-series|bold|D>=<matrix|<tformat|<cwith|1|-1|1|-1|cell-lsep|5pt>|<cwith|1|-1|1|-1|cell-rsep|5pt>|<cwith|1|-1|1|-1|cell-bsep|5pt>|<cwith|1|-1|1|-1|cell-tsep|5pt>|<table|<row|<cell|<frac|u<rsub|1,1>|h<rsub|1>>+<frac|u<rsub|2>|h<rsub|2>>*<frac|h<rsub|1,2>|h<rsub|1>>>|<cell|<frac|1|2>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)>>>|<row|<cell|symmetric>|<cell|<frac|u<rsub|2,2>|h<rsub|2>>+<frac|u<rsub|1>|h<rsub|1>>*<frac|h<rsub|2,1>|h<rsub|2>>>>>>>
  </equation*>

  and

  <\equation*>
    <with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)=<frac|2|h<rsub|1>*h<rsub|2>>*<matrix|<tformat|<table|<row|<cell|(\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(\<eta\>*h<rsub|1>*D<rsub|12>)<rsub|,2>>>|<row|<cell|(\<eta\>*h<rsub|1>*D<rsub|22>)<rsub|,2>+(\<eta\>*h<rsub|2>*D<rsub|21>)<rsub|,1>>>>>>+<frac|2*\<eta\>|h<rsub|1>*h<rsub|2>>*<matrix|<tformat|<table|<row|<cell|D<rsub|12>*h<rsub|1,2>-D<rsub|22>*h<rsub|2,1>>>|<row|<cell|D<rsub|21>*h<rsub|2,1>-D<rsub|11>h<rsub|1,2>>>>>>
  </equation*>

  The divergence-free condition is written

  <\equation*>
    h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|u>=(u<rsub|1>h<rsub|2>)<rsub|,1>+(u<rsub|2>*h<rsub|1>)<rsub|,2>=0
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>*D<rsub|12>)<rsub|,2>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)><right|]><rsub|,2>>>>>
  </eqnarray*>

  Identities useful in the general case

  <\equation*>
    (\<eta\>*u<rsub|2,1>)<rsub|,2>=\<eta\><rsub|,2>*u<rsub|2,1>+\<eta\>*u<rsub|2,1,2>
  </equation*>

  <\equation*>
    (\<eta\>*u<rsub|2,2>)<rsub|,1>=\<eta\><rsub|,1>*u<rsub|2,2>+\<eta\>*u<rsub|2,1,2>
  </equation*>

  <\equation*>
    (\<eta\>*u<rsub|2,1>)<rsub|,2>=\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>+(\<eta\>*u<rsub|2,2>)<rsub|,1>
  </equation*>

  <subsection|General case>

  <\equation*>
    h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|u>=(u<rsub|1>h<rsub|2>)<rsub|,1>+(u<rsub|2>*h<rsub|1>)<rsub|,2>=0
  </equation*>

  <\equation>
    h<rsub|2>*u<rsub|1,1>+h<rsub|1>*u<rsub|2,2>+u<rsub|1>*h<rsub|2,1>+u<rsub|2>*h<rsub|1,2>=0<label|incomp>
  </equation>

  Multiplying by <math|\<eta\>/h<rsub|1>> and taking the derivative gives

  <\equation*>
    <left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>>*<right|]><rsub|,1>+(\<eta\>*u<rsub|2,2>)<rsub|,1>+<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|2,1>+u<rsub|2>*h<rsub|1,2>|h<rsub|1>><right|)><right|]><rsub|,1>=0
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|-*<left|[>\<eta\>*h<rsub|1>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+(\<eta\>*u<rsub|2,1>)<rsub|,2>>>|<row|<cell|>|<cell|>|<cell|+2*<left|[>\<eta\>*h<rsub|1,2>*<frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>+(\<eta\>*u<rsub|2,2>)<rsub|,1>>>|<row|<cell|>|<cell|>|<cell|+2*<left|[>\<eta\>*h<rsub|1,2>*<frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>>>|<row|<cell|>|<cell|>|<cell|+2*<left|[>\<eta\>*h<rsub|1,2>*<frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|2,1>+u<rsub|2>*h<rsub|1,2>|h<rsub|1>><right|)><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|1,2>*<frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>-<left|[>\<eta\>**h<rsub|1,2>*<frac|u<rsub|1>*|h<rsub|2>><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|-u<rsub|1>*<left|[>\<eta\>*<frac|h<rsub|2,1>|h<rsub|1>><right|]><rsub|,1>-h<rsub|2,1>*<left|[>\<eta\>*<frac|u<rsub|1>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|-u<rsub|2>*<left|[>\<eta\>*<frac|h<rsub|2,1>|h<rsub|2>><right|]><rsub|,2>-h<rsub|2,1>*<left|[>\<eta\>*<frac|u<rsub|2>*|h<rsub|2>><right|]><rsub|,2>>>>>
  </eqnarray*>

  <subsection|Decoupled case (<math|h<rsub|1,2>=h<rsub|2,1>=0>)>

  <\equation*>
    h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|u>=(u<rsub|1>h<rsub|2>)<rsub|,1>+(u<rsub|2>*h<rsub|1>)<rsub|,2>=0
  </equation*>

  <\equation*>
    h<rsub|2>*u<rsub|1,1>+h<rsub|1>*u<rsub|2,2>=0
  </equation*>

  Multiplying by <math|\<eta\>/h<rsub|1>> and taking the derivative gives

  <\equation*>
    <left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>>*<right|]><rsub|,1>+(\<eta\>*u<rsub|2,2>)<rsub|,1>=0
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+(\<eta\>*u<rsub|2,1>)<rsub|,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>+(\<eta\>*u<rsub|2,2>)<rsub|,1>>>|<row|<cell|>|<cell|=>|<cell|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>>>|<row|<cell|>|<cell|=>|<cell|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,1>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>>+\<eta\><rsub|,2>**h<rsub|2>*<frac|u<rsub|2,1>|h<rsub|2>>>>>>
  </eqnarray*>

  <subsection|Orthonormal case (<math|h<rsub|1>=h<rsub|2>=1>)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>>|<cell|=>|<cell|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>>>|<row|<cell|>|<cell|=>|<cell|(\<eta\>**u<rsub|1,1>)<rsub|,1>+(\<eta\>**u<rsub|1,2>)<rsub|,2>+\<eta\><rsub|,1>*u<rsub|1,1>+\<eta\><rsub|,2>*u<rsub|2,1>>>>>
  </eqnarray*>

  <subsection|Viscous term>

  <\eqnarray*>
    <tformat|<table|<row|<cell|h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|1>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>+2*\<eta\>**(D<rsub|12>*h<rsub|1,2>-D<rsub|22>*h<rsub|2,1>)>>|<row|<cell|>|<cell|=>|<cell|<with|color|red|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>>+\<eta\><rsub|,2>*u<rsub|2,1>-\<eta\><rsub|,1>*u<rsub|2,2>>>|<row|<cell|>|<cell|>|<cell|+2*<left|[>\<eta\>*h<rsub|1,2>*<frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|<with|color|blue|-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|2,1>+u<rsub|2>*h<rsub|1,2>|h<rsub|1>><right|)><right|]><rsub|,1>><eq-number>>>|<row|<cell|>|<cell|>|<cell|<with|color|orange|-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>>|<row|<cell|>|<cell|>|<cell|+<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)>*\<eta\>*h<rsub|1,2>>>|<row|<cell|>|<cell|>|<cell|<with|color|pink|-2*<left|(><frac|u<rsub|2,2>|h<rsub|2>>+<frac|u<rsub|1>|h<rsub|1>>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>*\<eta\>*h<rsub|2,1>>>>>>
  </eqnarray*>

  Regrouping terms then gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|1>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>+2*\<eta\>**(D<rsub|12>*h<rsub|1,2>-D<rsub|22>*h<rsub|2,1>)>>|<row|<cell|>|<cell|=>|<cell|<with|color|red|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>>>>|<row|<cell|>|<cell|>|<cell|<with|color|blue|-\<eta\><rsub|,1>*<left|(><frac|u<rsub|1>*h<rsub|2,1>-u<rsub|2>*h<rsub|1,2>|h<rsub|1>>+u<rsub|2,2><right|)>><with|color|orange|-\<eta\><rsub|,2>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>>-u<rsub|2,1><right|)>>>>|<row|<cell|>|<cell|>|<cell|<with|color|orange|-\<eta\>*u<rsub|1>*<left|[><left|(><frac|h<rsub|2,1>|h<rsub|1>><right|)><rsub|,1>+<left|(><frac|h<rsub|1,2>|h<rsub|2>><right|)><rsub|,2>+<frac|h<rsup|2><rsub|1,2>+2*h<rsup|2><rsub|2,1>|h<rsub|1>*h<rsub|2>>*<right|]>>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|2>*<left|[><left|(><frac|h<rsub|2,1>|h<rsub|2>><right|)><rsub|,2>-<left|(><frac|h<rsub|1,2>|h<rsub|1>><right|)><rsub|,1>+<frac|h<rsub|2,1>*h<rsub|1,2>|h<rsub|1>*h<rsub|2>>*<right|]><eq-number>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*<left|(>u<rsub|1,1>*<frac|h<rsub|2,1>|h<rsub|1>><with|color|pink|+3*u<rsub|2,2>*<frac|h<rsub|2,1>|h<rsub|2>>*>-2*u<rsub|2,1><frac|h<rsub|1,2>*|h<rsub|1>><right|)>>>>>
  </eqnarray*>

  Using the incompressibility condition (<reference|incomp>) one gets

  <\eqnarray*>
    <tformat|<table|<row|<cell|u<rsub|1,1>*<frac|h<rsub|2,1>|h<rsub|1>><with|color|pink|+3*u<rsub|2,2>*<frac|h<rsub|2,1>|h<rsub|2>>*>-2*u<rsub|2,1><frac|h<rsub|1,2>*|h<rsub|1>>>|<cell|=>|<cell|-u<rsub|1>*<frac|h<rsup|2><rsub|2,1>|h<rsub|2>*h<rsub|1>>-u<rsub|2>*<frac|h<rsub|2,1>*h<rsub|1,2>|h<rsub|2>*h<rsub|1>><with|color|pink|>-2*u<rsub|2,1><frac|h<rsub|1,2>*|h<rsub|1>>+2*u<rsub|2,2>*<frac|h<rsub|2,1>|h<rsub|2>>>>>>
  </eqnarray*>

  Replacing in the above and regrouping terms then gives

  <\eqnarray*>
    <tformat|<table|<row|<cell|h<rsub|1>*h<rsub|2>*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|1>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|1>D<rsub|12>)<rsub|,2>+2*\<eta\>**(D<rsub|12>*h<rsub|1,2>-D<rsub|22>*h<rsub|2,1>)>>|<row|<cell|>|<cell|=>|<cell|<with|color|red|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>>>>|<row|<cell|>|<cell|>|<cell|<with|color|blue|-\<eta\><rsub|,1>*<left|(><frac|u<rsub|1>*h<rsub|2,1>-u<rsub|2>*h<rsub|1,2>|h<rsub|1>>+u<rsub|2,2><right|)>><with|color|orange|-\<eta\><rsub|,2>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>>-u<rsub|2,1><right|)>>>>|<row|<cell|>|<cell|>|<cell|<with|color|orange|-\<eta\>*u<rsub|1>*<left|[><left|(><frac|h<rsub|2,1>|h<rsub|1>><right|)><rsub|,1>+<left|(><frac|h<rsub|1,2>|h<rsub|2>><right|)><rsub|,2>+<frac|h<rsup|2><rsub|1,2>+h<rsup|2><rsub|2,1>|h<rsub|1>*h<rsub|2>>*<right|]>>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|2>*<left|[><left|(><frac|h<rsub|2,1>|h<rsub|2>><right|)><rsub|,2>-<left|(><frac|h<rsub|1,2>|h<rsub|1>><right|)><rsub|,1><right|]><eq-number>>>|<row|<cell|>|<cell|>|<cell|+2*\<eta\>*<left|(>u<rsub|2,1>*<frac|h<rsub|1,2>*|h<rsub|1>>-u<rsub|2,2>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>>>>>
  </eqnarray*>

  <subsection|Application to polar coordinates>

  <\eqnarray*>
    <tformat|<table|<row|<cell|y<rsup|1>=r>|<cell|>|<cell|y<rsup|2>=\<theta\>>>|<row|<cell|h<rsub|1>=1>|<cell|>|<cell|h<rsub|2>=r>>|<row|<cell|h<rsub|1,2>=0>|<cell|>|<cell|h<rsub|2,1>=1>>|<row|<cell|>|<cell|\<eta\>=cte>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|r*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|r>>|<cell|=>|<cell|<with|color|red|<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>>>>|<row|<cell|>|<cell|>|<cell|<with|color|orange|-\<eta\>*<frac|u<rsub|r>|r>>-2*\<eta\>**<frac|u<rsub|\<theta\>,\<theta\>>|r>>>|<row|<cell|>|<cell|=>|<cell|\<eta\>*<left|[>(r*u<rsub|r,r>)<rsub|,r>+<frac|1|r>*u<rsub|r,\<theta\>,\<theta\>>-2*<frac|u<rsub|\<theta\>,\<theta\>>|r>-*<frac|u<rsub|r>|r><right|]>>>|<row|<cell|r*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|\<theta\>>>|<cell|=>|<cell|<with|color|red|<left|[>\<eta\>*h<rsub|1>*<frac|u<rsub|2,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|2,1>|h<rsub|1>><right|]><rsub|,1>>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>**<frac|u<rsub|\<theta\>>|r>*+2*\<eta\>**<frac|u<rsub|r,\<theta\>>*|r>>>|<row|<cell|>|<cell|=>|<cell|\<eta\>*<left|[>(r*u<rsub|\<theta\>,r>)<rsub|,r>+<frac|1|r>*u<rsub|\<theta\>,\<theta\>,\<theta\>>+2*<frac|u<rsub|r,\<theta\>>|r>-*<frac|u<rsub|\<theta\>>|r><right|]>>>>>
  </eqnarray*>

  which matches the result in Germain and Muller, chapter XIII.8.1.

  \;
</body>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-2|<tuple|2|?>>
    <associate|auto-3|<tuple|3|?>>
    <associate|auto-4|<tuple|4|?>>
    <associate|auto-5|<tuple|5|?>>
    <associate|auto-6|<tuple|6|?>>
    <associate|incomp|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <with|par-left|<quote|1.5fn>|1<space|2spc>General case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1>>

      <with|par-left|<quote|1.5fn>|2<space|2spc>Decoupled case
      (<with|mode|<quote|math>|h<rsub|1,2>=h<rsub|2,1>=0>)
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|3<space|2spc>Orthonormal case
      (<with|mode|<quote|math>|h<rsub|1>=h<rsub|2>=1>)
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1.5fn>|4<space|2spc>Viscous term
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <with|par-left|<quote|1.5fn>|5<space|2spc>Application to polar
      coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>
    </associate>
  </collection>
</auxiliary>