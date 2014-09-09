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
    <with|math-font-series|bold|D>=<matrix|<tformat|<cwith|1|-1|1|-1|cell-lsep|5pt>|<cwith|1|-1|1|-1|cell-rsep|5pt>|<cwith|1|-1|1|-1|cell-bsep|5pt>|<cwith|1|-1|1|-1|cell-tsep|5pt>|<table|<row|<cell|<frac|u<rsub|1,1>|h<rsub|1>>+<frac|u<rsub|2>|h<rsub|2>>*<frac|h<rsub|1,2>|h<rsub|1>>+<frac|u<rsub|3>|h<rsub|3>>*<frac|h<rsub|1,3>|h<rsub|1>>>|<cell|<frac|1|2>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)>>|<cell|<frac|1|2>*<left|(><frac|u<rsub|1,3>|h<rsub|3>>+<frac|u<rsub|3,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|1>*h<rsub|3>><right|)>>>|<row|<cell|symmetric>|<cell|circular
    permutation>|<cell|circular permutation>>|<row|<cell|symmetric>|<cell|symmetric>|<cell|circular
    permutation>>>>>
  </equation*>

  and

  <\equation*>
    <with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)=<frac|2|h<rsub|1>*h<rsub|2>*h<rsub|3>>*<big|sum><left|[>\<eta\>*h<rsub|2>*h<rsub|3>*(D<rsub|11>*\<b-e\><rsub|1>+D<rsub|21>*\<b-e\><rsub|2>+D<rsub|31>*\<b-e\><rsub|3>)<right|]><rsub|,1>
  </equation*>

  The first component can then be written

  <\eqnarray*>
    <tformat|<table|<row|<cell|h<rsub|1>*h<rsub|2>*h<rsub|3>*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|1>>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*h<rsub|3>*D<rsub|11>)<rsub|,1>*+(2*\<eta\>*h<rsub|3>*h<rsub|1>*D<rsub|12>)<rsub|,2>*+(2*\<eta\>*h<rsub|1>*h<rsub|2>*D<rsub|13>)<rsub|,3>*>>|<row|<cell|>|<cell|>|<cell|+2*\<eta\>**(D<rsub|21>*h<rsub|3>*h<rsub|1,2>+D<rsub|31>*h<rsub|2>*h<rsub|1,3>-D<rsub|22>*h<rsub|3>*h<rsub|2,1>-D<rsub|33>*h<rsub|2>*h<rsub|3,1>)>>>>
  </eqnarray*>

  The divergence-free condition is written

  <\equation*>
    h<rsub|1>*h<rsub|2>*h<rsub|3>*<with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|u>=<big|sum>(u<rsub|1>*h<rsub|2>*h<rsub|3>)<rsub|,1>=0
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|A\<equiv\>(2*\<eta\>*h<rsub|2>*h<rsub|3>*D<rsub|11>)<rsub|,1>+(2*\<eta\>*h<rsub|3>*h<rsub|1>*D<rsub|12>)<rsub|,2>+(2*\<eta\>*h<rsub|1>*h<rsub|2>*D<rsub|13>)<rsub|,3>*=>|<cell|>>|<row|<cell|>|<cell|(2*\<eta\>*h<rsub|2>*h<rsub|3>*D<rsub|11>)<rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<left|(><frac|u<rsub|1,3>|h<rsub|3>>+<frac|u<rsub|3,1>|h<rsub|1>>-<frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|1>*h<rsub|3>><right|)><right|]><rsub|,3>>|<cell|>>|<row|<cell|>|<cell|>|<cell|>>>>
  </eqnarray*>

  Identities useful in the general case

  <\equation*>
    (\<eta\>*h<rsub|3>*u<rsub|2,1>)<rsub|,2>=(\<eta\>*h<rsub|3>)<rsub|,2>*u<rsub|2,1>+\<eta\>*h<rsub|3>*u<rsub|2,1,2>
  </equation*>

  <\equation*>
    (\<eta\>*h<rsub|3>*u<rsub|2,2>)<rsub|,1>=(\<eta\>*h<rsub|3>)<rsub|,1>*u<rsub|2,2>+\<eta\>*h<rsub|3>*u<rsub|2,1,2>
  </equation*>

  <\equation*>
    (\<eta\>*h<rsub|3>*u<rsub|2,1>)<rsub|,2>=(\<eta\>*h<rsub|3>)<rsub|,2>*u<rsub|2,1>-(\<eta\>*h<rsub|3>)<rsub|,1>*u<rsub|2,2>+(\<eta\>*h<rsub|3>*u<rsub|2,2>)<rsub|,1>
  </equation*>

  <subsection|General case>

  <\equation*>
    h<rsub|1>*h<rsub|2>*h<rsub|3>*<with|math-font-series|bold|\<nabla\>>\<cdot\><with|math-font-series|bold|u>=(u<rsub|1>h<rsub|2>*h<rsub|3>)<rsub|,1>+(u<rsub|2>*h<rsub|3>*h<rsub|1>)<rsub|,2>+(u<rsub|3>*h<rsub|1>*h<rsub|2>)<rsub|,3>=0
  </equation*>

  <\equation>
    u<rsub|1>(h<rsub|2>*h<rsub|3>)<rsub|,1>+h<rsub|2>*h<rsub|3>*u<rsub|1,1>+u<rsub|2*>*(h<rsub|3>*h<rsub|1>)<rsub|,2>+h<rsub|3>*h<rsub|1>*u<rsub|2,2>+u<rsub|3>*(h<rsub|1>*h<rsub|2>)<rsub|,3>+h<rsub|1>*h<rsub|2>*u<rsub|3,3>=0<label|incomp>
  </equation>

  Multiplying by <math|\<eta\>/h<rsub|1>> and taking the derivative gives

  <\equation*>
    <left|[>\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>>*<right|]><rsub|,1>+(\<eta\>*h<rsub|3>*u<rsub|2,2>)<rsub|,1>+(\<eta\>*h<rsub|2>*u<rsub|3,3>)<rsub|,1>+<left|[>\<eta\>*<left|(><frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>+u<rsub|2>*(h<rsub|3>*h<rsub|1>)<rsub|,2>+u<rsub|3>*(h<rsub|1>*h<rsub|2>)<rsub|,3>|h<rsub|1>><right|)><right|]><rsub|,1>=0
  </equation*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|A>|<cell|=>|<cell|(2*\<eta\>*h<rsub|2>*h<rsub|3>*D<rsub|11>)<rsub|,1>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<left|(><frac|u<rsub|1,2>|h<rsub|2>>+<frac|u<rsub|2,1>|h<rsub|1>><right|)><right|]><rsub|,2>-*<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|1>*h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<left|(><frac|u<rsub|1,3>|h<rsub|3>>+<frac|u<rsub|3,1>|h<rsub|1>><right|)><right|]><rsub|,3>-<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|1>*h<rsub|3>><right|)><right|]><rsub|,3>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|3>*h<rsub|1,2>*<frac||><frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|2>*h<rsub|1,3>*<frac|u<rsub|3>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+(\<eta\>*h<rsub|3>*u<rsub|2,1>)<rsub|,2>-*<left|[>\<eta\>*h<rsub|3>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>+(\<eta\>*h<rsub|2>*u<rsub|3,1>)<rsub|,3>-<left|[>\<eta\>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|3>><right|)><right|]><rsub|,3>>>|<row|<cell|>|<cell|=>|<cell|<left|[>2*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|3>*h<rsub|1,2>*<frac||><frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|2>*h<rsub|1,3>*<frac|u<rsub|3>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+(\<eta\>*h<rsub|3>)<rsub|,2>*u<rsub|2,1>-(\<eta\>*h<rsub|3>)<rsub|,1>*u<rsub|2,2>+(\<eta\>*h<rsub|3>*u<rsub|2,2>)<rsub|,1>-*<left|[>\<eta\>*h<rsub|3>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>+(\<eta\>*h<rsub|2>)<rsub|,3>*u<rsub|3,1>-(\<eta\>*h<rsub|2>)<rsub|,1>*u<rsub|3,3>+(\<eta\>*h<rsub|2>*u<rsub|3,3>)<rsub|,1>-<left|[>\<eta\>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|3>><right|)><right|]><rsub|,3>>>|<row|<cell|>|<cell|=>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|3>*h<rsub|1,2>*<frac||><frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>+<left|[>2*\<eta\>*h<rsub|2>*h<rsub|1,3>*<frac|u<rsub|3>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+(\<eta\>*h<rsub|3>)<rsub|,2>*u<rsub|2,1>-(\<eta\>*h<rsub|3>)<rsub|,1>*u<rsub|2,2>-*<left|[>\<eta\>*h<rsub|3>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|)><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>+(\<eta\>*h<rsub|2>)<rsub|,3>*u<rsub|3,1>-(\<eta\>*h<rsub|2>)<rsub|,1>*u<rsub|3,3>-<left|[>\<eta\>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|3>><right|)><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-<left|[>\<eta\>*<left|(><frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>+u<rsub|2>*(h<rsub|3>*h<rsub|1>)<rsub|,2>+u<rsub|3>*(h<rsub|1>*h<rsub|2>)<rsub|,3>|h<rsub|1>><right|)><right|]><rsub|,1>>>|<row|<cell|>|<cell|=>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|+\<eta\><rsub|,2>*h<rsub|3>*u<rsub|2,1>-(\<eta\>*h<rsub|3>)<rsub|,1>*u<rsub|2,2>+\<eta\><rsub|,3>*h<rsub|2>*u<rsub|3,1>-(\<eta\>*h<rsub|2>)<rsub|,1>*u<rsub|3,3>>>|<row|<cell|>|<cell|>|<cell|+<left|[>\<eta\>*h<rsub|3>*h<rsub|1,2><frac|u<rsub|2>|h<rsub|1>><right|]><rsub|,1>-(\<eta\>*h<rsub|3,2>*)<rsub|,1>*u<rsub|2>+<left|[>\<eta\>*h<rsub|2>*h<rsub|1,3>*<frac|u<rsub|3>|h<rsub|1>><right|]><rsub|,1>-(\<eta\>*h<rsub|2,3>*)<rsub|,1>*u<rsub|3>>>|<row|<cell|>|<cell|>|<cell|-*<left|[>\<eta\>*h<rsub|3>*<frac|u<rsub|1>*h<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>-*<left|[>\<eta\>*h<rsub|3>*<frac|u<rsub|2>*h<rsub|2,1>|h<rsub|2>><right|]><rsub|,2>>>|<row|<cell|>|<cell|>|<cell|-<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|1>*h<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>-<left|[>\<eta\>*h<rsub|2>*<frac|u<rsub|3>*h<rsub|3,1>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-<left|[>\<eta\>*<frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>><right|]><rsub|,1>>>|<row|<cell|>|<cell|=>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,1>*<left|(>h<rsub|3>*u<rsub|2,2>+*h<rsub|2>*u<rsub|3,3>+h<rsub|3,2>*u<rsub|2>+h<rsub|2,3>*u<rsub|3>-h<rsub|3>*h<rsub|1,2><frac|u<rsub|2>|h<rsub|1>>-h<rsub|2>*h<rsub|1,3>*<frac|u<rsub|3>|h<rsub|1>>+<frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>><right|)>+\<eta\><rsub|,2>*<left|(>h<rsub|3>*u<rsub|2,1>-<frac|h<rsub|3>|h<rsub|2>>**(u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>)<right|)>+\<eta\><rsub|,3>*<left|(>h<rsub|2>*u<rsub|3,1>-<frac|h<rsub|2>|h<rsub|3>>**(u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>)<right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|1>*<left|[><left|(><frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>><right|)><rsub|,1>+<left|(>h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|2>><right|)><rsub|,2>+<left|(>h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|3>><right|)><rsub|,3><right|]>>>|<row|<cell|>|<cell|>|<cell|+\<eta\>*u<rsub|2>*<left|[><left|(>h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|1>><right|)><rsub|,1>-*<left|(>h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)><rsub|,2>-h<rsub|3,2,1><right|]>>>|<row|<cell|>|<cell|>|<cell|+\<eta\>*u<rsub|3>*<left|[><left|(>h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|1>><right|)><rsub|,1>-<left|(>h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)><rsub|,3>-h<rsub|2,3,1><right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|1,1>*<frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>>-\<eta\>**u<rsub|2,2>*<left|(>h<rsub|3,1>+h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>-\<eta\>*u<rsub|3,3>*<left|(>h<rsub|2,1>+h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)>>>|<row|<cell|>|<cell|>|<cell|-*\<eta\>*u<rsub|1,2>*h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|2>>-\<eta\>*u<rsub|1,3>*h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|3>>>>|<row|<cell|>|<cell|>|<cell|+\<eta\>*u<rsub|2,1>*h<rsub|1,2><frac|h<rsub|3>|h<rsub|1>>+\<eta\>*u<rsub|3,1>*h<rsub|1,3>*<frac|h<rsub|2>|h<rsub|1>>>>>>
  </eqnarray*>

  <subsection|Viscous term>

  <with|mode|math|A+2*\<eta\>**(D<rsub|21>*h<rsub|3>*h<rsub|1,2>+D<rsub|31>*h<rsub|2>*h<rsub|1,3>-D<rsub|22>*h<rsub|3>*h<rsub|2,1>-D<rsub|33>*h<rsub|2>*h<rsub|3,1>)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|=>|<cell|>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,1>*<left|(><frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>-h<rsub|3>*h<rsub|1,2>*u<rsub|2>-h<rsub|2>*h<rsub|1,3>*u<rsub|3>|h<rsub|1>>+(h<rsub|3>*u<rsub|2>)<rsub|,2>+*(h<rsub|2>*u<rsub|3>)<rsub|,3><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,2>*h<rsub|3>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>>*-u<rsub|2,1><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,3>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|3>>**-u<rsub|3,1><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|1>*<left|[><left|(><frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>><right|)><rsub|,1>+<left|(>h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|2>><right|)><rsub|,2>+<left|(>h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|3>><right|)><rsub|,3>+<frac|h<rsub|3>*(h<rsup|2><rsub|1,2>+2*h<rsup|2><rsub|2,1>)|h<rsub|1>*h<rsub|2>>*+<frac|h<rsub|2>*(h<rsup|2><rsub|1,3>+2*h<rsup|2><rsub|3,1>)|h<rsub|1>*h<rsub|3>><right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|2>*<left|[><left|(>h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)><rsub|,2>-<left|(>h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|1>><right|)><rsub|,1>+h<rsub|3,2,1>+<frac|h<rsub|3>*h<rsub|2,1>*h<rsub|1,2>|h<rsub|1>*h<rsub|2>>*+2*<frac|h<rsub|3,2>*h<rsub|3,1>|h<rsub|3>>*<right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|3>*<left|[><left|(>h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)><rsub|,3>-<left|(>h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|1>><right|)><rsub|,1>+h<rsub|2,3,1>+<frac|h<rsub|2>*h<rsub|3,1>*h<rsub|1,3>|h<rsub|1>*h<rsub|3>>+2*<frac|h<rsub|2,3>*h<rsub|2,1>|h<rsub|2>><right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*<left|[>u<rsub|1,1>*<frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>>+u<rsub|2,2>*<left|(>h<rsub|3,1>+3*h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>+*u<rsub|3,3>*<left|(>h<rsub|2,1>+3*h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)>>>|<row|<cell|>|<cell|>|<cell|
    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ -2*u<rsub|2,1>*h<rsub|1,2><frac|h<rsub|3>|h<rsub|1>>-2*u<rsub|3,1>*h<rsub|1,3>*<frac|h<rsub|2>|h<rsub|1>><right|]>>>>>
  </eqnarray*>

  Using the incompressibility condition (<reference|incomp>) one gets

  <\eqnarray*>
    <tformat|<table|<row|<cell|>|<cell|>|<cell|u<rsub|1,1>*<frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>>+u<rsub|2,2>*<left|(>h<rsub|3,1>+3*h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>+*u<rsub|3,3>*<left|(>h<rsub|2,1>+3*h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)>-2*u<rsub|2,1>*h<rsub|1,2><frac|h<rsub|3>|h<rsub|1>>-2*u<rsub|3,1>*h<rsub|1,3>*<frac|h<rsub|2>|h<rsub|1>>>>|<row|<cell|>|<cell|=>|<cell|<frac|-u<rsub|1>(h<rsub|2>*h<rsub|3>)<rsub|,1>-u<rsub|2*>*(h<rsub|3>*h<rsub|1>)<rsub|,2>-h<rsub|3>*h<rsub|1>*u<rsub|2,2>-u<rsub|3>*(h<rsub|1>*h<rsub|2>)<rsub|,3>-h<rsub|1>*h<rsub|2>*u<rsub|3,3>|h<rsub|2>*h<rsub|3>>*<frac|(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>>>>|<row|<cell|>|<cell|>|<cell|+u<rsub|2,2>*<left|(>h<rsub|3,1>+3*h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>><right|)>+*u<rsub|3,3>*<left|(>h<rsub|2,1>+3*h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>><right|)>-2*u<rsub|2,1>*h<rsub|1,2><frac|h<rsub|3>|h<rsub|1>>-2*u<rsub|3,1>*h<rsub|1,3>*<frac|h<rsub|2>|h<rsub|1>>>>|<row|<cell|>|<cell|=>|<cell|-u<rsub|1>*<frac|<left|[>(h<rsub|2>*h<rsub|3>)<rsub|,1><right|]><rsup|2>|h<rsub|1>*h<rsub|2>*h<rsub|3>>*-u<rsub|2>*<frac|(h<rsub|3>*h<rsub|1>)<rsub|,2>*(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>*h<rsub|2>*h<rsub|3>>-u<rsub|3>*<frac|(h<rsub|1>*h<rsub|2>)<rsub|,3>*(h<rsub|2>*h<rsub|3>)<rsub|,1>|h<rsub|1>*h<rsub|2>*h<rsub|3>>>>|<row|<cell|>|<cell|>|<cell|+2*u<rsub|2,2>*h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>>+2*u<rsub|3,3>*h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>>-2*u<rsub|2,1>*h<rsub|1,2><frac|h<rsub|3>|h<rsub|1>>-2*u<rsub|3,1>*h<rsub|1,3>*<frac|h<rsub|2>|h<rsub|1>>>>>>
  </eqnarray*>

  Replacing in the above and regrouping terms then gives

  <with|mode|math|A+2*\<eta\>**(D<rsub|21>*h<rsub|3>*h<rsub|1,2>+D<rsub|31>*h<rsub|2>*h<rsub|1,3>-D<rsub|22>*h<rsub|3>*h<rsub|2,1>-D<rsub|33>*h<rsub|2>*h<rsub|3,1>)>

  <\eqnarray*>
    <tformat|<table|<row|<cell|=>|<cell|>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,1>*<left|(><frac|u<rsub|1>*(h<rsub|2>*h<rsub|3>)<rsub|,1>-h<rsub|3>*h<rsub|1,2>*u<rsub|2>-h<rsub|2>*h<rsub|1,3>*u<rsub|3>|h<rsub|1>>+(h<rsub|3>*u<rsub|2>)<rsub|,2>+*(h<rsub|2>*u<rsub|3>)<rsub|,3><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,2>*h<rsub|3>*<left|(><frac|u<rsub|1>*h<rsub|1,2>+u<rsub|2>*h<rsub|2,1>|h<rsub|2>>*-u<rsub|2,1><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\><rsub|,3>*h<rsub|2>*<left|(><frac|u<rsub|1>*h<rsub|1,3>+u<rsub|3>*h<rsub|3,1>|h<rsub|3>>**-u<rsub|3,1><right|)>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|1>*<left|[>h<rsub|3>*<left|(><frac|h<rsub|2,1>|h<rsub|1>><right|)><rsub|,1>+h<rsub|2>**<left|(><frac|h<rsub|3,1>|h<rsub|1>><right|)><rsub|,1>+<left|(>h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|2>><right|)><rsub|,2>+<left|(>h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|3>><right|)><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|
    \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ +h<rsub|3>*<frac|h<rsup|2><rsub|1,2>+h<rsup|2><rsub|2,1>|h<rsub|1>*h<rsub|2>>*+h<rsub|2>*<frac|h<rsup|2><rsub|1,3>+h<rsup|2><rsub|3,1>|h<rsub|1>*h<rsub|3>><right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|2>*<left|[>h<rsub|3>*<left|(><frac|h<rsub|2,1>|h<rsub|2>><right|)><rsub|,2>-h<rsub|3>*<left|(><frac|h<rsub|1,2>|h<rsub|1>><right|)><rsub|,1>+h<rsub|3,2,1>+h<rsub|3,1>*<left|(><frac|h<rsub|3,2>*|h<rsub|3>>*-2*<frac|h<rsub|1,2>*|h<rsub|1>><right|)><right|]>>>|<row|<cell|>|<cell|>|<cell|-\<eta\>*u<rsub|3>*<left|[>h<rsub|2>*<left|(><frac|h<rsub|3,1>|h<rsub|3>><right|)><rsub|,3>-h<rsub|2>*<left|(><frac|h<rsub|1,3>|h<rsub|1>><right|)><rsub|,1>+h<rsub|2,3,1>+h<rsub|2,1>*<left|(><frac|h<rsub|2,3>*|h<rsub|2>>-2*<frac|h<rsub|1,3>|h<rsub|1>><right|)><right|]>>>|<row|<cell|>|<cell|>|<cell|-2*\<eta\>*<left|[>u<rsub|2,2>*h<rsub|3>*<frac|h<rsub|2,1>|h<rsub|2>>+u<rsub|3,3>*h<rsub|2>*<frac|h<rsub|3,1>|h<rsub|3>>-u<rsub|2,1>*h<rsub|3>*<frac|h<rsub|1,2>|h<rsub|1>>-u<rsub|3,1>*h<rsub|2>*<frac|h<rsub|1,3>|h<rsub|1>><right|]>>>>>
  </eqnarray*>

  <subsection|Application to spherical coordinates>

  <\eqnarray*>
    <tformat|<table|<row|<cell|y<rsup|1>=r>|<cell|y<rsup|2>=\<phi\>>|<cell|y<rsup|3>=\<theta\>>>|<row|<cell|h<rsub|1>=1>|<cell|h<rsub|2>=r>|<cell|h<rsub|3>=r*sin\<phi\>>>|<row|<cell|h<rsub|1,2>=0>|<cell|h<rsub|1,3>=0>|<cell|>>|<row|<cell|h<rsub|2,1>=1>|<cell|h<rsub|2,3>=0>|<cell|>>|<row|<cell|h<rsub|3,1>=sin\<phi\>>|<cell|h<rsub|3,2>=r*cos\<phi\>>|<cell|>>|<row|<cell|>|<cell|\<eta\>=cte>|<cell|>>>>
  </eqnarray*>

  <\eqnarray*>
    <tformat|<table|<row|<cell|r<rsup|2>*sin\<phi\>*<with|math-font-series|bold|\<nabla\>>\<cdot\>(2*\<eta\>*<with|math-font-series|bold|D>)<rsub|r>>|<cell|=>|<cell|<left|[>*\<eta\>*h<rsub|2>*h<rsub|3>*<frac|u<rsub|1,1>|h<rsub|1>><right|]><rsub|,1>+<left|[>\<eta\>*h<rsub|3>*h<rsub|1>*<frac|u<rsub|1,2>|h<rsub|2>><right|]><rsub|,2>+<left|[>\<eta\>*h<rsub|1>*h<rsub|2>*<frac|u<rsub|1,3>|h<rsub|3>><right|]><rsub|,3>>>|<row|<cell|>|<cell|>|<cell|-2*\<eta\>*u<rsub|r>*sin\<phi\>>>|<row|<cell|>|<cell|>|<cell|-2*\<eta\>*u<rsub|\<phi\>>*cos\<phi\>>>|<row|<cell|>|<cell|>|<cell|-2*\<eta\>*(u<rsub|\<phi\>,\<phi\>>*sin\<phi\>+u<rsub|\<theta\>,\<theta\>>*)>>|<row|<cell|>|<cell|=>|<cell|\<eta\>*<left|[>(r<rsup|2>*sin\<phi\>*u<rsub|r,r>)<rsub|,r>+(sin\<phi\>*u<rsub|r,\<phi\>>)<rsub|,\<phi\>>+<left|(><frac|u<rsub|r,\<theta\>>|sin\<phi\>><right|)><rsub|,\<theta\>><right|]>>>|<row|<cell|>|<cell|>|<cell|-2*\<eta\>*(u<rsub|r>*sin\<phi\>+u<rsub|\<phi\>>*cos\<phi\>+u<rsub|\<phi\>,\<phi\>>*sin\<phi\>+u<rsub|\<theta\>,\<theta\>>*)>>>>
  </eqnarray*>

  which matches the result in Germain and Muller, chapter XIII.8.2.b.
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

      <with|par-left|<quote|1.5fn>|2<space|2spc>Viscous term
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1.5fn>|3<space|2spc>Application to spherical
      coordinates <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>
    </associate>
  </collection>
</auxiliary>