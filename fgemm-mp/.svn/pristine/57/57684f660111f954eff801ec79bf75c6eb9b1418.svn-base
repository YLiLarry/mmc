<TeXmacs|1.99.2>

<style|article>

<\body>
  <\hide-preamble>
    \;
  </hide-preamble>

  <section|Introduction>

  biblio ? checker Bernstein.

  <paragraph|Notations>

  <\itemize-dot>
    <item><math|<math-ss|M><around*|(|n|)>> polynomial arithmetic and
    <math|<math-ss|I><around*|(|n|)>> multiprecision integer arithmetic.\ 

    <item>We denote resp. by <math|<around*|(|a rem p|)>> and
    <math|<around*|(|a quo p|)>> the remainder and quotient of the Euclidean
    division of <math|a\<in\>\<bbb-Z\>> by <math|p\<in\>\<bbb-N\>> where
    <math|0\<leqslant\><around*|(|a rem b|)>\<less\>b>.

    <item>For any <math|a=n/d> with <math|d> coprime to <math|p>, we let
    <math|<around*|[|a|]><rsub|p>> be the unique representative in
    <math|<around*|{|0,\<ldots\>,p-1|}>> of <math|a> modulo <math|p>.

    <item>We call informally a pseudo-reduction of <math|a> modulo <math|p>
    the computation of <math|b> such that <math|a=b<bmod>p> and <math|b>
    ``not too big'' compared to <math|p>. In practice, we often will have
    that <math|b=\<cal-O\><around*|(|p<rsup|2>|)>>.

    <item>Say that our complexity model is bit complexity.

    <item>Should we take <math|\<beta\>> as a constant and simplify
    complexity ?

    <item>Only give costs for the typical case
    <math|p<rsub|i>\<simeq\>\<beta\>,s\<simeq\>B,r\<gg\>s> because of
    <math|\<beta\><rsup|B>\<less\>p<rsub|1>*\<cdots\>*p<rsub|s>>, which
    implies that <math|B\<less\>s> ?
  </itemize-dot>

  <paragraph|Bibliography> Bibliography on reductions and pseudo-reductions.
  Recall Barett, Montgomery results ?

  Cost of <math|a mod p> when :\ 

  <\enumerate-numeric>
    <item><math|log<around*|(|a|)>=\<Theta\><around*|(|log<around*|(|p|)>|)>>.
    **Cas classique, on en a vraiment besoin \U
    <math|\<cal-O\><around*|(|<math-ss|I><around*|(|log<around*|(|p|)>|)>|)>>
    ?**

    <item><math|log<around*|(|a|)>\<gg\>\<Theta\><around*|(|log<around*|(|p|)>|)>>
    **Servira à prouver l'algo naïf de multi-réduction \U
    <math|\<cal-O\><around*|(|log<around*|(|a|)>/log<around*|(|p|)>*<math-ss|I><around*|(|log<around*|(|p|)>|)>|)>>
    ?**

    <item><math|log<around*|(|a/p|)>\<ll\>\<Theta\><around*|(|log<around*|(|p|)>|)>>
    **Sert à montrer que la finalisation des pseudos-réductions est peu
    couteuse**

    **Polynomial analog suggests <math|\<cal-O\><around*|(|<math-ss|I><around*|(|c|)>*log<around*|(|p|)>/c|)>>
    where <math|c\<assign\>log<around*|(|a/p|)>>. Maybe we don't need to be
    so specific**
  </enumerate-numeric>

  <section|Conversions with Residue Number System>

  <subsection|Residue Number System>

  <subsection|Naive approach>

  In order to convert an integer <math|a> to a residue number system
  <math|<around|(|m<rsub|1>,\<ldots\>*m<rsub|k>|)>>, one can of course apply
  an Euclidean division of <math|a> by <math|m<rsub|i>> for
  <math|i\<in\><around|{|1,\<ldots\>,k|}>>.

  In this section, we want to reduce the integers <math|n\<in\>\<bbb-Z\>>
  modulo each of the positive integers <math|p<rsub|1>,\<ldots\>,p<rsub|s>\<in\>\<bbb-N\>>.
  Let us assume that <math|p<rsub|1>,\<ldots\>,p<rsub|s>\<less\>\<beta\>> and
  that <math|n<rsub|i>\<less\>\<beta\><rsup|B><rsup|>>. In practice,
  <math|\<beta\>> will be related to a certain number of machine words.

  <\algorithm>
    <strong|Input:> <math|n=<big|sum><rsub|j=0><rsup|B-1>c<rsub|j>*\<beta\><rsup|j>,p>

    <strong|Output:> <math|n<bmod>p>

    \;

    <strong|Algo:>

    <\indent>
      <math|c=c<rsub|B-1>>

      for <math|i=B-2\<ldots\>0> do

      <\indent>
        <math|r=c*\<beta\>><bmod>p

        <math|c=r+c<rsub|i>>
      </indent>

      <strong|return> <math|c<bmod>p>
    </indent>
  </algorithm>

  The bit complexity for computing <math|n<bmod>p> is
  <math|\<cal-O\><around*|(|B*<math-ss|I><around*|(|log \<beta\>|)>|)>>. The
  conversion to the RNS thus costs <math|\<cal-O\><around*|(|s*B*<math-ss|<samp|I>><around*|(|log
  \<beta\>|)>|)>>. **<math|\<cal-O\><around*|(|s<rsup|2>*<math-ss|<samp|I>><around*|(|log
  \<beta\>|)>|)>>**

  **Mention Barett/Montgomery optimizations ? Precomputation of floating
  number <math|\<beta\>/p> ?**\ 

  \;

  <subsection|Quasi-linear approach>

  Classic binary tree approach. Precomputation of binary tree :
  <math|\<cal-O\><around*|(|<math-ss|I><around*|(|s*log \<beta\>|)>*log
  s|)>>. [MCA, 3rd edition, Th 9.17]

  Cost in typical case : <math|\<cal-O\><around*|(|<math-ss|I><around*|(|s*log
  \<beta\>|)>*log s|)>>

  <section|Simultaneous RNS conversions>

  <subsection|Straighforward>

  Advantage of simultaneous reductions with naive algorithm : can benefit
  from (SIMD) vectorized instructions. **<math|\<cal-O\><around*|(|r*s<rsup|2>*<math-ss|<samp|I>><around*|(|log
  \<beta\>|)>|)>>**.

  Straighforward simultaneous reduction using quasi-linear approach :
  <math|\<cal-O\><around*|(|r*<math-ss|I><around*|(|s*log \<beta\>|)>*log
  s|)>>. However, do not benefit from SIMD.

  <subsection|Linear Algebra>

  <subsubsection|Linear algebra reductions>

  <paragraph|Simultaneous pseudo-reductions>In this section, we want to
  simultaneously reduce the integers <math|n<rsub|1>,\<ldots\>,n<rsub|r>\<in\>\<bbb-Z\>>
  modulo each of the positive integers <math|p<rsub|1>,\<ldots\>,p<rsub|s>\<in\>\<bbb-N\>>.

  Let us assume that <math|p<rsub|1>,\<ldots\>,p<rsub|s>\<less\>\<beta\>> and
  that <math|n<rsub|i>\<less\>\<beta\><rsup|B><rsup|>>. In practice,
  <math|\<beta\>> will be related to a certain number of machine words.

  The first thing to do is to write the expansion in base <math|\<beta\>> of

  <\equation*>
    n<rsub|i>=<big|sum><rsub|j=0><rsup|B-1>c<rsub|i,j>*\<beta\><rsup|j>
  </equation*>

  for <math|1\<leqslant\>i\<leqslant\>r>. Let's precompute the values
  <math|r<rsub|i,j>\<assign\>\<beta\><rsup|i> rem p<rsub|j>> for
  <math|0\<leqslant\>i\<less\>B> and <math|1\<leqslant\>j\<leqslant\>s>.

  Let <math|n<rsub|i,\<ell\>>\<assign\><big|sum><rsub|j=0><rsup|B-1>c<rsub|i,j>*r<rsub|j,\<ell\>>>
  then we have <math|n<rsub|i>=n<rsub|i,\<ell\>><bmod>p<rsub|\<ell\>>>. The
  value <math|n<rsub|i,\<ell\>> > is bounded by <math|B*\<beta\><rsup|2>>,
  whereas <math|n<rsub|i>> was of size <math|\<beta\><rsup|B>>. We say that
  <math|n<rsub|i,\<ell\>>> is a pseudo-reduction of <math|n<rsub|i>> modulo
  <math|p<rsub|\<ell\>>>.

  The values <math|n<rsub|i,\<ell\>>> can be computed by linear algebra :
  <math|<around*|(|n<rsub|i,\<ell\>>|)>\<in\>\<cal-M\><rsub|r\<times\>s><around*|(|k|)>>
  is the product of <math|<around*|(|c<rsub|i,j>|)>\<in\>\<cal-M\><rsub|r\<times\>B><around*|(|k|)>>
  and <math|<around*|(|r<rsub|j,\<ell\>>|)>\<in\>\<cal-M\><rsub|B\<times\>s><around*|(|k|)>>.

  <paragraph|Cost.>In the case where we want to represent our integers in the
  RNS representation, we will choose <math|p<rsub|1>,\<ldots\>,p<rsub|s>>
  such that <math|\<beta\><rsup|B>\<less\>p<rsub|1>*\<cdots\>*p<rsub|s>>,
  which implies that <math|B\<less\>s>. Then the matrix product to compute
  <math|<around*|(|n<rsub|i,\<ell\>>|)>> can be done in bit complexity
  <math|\<cal-O\><around*|(|r/s\<cdot\>s<rsup|\<omega\>>*<math-ss|I><around*|(|log<around*|(|\<beta\>|)>|)>|)>=\<cal-O\><around*|(|r*s<rsup|\<omega\>-1>*<math-ss|I><around*|(|log<around*|(|\<beta\>|)>|)>|)>>.
  The precomputation of the residues <math|<around*|(|r<rsub|i,j>|)>> costs
  <math|\<cal-O\><around*|(|s<rsup|2>*<math-ss|I><around*|(|log<around*|(|\<beta\>|)>|)>|)>>.

  \;

  <paragraph|Simultaneous reductions>Now the cost of computing the remainder
  <math|<around*|(|a rem p|)>> when <math|a\<less\>B*\<beta\><rsup|2>> and
  <math|p\<less\>\<beta\>> is <math|\<cal-O\><around*|(|<math-ss|I><around*|(|log<around*|(|\<beta\>*B|)>|)>|)>>.
  Therefore, our final step to compute our simultenous reductions costs
  <math|\<cal-O\><around*|(|r*s*<math-ss|I><around*|(|log<around*|(|\<beta\>*B|)>|)>|)>>.

  <subsubsection|Linear algebra reconstructions>

  Hypothesis for the reconstruction : <math|p<rsub|1>,\<ldots\>,p<rsub|s>>
  are pairwise coprime.\ 

  <paragraph|Simultaneous pseudo-reconstructions>Let
  <math|P=p<rsub|1>*\<cdots\>*p<rsub|s>>, <math|P<rsub|i>=P/p<rsub|i>> for
  <math|1\<leqslant\>i\<leqslant\>s>. Let
  <math|l<rsub|i>\<assign\><big|sum><rsub|j=1><rsup|s>n<rsub|i,j>P<rsub|j>*<around*|[|P<rsub|j><rsup|-1>|]><rsub|p<rsub|j>>>
  so that <math|n<rsub|i>=l<rsub|i><bmod>P> with
  <math|l<rsub|i>\<less\>P*\<beta\>>. Then <math|l<rsub|i>> are
  pseudo-reconstructions of <math|<around*|(|n<rsub|i,\<ell\>>|)>> modulo
  <math|p<rsub|1>,\<ldots\>,p<rsub|s>>.

  Once again, we perform the computation of <math|l<rsub|i>> using linear
  algebra. Let <math|P<rsub|j>*<around*|[|P<rsub|j><rsup|-1>|]><rsub|p<rsub|j>>=<big|sum><rsub|k=0><rsup|s-1>e<rsub|j,k>*\<beta\><rsup|k>>
  be the expansion in base <math|\<beta\>> of
  <math|P<rsub|\<ell\>>*<around*|[|P<rsub|\<ell\>><rsup|-1>|]><rsub|p<rsub|\<ell\>>>>.
  Put together, we have

  <\equation*>
    l<rsub|i>\<assign\><big|sum><rsub|j=1><rsup|s>n<rsub|i,j>P<rsub|j>*<around*|[|P<rsub|j><rsup|-1>|]><rsub|p<rsub|j>>=<big|sum><rsub|j=1><rsup|s>n<rsub|i,j><big|sum><rsub|k=0><rsup|s-1>e<rsub|j,k>*\<beta\><rsup|k>=<big|sum><rsub|k=0><rsup|s-1><around*|(|<big|sum><rsub|j=1><rsup|s>n<rsub|i,j>\<cdot\>e<rsub|j,k>|)>*\<beta\><rsup|k>.
  </equation*>

  Let <math|<around*|(|d<rsub|i,j>|)>\<in\>\<cal-M\><rsub|r\<times\>s>> be
  the product of the matrices <math|<around*|(|n<rsub|i,j>|)>\<in\>\<cal-M\><rsub|r\<times\>s>>
  and <math|<around*|(|e<rsub|j,k>|)>\<in\>\<cal-M\><rsub|s\<times\>s>>. Then
  <math|l<rsub|i>=<big|sum><rsub|k=0><rsup|s-1>d<rsub|i,j>*\<beta\><rsup|k>>.\ 

  Note that <math|<around*|(|d<rsub|i,j>|)>> are not the exact coefficients
  of the <math|\<beta\>>-expansion of <math|l<rsub|i>> since
  <math|d<rsub|i,j>\<leqslant\>s*\<beta\><rsup|2>>. But the correction's cost
  is <math|\<cal-O\><around*|(|s*<math-ss|I><around*|(|log \<beta\>|)>|)>>.

  **say something about <math|d<rsub|i,j>> not being the
  <math|\<beta\>>-expansion, but close**.

  <paragraph|Cost.>The matrix product to compute
  <math|<around*|(|d<rsub|i,j>|)>> can be done in bit complexity
  <math|\<cal-O\><around*|(|r/s\<cdot\>s<rsup|\<omega\>>*<math-ss|I><around*|(|log<around*|(|\<beta\>|)>|)>|)>=\<cal-O\><around*|(|r*s<rsup|\<omega\>-1>*<math-ss|I><around*|(|log<around*|(|\<beta\>|)>|)>|)>>.
  Precomputation of <math|<around*|(|e<rsub|j,k>|)>> costs
  <math|\<cal-O\><around*|(|s*<math-ss|I><around*|(|s*log \<beta\>|)>|)>>.

  <paragraph|Simultaneous reconstructions>The final step of the
  reconstruction consists in reducing <math|l<rsub|i>> modulo <math|P>. This
  step is relatively cheap since <math|l<rsub|i>> is almost reduced.

  Using the reduction when <math|log<around*|(|a/p|)>\<ll\>\<Theta\><around*|(|log<around*|(|p|)>|)>>
  in our case <math|l<rsub|i>=\<cal-O\><around*|(|s*\<beta\><rsup|s+1>|)>>
  and <math|P=\<cal-O\><around*|(|\<beta\><rsup|s>|)>>, the last reduction
  step costs ??<math|\<cal-O\><around*|(|s*log<around*|(|\<beta\>|)>*<math-ss|I><around*|(|log<around*|(|s*b|)>|)>/log<around*|(|s*b|)>|)>>**
  per <math|l<rsub|i>>. Thus a total cost of
  <math|<wide|\<cal-O\>|~><around*|(|r*s*log \<beta\>|)>>.

  \;

  <subsection|Hybrid approach ?>

  **Linear algebra up to intermediate sizes. Asymptotic complexity (binary
  tree) equivalent (not even a change of the constant).**

  <big-figure|<with|gr-mode|<tuple|edit|point>|gr-frame|<tuple|scale|1.41422cm|<tuple|0.0900389gw|0.0900355gh>>|gr-geometry|<tuple|geometry|0.50002par|0.393342par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|10|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|magnify|1.41421356236662|<graphics||<spline|<point|0.3|0.0>|<point|1.0|0.3>|<point|1.7|1.0>|<point|2.5|2.2>|<point|2.8|3.2>>|<math-at|Linear
  Algebra|<point|2.39999996930345|3.30000298417932>>|<math-at|Binary
  tree|<point|0.1|1.63970720641714>>|<math-at|hybrid|<point|2.4|1.1>>|<spline|<point|0.3|0.6>|<point|0.7|0.9>|<point|1.3|1.3>|<point|2.0|1.6>|<point|3.0|1.8>|<point|4.0|1.9>>|<spline|<point|0.3|0.0>|<point|1.0|0.3>|<point|1.7|1.0>|<point|2.2825154055729|1.4602963014932>|<point|3.2|1.7>|<point|4.0|1.8>>>>|>

  <section|Matrix Multiplication with multi-precision integer coefficients>

  <section|Implementation>

  <subsection|Reduction to word-size matrix multiplication>

  <subsection|Kronecker substitution>

  <subsubsection|From integer to <math|\<beta\>>-adic>

  <subsection|Linear storage for multi-modular matrix>

  <section|Benchmarks>

  <subsection|Conversion to and from the residue number system>

  \;

  \;

  \;
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|addrns|<tuple|1|?>>
    <associate|auto-1|<tuple|1|?>>
    <associate|auto-10|<tuple|3.2|?>>
    <associate|auto-11|<tuple|3.2.1|?>>
    <associate|auto-12|<tuple|3.2.1.1|?>>
    <associate|auto-13|<tuple|3.2.1.2|?>>
    <associate|auto-14|<tuple|3.2.1.3|?>>
    <associate|auto-15|<tuple|3.2.2|?>>
    <associate|auto-16|<tuple|3.2.2.1|?>>
    <associate|auto-17|<tuple|3.2.2.2|?>>
    <associate|auto-18|<tuple|3.2.2.3|?>>
    <associate|auto-19|<tuple|3.3|?>>
    <associate|auto-2|<tuple|1|?>>
    <associate|auto-20|<tuple|1|?>>
    <associate|auto-21|<tuple|4|?>>
    <associate|auto-22|<tuple|5|?>>
    <associate|auto-23|<tuple|5.1|?>>
    <associate|auto-24|<tuple|5.2|?>>
    <associate|auto-25|<tuple|5.2.1|?>>
    <associate|auto-26|<tuple|5.3|?>>
    <associate|auto-27|<tuple|6|?>>
    <associate|auto-28|<tuple|6.1|?>>
    <associate|auto-3|<tuple|2|?>>
    <associate|auto-4|<tuple|2|?>>
    <associate|auto-5|<tuple|2.1|?>>
    <associate|auto-6|<tuple|2.2|?>>
    <associate|auto-7|<tuple|2.3|?>>
    <associate|auto-8|<tuple|3|?>>
    <associate|auto-9|<tuple|3.1|?>>
    <associate|mulrns|<tuple|2|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|figure>
      <tuple|normal||<pageref|auto-20>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Notations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Bibliography
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Conversions
      with Residue Number System> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Residue Number System
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>Naive approach
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>Quasi-linear approach
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Simultaneous
      RNS conversions> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Straighforward
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Linear Algebra
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|2tab>|3.2.1<space|2spc>Linear algebra reductions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|4tab>|Simultaneous pseudo-reductions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Cost. <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Simultaneous reductions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.15fn>>

      <with|par-left|<quote|2tab>|3.2.2<space|2spc>Linear algebra
      reconstructions <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|4tab>|Simultaneous pseudo-reconstructions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Cost. <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Simultaneous reconstructions
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Hybrid approach ?
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Matrix
      Multiplication with multi-precision integer coefficients>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Implementation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22><vspace|0.5fn>

      <with|par-left|<quote|1tab>|5.1<space|2spc>Reduction to word-size
      matrix multiplication <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|1tab>|5.2<space|2spc>Kronecker substitution
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|2tab>|5.2.1<space|2spc>From integer to
      <with|mode|<quote|math>|\<beta\>>-adic
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25>>

      <with|par-left|<quote|1tab>|5.3<space|2spc>Linear storage for
      multi-modular matrix <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Benchmarks>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1<space|2spc>Conversion to and from the
      residue number system <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28>>
    </associate>
  </collection>
</auxiliary>