SeqSGPV Package
================

# Purpose

The SeqSGPV package is used to design a study with sequential monitoring
of scientifically meaningful hypotheses using the second generation
p-value (SGPV).

It supports the paper [Sequential monitoring using the Second Generation
P-Value with Type I error controlled by monitoring
frequency](https://arxiv.org/pdf/2204.10678.pdf) which advances how to:

1.  Specify scientifically meaningful hypotheses using a constrained
    Region of Equilvance (Freedman et al., 1984). The constrained set of
    hypotheses are called Pre-Specified Regions Indicating Scientific
    Merit (PRISM).
2.  Sequentially monitor the SGPV (SeqSGPV) for scientifically
    meaningful hypotheses.
3.  Control error rates uses monitoring frequency mechanisms, including
    an affirmation step.

# Why PRISM

Monitoring until establishing statistical significance does not provide
information on whether an effect is scientifically meaningful.
Establishing scientific relevance requires pre-specifying which effects
are scientifically meaningful and an end-study inference that evaluates
these effects.

For intervention studies, Freedman et al. (1984) divides effects as
being: (1) universally acceptable for adopting an intervention, (2)
universally unacceptable, or (3) debatably acceptable. The latter set of
effects form the ROE. Strategies for specifying ROE include setting the
point null as a ROE boundary \[Hobbs & Carlin (2008); Section 2.2\],
setting the ROE away from the point null \[Freedman et al. (1984);
Figure 1\], and surrounding the point null (Kruschke, 2013). Kruschke
(2013) calls the latter strategy a Region of Practical Equivalence
(ROPE); see footnote[^1] for clarification between the ROPE and ROE.

In a 2-sided study, the PRISM includes a ROPE and Region of Meaningful
effects (ROME). In a 1-sided study, the PRISM includes a Region of Worse
or Practically Equivalent effects (ROWPE) and ROME.

<figure>
<img src="images/AMwithSGPV_Figsv09_1.jpg"
alt="Pre-Specified Regions Indicating Scientific Merit (PRISM) for one- and two-sided hypotheses. The PRISM always includes an indifference zone that surrounds the point null hypothesis (i.e. ROPE/ROWPE)." />
<figcaption aria-hidden="true">Pre-Specified Regions Indicating
Scientific Merit (PRISM) for one- and two-sided hypotheses. The PRISM
always includes an indifference zone that surrounds the point null
hypothesis (i.e. ROPE/ROWPE).</figcaption>
</figure>

In the context of interval monitoring, error rates and sample size are
impacted by ROE specifation. The PRISM is a constrained version of the
ROE in which the ROE boundary is set away from the point null. This
specification induces two adjacent regions to the ROE: the ROPE and a
Region of Meaningful Effects (ROME).

Compared to ROPE monitoring, PRISM monitoring also reduces the risk of
type I error yet resolves the issue of indefinite monitoring at ROPE
boundaries.

Compared to null-bound ROE monitoring, PRISM monitoring reduces the risk
of Type I error for the same monitoring frequency and allows for earlier
monitoring and yields smaller average sample size to achieve the same
Type I error.

# Why SGPV

The SGPV is an evidence-based metric that measures the overlap between
an inferential interval and scientifically meaningful hypotheses.
Described as ‘method agnostic’ (Stewart & Blume, 2019), the SGPV may be
calculated for any inferential interval (ex: bayesian, frequentist,
likelihood, etc.). For an interval hypothesis, $H$, the SGPV is denoted
as p$_H$ and calculated as:

The SGPV (Blume et al., 2018) is an evidence-based metric that
quantifies the overlap between an interval $I$ and the set of effects in
a composite hypothesis *H*. The interval includes \[*a*, *b*\] where *a*
and *b* are real numbers such that *a*\<*b*, and the length of the
interval is *b-a* and denoted $I$. The overlap between the interval and
the set *H* is $\left| I \cap H \right|$. The SGPV is then calculated as

$$
p_{H} = \frac{| I \cap H |}{|I|} \times \max \left\{\frac{| I |}{2 | H |}, 1 \right\}.
$$

The adjustment, $max\left\{\frac{|I|}{2|H|},1\right\}$, is a small
sample size correction – setting $p_H$ to half of the overlap when the
inferential interval overwhelms $H$ by at least twice the length.

In [Sequential monitoring using the Second Generation P-Value with Type
I error controlled by monitoring
frequency](https://arxiv.org/pdf/2204.10678.pdf), foundational
likelihood, frequentist, and Bayesian metrics are compared and
contrasted in terms of minimal assumptions, handling of composite
hypotheses, and conclusions that can be drawn. The SGPV makes no further
assumptions beyond those inherited by the inferential interval. It does
not require a likelihood, prior, study design, or error rates.

# Why change monitoring frequency

Controlling the design-based Type I error is generally considered an
important metric for reducing the risk of false discoveries. In SeqSGPV,
error rates can be controlled through PRISM specification and/or
monitoring frequency.[^2]

Since scientific relevance (i.e., PRISM) is considered fixed, monitoring
frequency is a targettable means for controlling error rates. These
include a wait time until evaluating stopping rules, the frequency of
evaluations, a maximum sample size, and an affirmation rule. The
affirmation rule is used in dose-escalation trials once a number of
patients have consecutively been enrolled at the recommended maximum
tolerable dose. We use it here to further control error rates.

**Synergy between 1-sided PRISM and monitoring frequency**: On their
own, both the PRISM and monitoring frequency help reduce the risk of
Type I error. When used together, the 1-sided PRISM and monitoring
frequency can dramatically reduce the average sample size to achieve a
Type I error. When outcomes are delayed, the risk of reversing a
decision on the null hypothesis decreases when monitoring a 1-sided
PRISM more so than under a 1-sided null-bound ROE. Additional
strategies, such as posterior predictive probabilities could be
considered to further inform decisions under delayed outcomes.

**Comment on monitoring confidence intervals**: When using confidence
intervals, the investigator should determine how to address issues of
bias and coverage which is common to sequential monitoring. This aspect
is beyond the scope of the paper [Sequential monitoring using the Second
Generation P-Value with Type I error controlled by monitoring
frequency](https://arxiv.org/pdf/2204.10678.pdf).

# Package overview

Trial designs are evaluated through simulation via the SeqSGPV function
until achieving desirable operating characteristics such as error rates,
sample size, bias, and coverage. The SeqSGPV function can be used to
assess departures from modelling assumptions.

Outcomes may be generated from any r\[dist\] distribution, a
user-supplied data generation function, or pre-existing data. Study
designs of bernoulli and normally distributed outcomes have been more
extensively evaluated and extra care should be provided when designing a
study with outcomes of other distribution families.

The user provides a function for obtaining interval of interest. Some
functions have been built for common interval estimations: binomial
credible and confidence intervals using binom::binom.confint, wald
confidence intervals using lm function for normal outcomes, and wald
confidence intervals using glm function with binomial link for bernoulli
outcomes.

Depending on computing environment, simulations may be time consuming to
obtain many (10s of thousands) replicates and more so for bernoulli
outcomes. The user may consider starting with a small number of
replicates (200 - 1000) to get a sense of design operating
characteristics. Sample size estimates of a single look trial may also
inform design parameters.

# Study design examples

Study designs and interpretations of a single trial are provided below
for 1-2 arm trials with bernoulli or normally distributed outcomes.

- [one arm trial, bernoulli
  outcomes](examples/one-arm-bernoulli/README.md)
- [one arm trial, bernoulli outcomes, effect generated from prior
  distribution](one-arm-bernoulli-prior-generated-effect/README.md)
- [two arm trial, continous
  outcomes](examples/two-arm-trial-continuous/two-arm-continuous.md)
- [two arm trial, bernoulli
  outcomes](examples/two-arm-bernoulli/README.md)

# References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-Blume:SGPV" class="csl-entry">

Blume, J. D., DAgostino McGowan, L., Dupont, W. D., & Greevy Jr., R. A.
(2018). <span class="nocase">Second-generation p-values: Improved rigor,
reproducibility, $\backslash$& transparency in statistical
analyses</span>. *PloS One*, *13*(3), e0188299 EP —–.

</div>

<div id="ref-Freedman:1984wz" class="csl-entry">

Freedman, L. S., Lowe, D., & Macaskill, P. (1984).
<span class="nocase">Stopping rules for clinical trials incorporating
clinical opinion.</span> *Biometrics*, *40*(3), 575–586.

</div>

<div id="ref-Hobbs:2008ce" class="csl-entry">

Hobbs, B. P., & Carlin, B. P. (2008). <span class="nocase">Practical
Bayesian design and analysis for drug and device clinical trials.</span>
*Journal of Biopharmaceutical Statistics*, *18*(1), 54–80.

</div>

<div id="ref-jennison1989interim" class="csl-entry">

Jennison, C., & Turnbull, B. W. (1989). <span class="nocase">Interim
analyses: the repeated confidence interval approach</span>. *Journal of
the Royal Statistical Society: Series B (Methodological)*, *51*(3),
305–334.

</div>

<div id="ref-Kruschke:2013jy" class="csl-entry">

Kruschke, J. K. (2013). <span class="nocase">Bayesian estimation
supersedes the t test.</span> *Journal of Experimental Psychology.
General*, *142*(2), 573–603.

</div>

<div id="ref-stewart2019second" class="csl-entry">

Stewart, T. G., & Blume, J. (2019).
<span class="nocase">Second-generation p-values, shrinkage, and
regularized models</span>. *Frontiers in Ecology and Evolution*, *7*,
486.

</div>

</div>

[^1]: ROPE as effects practically equivalent to the point null whereas
    theFreedman refers to the ROE as effects where the scientific
    communitymay disagree on whether an intervention has enough value to
    beimplemented. Hence, ROE is a broader term.

[^2]: Jennison & Turnbull (1989) provide error rate control for
    intervals that adjust for frequency properties using group
    sequential methods. For intervals which do not adjust for frequency
    properties, Type I error can be controlled through a combination of
    the PRISM’s ROPE and monitoring frequency mechanisms of wait time
    (W), frequency of evaluation (S for steps between evaluations),
    maximum sample size (N), and affirmation steps (A).
