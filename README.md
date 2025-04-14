# CRModel

## `Julia` routine to calculate Locally Optimal Restricted Designs for Continuation-Ratio Model. 

Here, we have a collection of `Julia` scripts to calculate locally optmal designs for the _**continuation-ratio model**_.
The model is a four-parameter model that might be used to model _**multinomial probabilites**_ in the _**dose-efficacy-toxiciity**_ studies.

$$\boldsymbol{\eta}(x; \boldsymbol{\theta}) = \left\\{
\begin{aligned}
\pi_0(x; \boldsymbol{\theta}) &= \frac{1}{1+\exp{(\theta_1+ + \theta_2 x)}}\frac{1}{1+\exp{(\theta_3+ + \theta_3 x)}} \\
\pi_1(x; \boldsymbol{\theta}) &= \frac{\exp{(\theta_1+ + \theta_2 x)}}{1+\exp{(\theta_1+ + \theta_2 x)}}\frac{1}{1+\exp{(\theta_3+ + \theta_3 x)}} \\
\pi_2(x; \boldsymbol{\theta}) &= \frac{\exp{(\theta_3+ + \theta_4 x)}}{1+\exp{(\theta_3+ + \theta_3 x)}}.
\end{aligned}\right.$$

Here, $x$ is a teratment dose, and 

- $\pi_0(x; \theta)$ is a _**probability of no response**_ (neither efficacy no toxicity).
- $\pi_1(x; \theta)$ is a _**probability of efficacy without toxicity**_.
- $\pi_2(x; \theta)$ is a _**probability of toxicity**_.

One might be interested in finding _**optimal experimental design**_ for such a model:

$$
\xi^\* = \begin{pmatrix}
x^\*_1 & \ldots & x^\*_K \\
\rho^\*_1 & \ldots & \rho^\*_K \\
\end{pmatrix},
$$

where $x^\*_1, \ldots , x^\*_K$ are _**optimal doses**_ and $\rho^\*_1, \ldots , \rho^\*_K$ are subjects' _**optimal allocation proportions**_.

Given study objectives, one can set up an optimal design problem, where $\left(x^\*_1, \ldots , x^\*_K, \rho^\*_1, \ldots , \rho^\*_K\right)$ represents 
an _**optimal solutions**_ of the optimization problem.

For instance, if one wants to esy´timate model parameters as precisely as possible, then it a _**D-optimal deisgn problem**_.

If the objective is estimating an _**optimal biological dose**_ (OBD), then it is a _**c-optimal design problem**_.

In this project, we use a _**Particle Swarm Optimization**_ (PSO) to find optimal design (_**D-**_ and _**c-optimal**_).

Scripts contain neccessary functionality needed to define a model, _**Fisher information matrix**_, optimality criteria, objective functions used in the 
optimization process, PSO algorithm for both finding designs and maximum likelihood estimations.

- Script `98-simulation-paper1.jl` has to be run to reproduce the results related to the CR model from the pablication (**"Name to be here"**).



