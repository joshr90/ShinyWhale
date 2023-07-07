
# ShinyWhale

------------------------------------------------------------------------

Estimate abundance and demographic parameters for North Atlantic right
whales.

### Overview

------------------------------------------------------------------------

This app is designed to provide simplified tool for estimating abundance
and demographic information for North Atlantic right whales. With the
use of this app, we aim to provide managers, NGOs, government
organisations and conservationists a tool to assess the changes in the
NARW population, without the need to understand the complete workings of
Bayesian mark-recapture techniques, or strong programming skills. This
app is designed to use data provided by the North Atlantic Right Whale
Consortium, and uses their data structure to conduct the model. The
columns required for this app are as below, please ensure that the
column headings match:

| SightingEGNo                                | SightingYear | SightingMonth | SightingDay     | Behaviors                                                            |
|---------------|---------------|---------------|---------------|---------------|
| Unique identifier for each individual whale | year(YYYY)   | month (MM)    | day of the year | list of behaviors for each sighting in capitals, separated by commas |

------------------------------------------------------------------------

**Before you start you must download [JAGS (Just Another Gibbs
Sampler)](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/).**

## Setup
ShinyWhale relies on the following packages:
- data prep: `unix`, `dplyr`,`reshape2`,`MASS`,`stringr`
- Shiny: `shiny`, `shinybusy`
- Model: `R2jags`
- Outputs: `MCMCvis`,`ggplot2`


### Data

The user supplied data needs to be converted into an observation matrix
that can be used in the Open-population Jolly-Seber framework.

#### Capture histories

Observations of individuals (i) in a given year (t) are compressed into
a single value:

1.  seen alive
2.  recovered dead, or
3.  neither seen nor recovered

The culmination of observations for each whale during the sampling
periods for the study forms the capture history. The capture history is
then combined with an additional period at the beginning of the study
period where all individuals are assigned a **3** (neither seen nor
recovered) which results in the observation matrix.

#### Data augmentation

The data set can also be adjusted to include data augmentation. This is
the inclusion of additional individuals into the observation matrix, who
are in state 3 (neither seen nor recovered). This is to account for
individuals who could exist in the population over the study period,
though have never been sighted. Individuals added in data augmentation
are not all going to enter the population, but can provide more
realistic estimates for abundance if it is believed that not every
individual has been sighted. More information on data augmentation can
be found here:

### Model

The model used in NARWPop is a Bayesian multi-event Jolly-Seber
framework [@modeling2009a; @schaub2013] fit with a
mark-recapture-recovery model [@liljestrand2019]{Lebreton, 2001 #162}.
In the model, we considered five true biological states:

1.  not yet entered the population (NE)
2.  alive within the study area (AI)
3.  alive outside the study area (AO)
4.  recovered dead (RD) and
5.  dead (D)

The focal temporal parameters that are estimated in this model are:

-   true survival (*s*)
-   recapture (*p*)
-   site fidelity (*F*)
-   dead recovery (*r*) and
-   abundance (*N*)

The state of all individuals *i* in the first occasion ($z_{i,1}$) was
set so that they were considered to not have entered the population
(*NE*):

$$
z_{i,1} = 1
$$

The state of individual *i* from the second occasion ($z_{i,2}$) until
they enter the population is:

$$
z_{i,2} = \varphi
$$

While the subsequent states are dependent on that state of the previous
time period, thus the state model is denoted as:

$$z_{i,t}|z_{i,t-1} = categorical(\Omega_{z~i,t},1...s,i,t)$$

where $z_{i,t}$ denotes the state of individual *i* at time *t*, and
$\Omega_{z_{i,t}},1....s,i,t$ denotes state membership over time, where
*s* is the number of true states.

The observation model is denoted as:

$$y_{i,t}|z_{i,t-1} = categorical(\Theta_{z~i,t},1...o,i,t)$$

where $y_{i,t}$ denotes the observation of individual *i* at time *t*, and
$\Theta_{z_{i,t}},1...o,i,t$ denotes the observational process, linking
the true states z~i,t~ to the observed states $y_{i,t}$, where *o* is the
number of observed states.

The number of individuals entering the population at time *t* ($B_{t}$) can
be calculated as:

$$
B_t = \sum_{i=1}^M (1 - z_{i,t-1})z_{i,t}
$$

while the total population size at time *t* is:

$$
N_t = \sum_{i=1}^M z_{i,t}
$$

Therefore, the full conditional probability of the model can be
expressed as:

$$
[z,s,F,r,p,\varphi|y] ∝
$$

$$
P\{z,y|s,F,r,p,\varphi\} =
$$

$$
\prod_{i=1}^M \prod_{t=2}^T \{z_{i,t}|z_{i,t-1},s_{i,t},F_{i,t},r_{i,t}, \varphi_{i,t}\} \times \{y_{i,t}|z_{i,t},p_{i,t}\}
$$

$$
\times [s][p][F][r][\varphi]
$$
