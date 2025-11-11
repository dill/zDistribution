# Z-variabilty in POD of rare taxa
This repo contains all code and analyses used to examine depth variability in the probability of detecting cetaceans in eDNA samples. A living version of the manuscript in progress can be found [here](https://mmarinedna.github.io/zDistribution/). Herein, we analyse marine mammal detections from water samples collected throughout the California Current and at depths of 0-500m with the goal of answering the following three questions:

1. Does the probability of a detecting cetaceans in eDNA samples vary with sample depth?
   - H0: Probability of detection does not vary with depth.
   - H1: Probability of detection varies across depth agnostic to species or functional group.
   - H2: Probability of detection varies across depth according to species or functional group.
2. Does the probability of a detecting cetaceans in eDNA samples change with the number of technical replicates?
   * NOTE should we wrap dilution into this question as well?
   - H0: Detection does not vary with # technical replicates.
   - H1: Detection varies with # technical replicated agnostic to species/fuctional group or depth.
   - H2: Detection varies with # technical replicates according to species/functional group, depth, or a combination of the two.
3. Does depth distribution of detections vary across xy spatial distribution?
   - H0: Depth distribution of detections does not vary across xy spatial distribution.
   - H1: Depth distribution of detection does vary across xy spatial distribution agnostic to oceanography (e.g. upwelling).
   - H2: Depth distribution of detection does vary across xy spatial distribution according to oceanography (e.g. upwelling).

# The Plan

All raw data are in "Data", and intermediate data products are in "ProcessedData." Dive data come from [here](https://apps.dtic.mil/sti/tr/pdf/ADA560975.pdf) and [here](https://www.nepa.navy.mil/Portals/20/Documents/aftteis4/Dive%20Profile%20and%20Group%20Size_TR_2017_05_22.pdf).

# Analysis 1: POD by depth, collapse data across X and Y

- Model 1: POD ~ z
- Model 2: POD ~ z with intercept by species
- Model 3: POD ~ z with smooth by species
- Model 4: POD ~ z with smooth by family
- Model 5: POD ~ z with smooth by prey category (invert, fish, squid)
- Run model diagnostics, select best model, interpret results within ecological context, develop recommendations for eDNA monitoring of marine mammals.

# Analysis 2: POD by primer, freeze/thaw, tech rep. and vary with depth?

- Using GAM splines from Analysis 1, retest varying number of replicates with GAMs
- Compare to alternative Bayesian occupancy model: [Brice's replication model](https://github.com/BriceSemmens/eDNA_patch) without assuming species' presence and adding depth as a covariate (see "Taking Brice's approach" below).
- Run model diagnostics, select best model, interpret results within ecological context, develop recommendations for eDNA monitoring of marine mammals.

# Analysis 3: depth variability across xy or oceanography

- Model 1: POD ~ xy + xy x species + z + z x species + xy x species x z
- Run model diagnostics, report significant covariates, interpret results within ecological context, develop recommendations for eDNA monitoring of marine mammals.

# Taking Brice's approach

We don't have observations of the presence or absence of marine mammals at each site, but we can still treat site occupancy as a latent variable (per species). 

We have: multiple sites, each of which contain biological samples at multiple depths (not replicates), then multiple primers x tech reps. Each has an associated volume and dilution. 

* need to add a dilution coefficient, or does this happen at the tech rep level?

### Hierarchical Model Structure

1. **Site-Level Occurrence**
   The presence of a marine mammal species at site $s$ is modeled as a Bernoulli random variable:

   $Z_s \sim \text{Bernoulli}(\psi)$

   where $Z_s$ is the site-level occurrence indicator, and $\psi$ is the overall occurrence probability, drawn from a Beta prior:

   $\psi \sim \text{Beta}(1,1)$

   --> note that in our case, particularly since we are working at the species level, psi might come from some field over X, Y

2. **Capture within a depth x primer x tech rep reaction **
   The logit-linear model for capture probability at depth $d$ is:

   $\text{logit}(p_{\text{capture},d}) = \beta_0 + \beta_{\text{vol}} \cdot X_{\text{vol},d} + \beta_{\text{depth}} \cdot X_{\text{depth},d} + \gamma_{s} + \delta_{\text{primer}[b]}$

   Where:
   - $p_{\text{capture},d}$ is the capture probability for depth $d$
   - $\beta_0$ is the intercept (site-level capture probability hyperparameter)
   - $\beta_{\text{vol}}$ is the volume coefficient
   - $\beta_{\text{depth}}$ is the depth coefficient
   - $X_{\text{vol},d}$ is the centered water volume
   - $X_{\text{depth},d}$ is the centered sampling depth
   - $\gamma_{s}$ is the site-specific random effect # EKJ note is this what we want?
   - $\delta_{\text{primer}[b]}$ is the primer-specific fixed effect, with method effects constrained such that:
   $\delta_{\text{MFU}} = 0$ and 
   $\delta_{\text{MV1}} \sim \text{Normal}(0, 1.7)$ # EKJ note may need to change these
   $\delta_{\text{DLP}} \sim \text{Normal}(0, 1.7)$
   
   The depth capture for a given site, volume, primer  is then modeled as:

   $Y_{\text{capture},d} \sim \text{Bernoulli}(Z_{s} \cdot p_{\text{capture},d})$

   Need to rewrite this bc it's only one Bernoulli trial in our case

   Conditional on capture in that volume x primer reaction, technical replicates are modeled as:

   $Y_{\text{detect},i} \sim \text{Bernoulli}(p_{\text{detect}} \cdot Y_{\text{capture},b[i]})$

   where $p_{\text{detect}}$ is the detection probability, drawn from a Beta prior:

   $p_{\text{detect}} \sim \text{Beta}(1,1)$
