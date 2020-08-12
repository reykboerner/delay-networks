# Delay master stability of inertial oscillator networks

An analytical approach to determine the stability of synchronized oscillators on complex networks with delayed dynamics.
* **Necessary and sufficient conditions** for asymptotic stability as a function of delay and network parameters
* Generalization of the master stability formalism to inertial oscillator networks with **arbitrary processing-type delay**
* Application to future **power grid** models.

This repository bundles work developed in the context of my bachelor's [thesis](https://github.com/reykboerner/delay-networks/blob/master/boerner_BA-thesis.pdf) and provides supplemental material for the corresponding [article](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023409) published in *Physical Review Research*.

## Where to start
* **What's this about?** To learn more about this research, have a look at the [plain language summary](https://github.com/reykboerner/delay-networks/blob/master/info/plain-summary.md) (0.4 pages), the peer-reviewed [paper](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023409) (4 pages), or the [thesis](https://github.com/reykboerner/delay-networks/blob/master/boerner_BA-thesis.pdf) (40 pages; see below for [more info on the thesis](#about-the-thesis)).

* **I've read the paper/thesis. Where can I find additional material?** Providing more detail than the paper, the thesis adds physical background, derivations, an in-depth treatment of the power grid application, and simulation results ([read more about the thesis](#about-the-thesis)). Additionally, this repo includes [supplemental figures](https://github.com/reykboerner/delay-networks/blob/master/figures) as well as an [implementation](#code) in the Julia language (dMSF calculations and DDE simulations).

* **I want to try the simulations. How do I get started?** Read more [here](#code).

<br/><br/>


## About the thesis
I conducted my bachelor's thesis in Dr. Frank Hellmann's group ["Dynamics, stability and resilience of complex hybrid infrastructure networks"](https://www.pik-potsdam.de/research/complexity-science/research/dynamics-stability-and-resilience-of-complex-hybrid-infrastructure-networks), which is part of *Research Domain 4 - Complexity Science* at the Potsdam Institute for Climate Impact Research (PIK). The thesis was formally supervised by Prof. Dr. Petra Imhof from the Institute of Theoretical Physics at Freie Universität Berlin (FU).

**The thesis includes**
- The **theoretical basics**: asymptotic stability, delay differential equations, synchronization in complex networks (chapter 2)
- A step-by-step **derivation** of the general method, including both processing and communication delay (chapter 3)
- A detailed discussion of the **application** to two future power grid models, including a comparison with simulation results (chapters 5 and 6)
- **Context** information on power grids and their study using complex networks tools (chapter 4)
- **Discussion** of the results and outlook (chapters 7 and 8).


## Further material
### Code
> Coming soon...
### Supplemental figures
> Coming soon...

<br/><br/>

### How to cite
Preferably, please cite the paper as follows when relating to this research:
> R. Börner, P. Schultz, B. Ünzelmann, D. Wang, F. Hellmann, J. Kurths, *Delay master stability of inertial oscillator networks*, Phys. Rev. Research **2**, 023409 (2020).

```
@article{boerner2020delay,
    title = {Delay master stability of inertial oscillator networks},
    author = {B\"orner, Reyk and Schultz, Paul and \"Unzelmann, Benjamin and Wang, Deli and Hellmann, Frank and Kurths, J\"urgen},
    journal = {Phys. Rev. Research},
    volume = {2},
    issue = {2},
    pages = {023409},
    numpages = {8},
    year = {2020},
    month = {Jun},
    publisher = {American Physical Society},
    doi = {10.1103/PhysRevResearch.2.023409},
    url = {https://link.aps.org/doi/10.1103/PhysRevResearch.2.023409}
}
```

To specifically cite the thesis:
> R. Börner, *Master Stability of Inertial Oscillator Networks with Delay*, Bachelor's thesis, Freie Universität Berlin (2019).

```
@mastersthesis{boerner2019thesis,
    title = {Master Stability of Inertial Oscillators with Delay -- An analytical approach applied to renewable power grids},
    author = {Reyk Börner},
    year = {2019},
    type = {Bachelor's thesis},
    school = {Freie Universit\"at Berlin},
    url = {https://github.com/reykboerner/delay-networks}
}
```

<br/><br/>

# Collaborators
#### Frank Hellmann, Jürgen Kurths, Anton Plietzsch, Paul Schultz, Benjamin Ünzelmann, Deli Wang

This work is part of the [CoNDyNet 2](condynet.de) project, sponsored by the German Federal Ministry of Education and Research (BMBF).

<p align = "center"><img src="https://github.com/reykboerner/delay-networks/blob/master/info/logo-banner.png" alt="logo-banner" width="70%"/></p>
