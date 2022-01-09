# Delay master stability of inertial oscillator networks

An analytical approach to determine the stability of synchronized oscillators on complex networks with delayed dynamics.

> This repository bundles work developed in the context of my bachelor's [thesis](https://github.com/reykboerner/delay-networks/blob/master/boerner_BA-thesis.pdf) and provides supplemental material for the corresponding peer-reviewed [article](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023409).

> Click [here](https://github.com/reykboerner/delay-networks/blob/master/info/220107_talk.pdf) for my **presentation slides** shown at the SIAM UKIE Annual Meeting 2022 on January 7. Slides from a longer talk are available [here](https://github.com/reykboerner/delay-networks/blob/master/info/210127_talk.pdf). To play around with the simulation slider, [read more here](#running-the-code).

## Where to start
* **What's this about?** To learn more about this research, have a look at the [plain language summary](https://github.com/reykboerner/delay-networks/blob/master/info/plain-summary.md) (0.4 pages), the peer-reviewed [paper](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.023409) (4 pages), or the [thesis](https://github.com/reykboerner/delay-networks/blob/master/boerner_BA-thesis.pdf) (40 pages).

* **Where can I find supplemental meterial for the paper?** Providing more detail than the paper, the thesis adds physical background, derivations, an in-depth treatment of the power grid application, supplemental figures, and simulation results ([read more about the thesis](#about-the-thesis)). <br/>
Some [code](#running-the-code) for DDE simulations is available; more will be added soon or can be requested (contact [me](mailto:reyk.boerner@fu-berlin.de)).

* **I want to run the code. How do I get started?** Read more [here](#running-the-code).


<br/>

## About the thesis
I conducted my bachelor's thesis in Dr. Frank Hellmann's group ["Dynamics, stability and resilience of complex hybrid infrastructure networks"](https://www.pik-potsdam.de/research/complexity-science/research/dynamics-stability-and-resilience-of-complex-hybrid-infrastructure-networks), which is part of *Research Domain 4 - Complexity Science* at the Potsdam Institute for Climate Impact Research (PIK). The thesis was formally supervised by Prof. Dr. Petra Imhof from the Institute of Theoretical Physics at Freie Universität Berlin (FU).

#### Content
- The **theoretical basics**: asymptotic stability, delay differential equations, synchronization in complex networks (chapter 2)
- A step-by-step **derivation** of the general method, including both processing and communication delay (chapter 3)
- A detailed discussion of the **application** to two future power grid models, including a comparison with simulation results (chapters 5 and 6)
- **Context** information on power grids and their study using complex networks tools (chapter 4)
- **Discussion** of the results and outlook (chapters 7 and 8).

<br/>

## Running the code

Currently, this repository contains a **slider application** which allows DDE simulations of the Decentral Smart Grid Control model on a 4-node network (see chapter 6 of the [thesis](https://github.com/reykboerner/delay-networks/blob/master/boerner_BA-thesis.pdf) for details). For a given delay, the application shows the position on the delay master stability function (left panel) and simulation results of the four nodal frequency deviations, after a random initial perturbation, as a function of time (right panel).

<p align = "center"><img src="https://github.com/reykboerner/delay-networks/blob/master/figures/dsgc_star_slider_snap.png" alt="slider-snap" width="90%"/></p>

### Prerequisites
To play around with the code, you need to have the [Julia](https://julialang.org/) programming language installed (`v1.1.0` or higher, last tested with `v1.5.3`). The following Julia packages will be required:
- `NLsolve.jl` (to obtain the synchronous fixed point)
- `DifferentialEquations.jl` (to perform the DDE simulations)
- `Plots.jl` (for plotting)
- `LaTeXStrings.jl` (for plotting)
- `Interact` (to make the slider)
- `Mux` (to make the slider)
- `DelimitedFiles` (to read data from files)

### Slider
* Run the file `dsgc_star_slider.jl` located in the `code` folder.
* If you wish, you can customize the set of slider values (delay in seconds) by modifying `tau_list` in line 12.
* Once the code is executed, go to your webbrowser and enter `localhost:8001/` in the search bar. The slider should appear.

*Note:* The precomputed datapoints of the delay master stability function (the curves in the left panel) are stored in the `.txt` files at the path `code/data`. To ensure that the Julia program can read these files, you must be in the working directory of the Julia file (i.e. `code`). To find out your current working directory, type the command `pwd()` into the REPL.

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
