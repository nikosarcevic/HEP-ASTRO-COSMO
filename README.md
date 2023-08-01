# HEP-ASTRO-COSMO

**This is a community effort to collect all HEP/ASTRO/COSMO open source packages/libraries/tools in one place. Everyone welcome to contribute!** :shipit:

**Package/library/tool descriptions copied from respective pages.**


***

## Table of contents
- [HEP](#hep)
  * [Tools](#tools)
  * [Event Generators](#event-generators)
  * [Event Analysers](#event-analysers)
  * [Global Fitters](#global-fitters)
  * [Spectrum Generators](#spectrum-generators)
  * [Model Building](#model-building)
  * [Effective Field Theories](#effective-field-theories)
  * [Direct Detection](#direct-detection)
  * [Feyn Family](#feyn-family)
  * [Statistics](#statistics)
  * [Neutrino Code](#neutrino-code)
- [COSMO](#cosmo)
  * [General Cosmology](#general-cosmology-resources)
  * [Inflationary Cosmology](#inflationary-cosmology)
  * [Big Bang Nucleosynthesis](#big-bang-nucleosynthesis)
  * [Einstein-Boltzmann Numerical Solvers](#einstein-boltzmann-numerical-solvers)
  * [Cosmological Parameter Estimation](#cosmological-parameter-estimation-and-statistical-analysis)
  * [Correlation function and higher order statistics](#correlation-function-and-higher-order-statistics)
  * [Simulations](#simulations)
  * [Extended Cosmologies](#extended-cosmologies)
  * [Perturbation Theory](#perturbation-theory)
- [ASTRO](#astro)
  * [Black hole perturbation theory](#black-hole-perturbation-theory)
  * [Cosmic Rays](#cosmic-rays)
  * [Gamma-Ray Astronomy](#gamma-ray-astronomy)
  * [N-body Simulations](#n-body-simulation)
  * [MHD](#mhd)
  * [Numerical Relativity](#numerical-relativity)
  * [Primordial Black Holes](#primordial-black-holes)
  * [Stellar Modelling](#stellar-modelling)
  * [Modelling of Active Galactic Nuclei](#Modelling-of-Active-Galactic-Nuclei)
  * [Multi-Messenger Analysis](#multi-messenger-analysis)
  * [Population Synthesis](#population-synthesis)
  * [Neutron Stars Equations of State](#Neutron-Stars-Equations-of-State)
  * [Pular Timing](#pulsar-timing)
- [Related List of Tools](#related-list-of-tools)
  * [Neutrino Code](#neutrino-code)

***
# HEP 

## Tools

### Root

ROOT is a framework for data processing, born at CERN, at the heart of the research on high-energy physics.

https://root.cern.ch/

### Geant4

Geant4 is a toolkit for the simulation of the passage of particles through matter. Its areas of application include high energy, 
nuclear and accelerator physics, as well as studies in medical and space science.

https://geant4.web.cern.ch/node/1

***

## Event Generators

### Pythia8

PYTHIA is a program for the generation of high-energy physics events, i.e. for the description of collisions at high energies between elementary particles such as e+, e-, p and pbar in various combinations. It contains theory and models for a number of physics aspects, including hard and soft interactions, parton distributions, initial- and final-state parton showers, multiparton interactions, fragmentation and decay. It is largely based on original research, but also borrows many formulae and other knowledge from the literature.

http://home.thep.lu.se/~torbjorn/Pythia.html

### hepmc

The HepMC package is an object oriented event record written in C++ for High Energy Physics Monte Carlo Generators. Many extensions from HEPEVT, the Fortran HEP standard, are supported: the number of entries is unlimited, spin density matrices can be stored with each vertex, flow patterns (such as color) can be stored and traced, integers representing random number generator states can be stored, and an arbitrary number of event weights can be included. Particles and vertices are kept separate in a graph structure, physically similar to a physics event. The added information supports the modularisation of event generators. The package has been kept as simple as possible with minimal internal/external dependencies. Event information is accessed by means of iterators supplied with the package.

https://ep-dep-sft.web.cern.ch/project/hepmc

### herwig

Herwig is a multi-purpose particle physics event generator. It is built based on the experience gained with both the HERWIG 6 and Herwig++ 2 event generators. Continuing the Herwig++ 2 development, Herwig 7.0 (Herwig++ 3.0) replaces any prior HERWIG or Herwig++ versions. Herwig provides significantly improved and extended physics capabilities when compared to both its predecessors, HERWIG 6 and Herwig++ 2, while keeping the key model motivations such as coherent parton showers (including angular-ordered and dipole-based evolution), the cluster hadronization model, an eikonal multiple-interaction model, highly flexible BSM capabilities and improved perturbative input using next-to-leading order QCD.

https://herwig.hepforge.org/

### sherpa

Sherpa is a Monte Carlo event generator for the Simulation of High-Energy Reactions of PArticles in lepton-lepton, lepton-photon, photon-photon, lepton-hadron and hadron-hadron collisions. Simulation programs - also dubbed event generators - like Sherpa are indispensable work horses for current particle physics phenomenology and are (at) the interface between theory and experiment.

https://sherpa.hepforge.org/trac/wiki

### powheg

The POWHEG BOX is a general computer framework for implementing NLO calculations in shower Monte Carlo programs according to the POWHEG method. It is also a library, where previously included processes are made available to the users. It can be interfaced with all modern shower Monte Carlo programs that support the Les Houches Interface for User Generated Processes.

http://powhegbox.mib.infn.it/

### MadGraph

MadGraph5_aMC@NLO is a framework that aims at providing all the elements necessary for SM and BSM phenomenology, such as the computations of cross sections, the generation of hard events and their matching with event generators, and the use of a variety of tools relevant to event manipulation and analysis. Processes can be simulated to LO accuracy for any user-defined Lagrangian, an the NLO accuracy in the case of models that support this kind of calculations -- prominent among these are QCD and EW corrections to SM processes. Matrix elements at the tree- and one-loop-level can also be obtained.

https://launchpad.net/mg5amcnlo

### EVTGen

EvtGen is a Monte Carlo event generator that simulates the decays of heavy flavour particles, primarily B and D mesons. It contains a range of decay models for intermediate and final states containing scalar, vector and tensor mesons or resonances, as well as leptons, photons and baryons. Decay amplitudes are used to generate each branch of a given full decay tree, taking into account angular and time-dependent correlations which allows for the simulation of CP-violating processes. Originally written by Anders Ryd and David Lange, this package is used by many particle physics experiments worldwide, including ATLAS, BaBar, Belle(-II), BES III, CDF, CLEO(-c), CMS, D0, and LHCb. The maintenance and development of the package is now performed by the particle physics group at the University of Warwick (in particular by John Back, Michal Kreps, and Thomas Latham).

https://evtgen.hepforge.org/

### Photos

Photos is a Monte Carlo for bremsstrahlung in the decay of particles and resonances is available with an interface to the C++ HepMC event record. It is usually used in conjuction with EvtGen. Authors: Davidson, T. Przedzinski, Z. Was. 

http://photospp.web.cern.ch/photospp/

### Alpaca

Alpaca is a Fortran based Monte Carlo event generator for ALP production in coherent proton-nucleus (and electron-nucleus) collisions. Arbitrary user-defined histograms and cuts may be made, as well as unweighted events in the HEPEVT, HEPMC and LHE formats.

https://alpaca.hepforge.org/

### WHIZARD 

Monte Carlo Event Generator for Tevatron, LHC, ILC, CLIC, CEPC, FCC-ee, FCC-hh, SppC, the muon collider and other High Energy Physics Experiments. WHIZARD is a program system designed for the efficient calculation of multi-particle scattering cross sections and simulated event samples. Tree-level matrix elements are generated automatically for arbitrary partonic processes by using the Optimized Matrix Element Generator O'Mega. The program is able to calculate numerically stable signal and background cross sections and generate unweighted event samples with reasonable efficiency for processes with up to eight final-state particles; more particles are possible. For more particles, there is the option to generate processes as decay cascades including complete spin correlations. Different options for QCD parton showers are available. Polarization is treated exactly for both the initial and final states. Final-state quark or lepton flavors can be summed over automatically where needed. For hadron collider physics, an interface to the standard LHAPDF is provided. For Linear Collider physics, beamstrahlung (CIRCE) and ISR spectra are included for electrons and photons. WHIZARD supports the Standard Model and a huge number of BSM models. There are also interfaces to FeynRules and SARAH.

https://whizard.hepforge.org/

### OpenLoops 2

The OpenLoops 2 program is a fully automated implementation of the Open Loops algorithm combined with on-the-fly reduction methods, which allows for the fast and stable numerical evaluation of tree and one-loop matrix elements for any Standard Model process at NLO QCD and NLO EW. 

https://openloops.hepforge.org/

### Caravel

The Caravel C++ framework provides an implementation of many algorithms necessary to perform computations of multi-loop scattering amplitudes within the multi-loop numerical unitarity method. It is based on the (generalized) unitarity approach, which was first developed for the analytic computation of one-loop amplitudes and later adapted for numeric calculations. In a nutshell, in this framework the amplitude is computed starting from a parametrization of its integrand. The corresponding free parameters are numerically computed at each phase-space point by constructing systems of linear equations in which the parameters are the unknowns and the numerical entries are associated to products of tree-level amplitudes. With a suitable choice of integrand parametrization, this directly gives a decomposition of the amplitude in terms of master integrals. Finally, after inserting the value of the integrals at the required phase-space point we obtain the value of the amplitude.

The current release of Caravel includes a module for computing products of tree-level amplitudes in several theories through off-shell recursion relations, and tools that allow the efficient construction and solution of the systems of linear equations that determine the integrand. Whilst these components work for generic multi-loop amplitudes, other components such as the construction of the parametrization are required as input.

https://gitlab.com/caravel-public/caravel

### MadFlow

MadFlow is a framework for Monte Carlo simulation of particle physics processes designed to take full advantage of hardware accelerators. Processes can be generated using [MadGraph5_aMC@NLO](#Madgraph) and are then output in vectorized (or tensorized) form by the madflow-provided plugin.

The vectorized output is compiled using the TensorFlow library (hence, tensorized) and then integrated using the [VegasFlow](https://vegasflow.readthedocs.io) library. The PDF interpolation is provided by [PDFflow](#PDFflow). All tools are capable of running hardware with different hardware acceleration capabilities, such as multi-threading CPU, single-GPU and multi-GPU setups.

https://github.com/N3PDF/madflow

***

## Event Analysers

### RIVET

The Rivet toolkit (Robust Independent Validation of Experiment and Theory) is a system for validation of Monte Carlo event generators. It provides a large (and ever growing) set of experimental analyses useful for MC generator development, validation, and tuning, as well as a convenient infrastructure for adding your own analyses. Rivet is the most widespread way by which analysis code from the LHC and other high-energy collider experiments is preserved for comparison to and development of future theory models. It is used by phenomenologists, MC generator developers, and experimentalists on the LHC and other facilities.

https://rivet.hepforge.org/

### SModelS

SModelS is based on a general procedure to decompose Beyond the Standard Model (BSM) collider signatures presenting a Z2 symmetry into Simplified Model Spectrum (SMS) topologies. Our method provides a way to cast BSM predictions for the LHC in a model independent framework, which can be directly confronted with the relevant experimental constraints.

https://smodels.github.io/

### CheckMATE

CheckMATE (Check Models At Terascale Energies) is a program package which accepts simulated event files in many formats for any given model. The program then determines whether the model is excluded or not at 95% C.L. by comparing to many recent experimental analyses. Furthermore the program can calculate confidence limits and provide detailed information about signal regions of interest. It is simple to use and the program structure allows for easy extensions to upcoming LHC results in the future.

https://checkmate.hepforge.org/

### Delphes

Delphes is a C++ framework, performing a fast multipurpose detector response simulation. The simulation includes a tracking system, embedded into a magnetic field, calorimeters and a muon system. The framework is interfaced to standard file formats (e.g. Les Houches Event File or HepMC) and outputs observables such as isolated leptons, missing transverse energy and collection of jets which can be used for dedicated analyses. The simulation of the detector response takes into account the effect of magnetic field, the granularity of the calorimeters and sub-detector resolutions. Visualisation of the final state particles is also built-in using the corresponding ROOT library.

https://cp3.irmp.ucl.ac.be/projects/delphes

### MadAnalysis

MadAnalysis 5 is a framework for phenomenological investigations at particle colliders. Based on a C++ kernel, this program allows to efficiently perform, in a straightforward and user-friendly fashion, sophisticated physics analyses of event files such as those generated by a large class of Monte Carlo event generators.

MadAnalysis 5 can also be used for the recasting of existing LHC analyses. These features are documented on the MA5 PAD (public analysis database), together with instructions to implement new analyses ([see this link](http://madanalysis.irmp.ucl.ac.be/wiki/PublicAnalysisDatabase)).

https://launchpad.net/madanalysis5

***

## Global Fitters

### HEPfit

A Code for the Combination of Indirect and Direct Constraints on High Energy Physics Models.

𝙷𝙴𝙿𝚏𝚒𝚝 is a flexible open-source tool which, given the Standard Model or any of its extensions, allows to:  
- fit the model parameters to a given set of experimental observables  
- obtain predictions for observables  

𝙷𝙴𝙿𝚏𝚒𝚝 can be used either in Monte Carlo mode, to perform a Bayesian Markov Chain Monte Carlo analysis of a given model, or as a library, to obtain predictions of observables for a given point in the parameter space of the model, allowing 𝙷𝙴𝙿𝚏𝚒𝚝 to be used in any statistical framework. In the present version, around a thousand observables have been implemented in the Standard Model and in several new physics scenarios. In this paper, we describe the general structure of the code as well as models and observables implemented in the current release.

[𝙷𝙴𝙿𝚏𝚒𝚝 website](https://hepfit.roma1.infn.it/)  
[𝙷𝙴𝙿𝚏𝚒𝚝 repository](https://github.com/silvest/HEPfit)  
[𝙷𝙴𝙿𝚏𝚒𝚝 paper/manual](https://arxiv.org/abs/1910.14012)  
[𝙷𝙴𝙿𝚏𝚒𝚝 documentation](https://hepfit.roma1.infn.it/doc/latest-release/index.html)  

### GAMBIT

Welcome to the GAMBIT homepage. GAMBIT is a global fitting code for generic Beyond the Standard Model theories, designed to allow fast and easy definition of new models, observables, likelihoods, scanners and backend physics codes.


https://gambit.hepforge.org/

### HiggsBounds

HiggsBounds takes a selection of Higgs sector predictions for any particular model as input and then uses the experimental topological cross section limits from Higgs searches at LEP, the Tevatron and the LHC to determine if this parameter point has been excluded at 95% C.L..

[HiggsBounds Git Repo](https://gitlab.com/higgsbounds/higgsbounds)

https://higgsbounds.hepforge.org/

### HiggsSignals

HiggsSignals performs a statistical test of the Higgs sector predictions of arbitrary models (using the HiggsBounds input routines) with the measurements of Higgs boson signal rates and masses from the Tevatron and the LHC.

[HiggsSignals Git Repo](https://gitlab.com/higgsbounds/higgssignals)

https://higgsbounds.hepforge.org/

### GFitter

The software package consists of abstract object-oriented code in C++ using ROOT functionality. Tools for the handling of the data, the fitting, and statistical analyses such as toy Monte Carlo sampling are provided by a core package, where theoretical errors, correlations, and inter-parameter dependencies are consistently dealt with. Theoretical models are inserted as plugin packages, which may be hierarchically organised. The use of dynamic parameter caching avoids the recalculation of unchanged results between fit steps, and thus significantly reduces the amount of computing time required for a fit.

http://gfitter.desy.de/Standard_Model/

### Contur

Exploring the sensitivity of unfolded collider measurements to BSM models. The manual for Contur 2.0, the first general user release, is available here: [arxiv link](https://arxiv.org/abs/2102.04377).

https://hepcedar.gitlab.io/contur-webpage/

### xFitter 

(former HERAFitter)
Proton parton distribution functions (PDFs) are essential for precision physics at the LHC and other hadron colliders. The determination of the PDFs is a complex endeavor involving several physics process. The main process is the lepton proton deep-inelastic scattering (DIS), with data collected by the HERA ep collider covering a large kinematic phase space needed to extract PDFs. Further processes (fixed target DIS, ppbar collisions etc.) provide additional constraining powers for flavour separation. In particular, the precise measurements obtained or to come from LHC will continue to improve the knowledge of the PDF.
The xFitter project is an open source QCD fit framework ready to extract PDFs and assess the impact of new data. The framework includes modules allowing for a various theoretical and methodological options, capable to fit a large number of relevant data sets from HERA, Tevatron and LHC. This framework is already used in many analyses at the LHC.

https://www.xfitter.org/xFitter/

### NNPDF fitting code

The NNPDF collaboration utilizes machine learning techniques to perform fits of Parton Distribution Functions.
All the code used to produce the latest family of PDFs from NNPDF, [NNPDF4.0](https://inspirehep.net/literature/1918284) has been made available open source as an [open-source package](https://github.com/NNPDF/nnpdf) together with the user-friendly examples and an extensive [documentation](https://docs.nnpdf.science/).

The code, based on Python and TensorFlow, can be used to produce the ingredients needed for PDF fits, to run the fits themselves, and to analyse the results. This framework, used to produce global PDF fits, allows for detailed external validation and reproducibility of the NNPDF4.0 analysis. Moreover, the code enables the user to explore a number of phenomenological applications, such as the assessment of the impact of new experimental data on PDFs, the effect of changes in theory settings on the resulting PDFs and a fast quantitative comparison between theoretical predictions and experimental data over a broad range of observables.

https://docs.nnpdf.science/

### LHAPDF

LHAPDF is a general purpose C++ interpolator, used for evaluating PDFs from discretised data files. Previous versions of LHAPDF were written in Fortran 77/90 and are documented at http://lhapdf.hepforge.org/lhapdf5/. LHAPDF6 vastly reduces the memory overhead of the Fortran LHAPDF (from gigabytes to megabytes!), entirely removes restrictions on numbers of concurrent PDFs, allows access to single PDF members without needing to load whole sets, and separates a new standardised PDF data format from the code library so that new PDF sets may be created and released easier and faster. The C++ LHAPDF6 also permits arbitrary parton contents via the standard PDG ID code scheme, is computationally more efficient (particularly if only one or two flavours are required at each phase space point, as in PDF reweighting), and uses a flexible metadata system which fixes many fundamental metadata and concurrency bugs in LHAPDF5. Compatibility routines are provided as standard for existing C++ and Fortran codes using the LHAPDF5 and PDFLIB legacy interfaces, so you can keep using your existing codes. But the new interface is much more powerful and pleasant to work with, so we think you'll want to switch once you've used it! LHAPDF6 is documented in more detail in http://arxiv.org/abs/1412.7420

https://lhapdf.hepforge.org

### PDFflow

PDFflow is a parton distribution function interpolation library written in Python and based on the TensorFlow framework.
It implements the same interpolation algorithms as LHAPDF with a focus on speed and parallelization.
The usage of the TensorFlow library allows for the usage of hardware accelerators (like GPUs) while writing hardware-agnostic code.

https://pdfflow.readthedocs.io/

### flavio

flavio is a Python package to compute observables in flavour physics, electroweak precision tests, Higgs physics, and other precision tests of the Standard Model, both in the Standard Model and in the presence of new physics encoded in Wilson coefficients of dimension-6 operators.

[flavio website](https://flav-io.github.io/)

[flavio repository](https://github.com/flav-io/flavio)

[flavio paper](https://arxiv.org/abs/1810.08132)

***

## Spectrum Generators

### Spheno

SPheno stands for S(upersymmetric) Pheno(menology). The code calculates the SUSY spectrum using low energy data and a user supplied high scale model as input. The spectrum is used to calculate two- and three body decay modes of supersymmetric particle as well as of Higgs bosons. In addition the production cross sections for supersymmetric particle and Higgs bosons in e^+ e^- annihilation is calculated. Moreover, the branching of the decay $b \to s \gamma$, the SUSY contribution to anomalous magnetic moment of the muon as well as the SUSY contributions to the rho parameter due to sfermions are calculated. The code is written in F90 with an emphasis on easy generalisability. The structure is set such that complex phases as well as the extension to include the flavour structure can be done in a straight forward way. The 2-loop renormalization group equations as well as the one-loop finite corrections a la Bagger, Matchev, Pierce and Zhang are included. In addition the two-loop corrections to the neutral Higgs boson masses (a la Brignole, Degrassi, Slavich and Zwirner) and to the mu-parameter (a la Dedes and Slavich) are included. Starting with version 2.2.2 the SUSY Les Houches Accord is supported as well as the SPA conventions (for details see hep-ph/0511344).

https://spheno.hepforge.org/

### SoftSUSY

This program provides a SUSY spectrum in the NMSSM, or the MSSM including flavour violation and with or without R-parity consistent with input Standard Model fermion mass/mixings and electroweak/strong coupling data. The R-parity violating mode can calculate neutrino masses and mixings to 1 loop. SOFTSUSY can be used in conjunction with other programs for many different particle physics calculations: see a SUSY tools review.

https://softsusy.hepforge.org/

### SuSpect

The public Fortran code SuSpect [79] calculates the supersymmetric and Higgs particle spec- trum in the Minimal Supersymmetric Standard Model (MSSM). In its present version (latest 2.41), it can deal with specific supersymmetry-breaking models with universal boundary con- ditions at high scales, such as the gravity (mSUGRA), anomaly (AMSB) or gauge (GMSB) mediated supersymmetry breaking models, as well as non-universal MSSM (restricted however to R–parity and CP conservation). Input and Output can be driven from the standard SLHA format files [15]. The algorithm includes the main mandatory ingredients such as the renormal- ization group evolution (RGE) of parameters between low and high energy scales, the consistent implementation of radiative electroweak symmetry breaking, and the calculation of the physi- cal masses of the Higgs bosons and supersymmetric particles including the full one-loop and dominant two-loop radiative corrections. In addition a control of important theoretical and ex- perimental features is available, such as the absence of non-physical minima, the amount of fine-tuning in the electroweak symmetry breaking condition, or the agreement with some preci- sion observables. Although SuSpect2 is still considered essentially up-to-date and will continue to be maintained in the future, a major upgrade is timely for several reasons. (copied from [this page](http://suspect.in2p3.fr/suspect3-LH.pdf))

http://suspect.in2p3.fr/updates.html

### ISAJET

ISAJET is a Monte Carlo program which simulates p p, pbar p, and e+ e- interactions at high energies. It is based on perturbative QCD plus phenomenological models for parton and beam jet fragmentation. Link to documentation can be found [here](http://www.nhn.ou.edu/~isajet/isajet788.pdf).

http://www.nhn.ou.edu/~isajet/


### TOP++

Purpose: The TOP++ program calculates the total inclusive cross-section for top-pair production at hadron colliders like the Tevatron and LHC. The program is capable of calculating the cross-section in fixed order QCD with exact NNLO. The program can also perform full NNLL soft gluon resummation. The resummation is done in Mellin space and then inverted numerically to x-space via the so-called Minimal Prescription.

http://www.precision.hep.phy.cam.ac.uk/top-plus-plus/

***

## Direct Detection

### DDCalc

Dark matter direct detection phenomenology package (DDCalc) is a software package for performing various dark matter direct detection calculations, including signal rate predictions and likelihoods for several experiments.

[DarkBit: A GAMBIT module for computing dark matter observables and likelihoods](https://arxiv.org/abs/1705.07920)

[Global analyses of Higgs portal singlet dark matter models using GAMBIT](https://arxiv.org/abs/1808.10465)

https://ddcalc.hepforge.org/

### WIMP Rates

Differential rates of WIMP-nucleus scattering in the standard halo model, for liquid xenon detectors.

https://github.com/JelleAalbers/wimprates

### NEST

NEST (Noble Element Simulation Technique) provides a simulation of the energy deposition-to-detector variable microphysics for liquid noble gas targets.

https://github.com/NESTCollaboration/nest
python bindings: 
https://github.com/NESTCollaboration/nestpy

### Flamedisx

Flamedisx aims to increase the practical number of dimensions and parameters in likelihoods for liquid-xenon (LXe) detectors. 
The LXe model is be computed with a series of (large) matrix multiplications, equivalent to the integral approximated by the MC simulation. Using TensorFlow makes this computation differentiable and GPU-scalable, so it can be used practically for fitting and statistical inference.

https://github.com/FlamTeam/flamedisx

### xenon1t_approximate_nuclear_recoil_likelihood

Binwise approximations of the XENON1T likelihood and XENONnT projections for fast inference on arbitrary models.

https://github.com/XENON1T/xenon1t_approximate_nuclear_recoil_likelihood

***

## Model Building
### SARAH

SARAH is a Mathematica package for building and analyzing SUSY and non-SUSY models. It calculates all vertices, mass matrices, tadpoles equations, one-loop corrections for tadpoles and self-energies, and two-loop RGEs for a given model. SARAH writes model files for FeynArts, CalcHep/CompHep, which can also be used for dark matter studies using MicrOmegas, the UFO format which is supported by MadGraph 5 and for WHIZARD and OMEGA. 
SARAH was also the first available spectrum-generator-generator: based on derived analytical expressions it creates source code for SPheno. It is therefore possible to implement new models in SPheno without the need to write any Fortran code by hand. The output for Vevacious can be used to check for the global minimum for a given model and parameter point.
Running SARAH is fast, it already includes a long list of SUSY and non-SUSY models, and the implementation of new models is efficient and straightforward.

[SARAH website](https://sarah.hepforge.org/)

### Sym2int

This page describes the Mathematica code Sym2Int (Symmetries to Interactions) which lists all valid interactions given the model's gauge group and fields (specified by their gauge and Lorentz representations). The program is valid for renormalizable interactions (mass dimension ≤4) as well as the ones which are not renormalizable (mass dimension >4). Since version 2, terms with derivatives and gauge bosons are also accounted for. More details can be found below. 

[Sym2int website](https://renatofonseca.net/sym2int)
[Sym2int paper 1](https://arxiv.org/abs/1703.05221)
[Sym2int paper 1](https://arxiv.org/abs/1907.12584)

## Susyno
Susyno is a Mathematica package which calculates the 2-loop renormalisation group equations of generic supersymmetric models, based on any gauge group (the only exception being multiple U(1) groups) and with any field content. 

[Susyno website](https://renatofonseca.net/susyno)

## Effective Field Theories

### DsixTools
DsixTools is a Mathematica package for the handling of the Standard Model Effective Field Theory (SMEFT) and the Low-energy Effective Field Theory (LEFT) with operators up to dimension six, both at the algebraic and numerical level.

[DsixTools website](https://dsixtools.github.io/)
[DsixTools Latest Paper](https://arxiv.org/pdf/2010.16341.pdf)

### smelli

smelli is a Python package providing a global likelihood function in the space of dimension-six Wilson coefficients in the Standard Model Effective Field Theory (SMEFT). The likelihood includes contributions from quark and lepton flavour physics, electroweak precision tests, and other precision observables.

The package is based on flavio for the calculation of observables and statistical treatment and wilson for the running, translation, and matching of Wilson coefficients.

[smelli website](https://smelli.github.io/)

[smelli repository](https://github.com/smelli/smelli)

[smelli paper](https://arxiv.org/abs/1810.07698)


## Feyn Family
(better name welcome)

### FeynRules

FeynRules is a Mathematica® package that allows the calculation of Feynman rules in momentum space for any QFT physics model. The user needs to provide FeynRules with the minimal information required to describe the new model, contained in the so-called model-file. This information is then used to calculate the set of Feynman rules associated with the Lagrangian. The Feynman rules calculated by the code can then be used to implement the new physics model into other existing tools, such as MC generators. This is done via a set of interfaces which are developed together and maintained by the corresponding MC authors.

https://feynrules.irmp.ucl.ac.be/

### FeynArts

FeynArts is a Mathematica package for the generation and visualization of Feynman diagrams and amplitudes.

http://www.feynarts.de/

### LoopTools

LoopTools is a package for evaluation of scalar and tensor one-loop integrals based on the FF package by G.J. van Oldenborgh. It features an easy Fortran, C++, and Mathematica interface to the scalar one-loop functions of FF and in addition provides the 2-, 3-, and 4-point tensor coefficient functions.

http://www.feynarts.de/looptools/

***

## Statistics

### pyhf

The HistFactory p.d.f. template [CERN-OPEN-2012-016] is per-se independent of its implementation in ROOT and sometimes, it’s useful to be able to run statistical analysis outside of ROOT, RooFit, RooStats framework.

https://github.com/scikit-hep/pyhf

### blueice

This package allows you to do parametric inference using likelihood functions, in particular likelihoods derived from Monte-Carlo or calibration sources. Especially when connected to a Monte Carlo, blueice lets you make likelihood functions which measure agreement between data and theory with flexibility: you choose which settings to vary (which parameters the likelihood functions has) and in which space the agreement is measured.

This package contains only generic code: you'll need a few things to make it useful for a particular experiment. Originally this code was developed for XENON1T only; the XENON1T models have since been split off to the laidbax repository.


https://github.com/JelleAalbers/blueice

python-based likelihood (in particular unbinned) construction/fitting framework used in XENON1T analyses. Linear template morphing, cached PDF generation

## Neutrino Code

Peter B. Denton's list for Neutrino code

https://peterdenton.github.io/Code/index.html


***

# COSMO

Note: Can be sorted better as the list grows, especially the general resources section

## General Cosmology Resources

### CosmoloPy

CosmoloPy is a package of cosmology routines built on NumPy/SciPy. Capabilities include: various cosmological densities, cosmological distance measures, galaxy luminosity functions (Schecter functions), conversion in and out of the AB magnitude system, pre-defined sets of cosmological parameters (e.g. from WMAP), perturbation theory and the power spectrum and reionization of the IGM.

http://roban.github.io/CosmoloPy/

### cosmopit

[Cosmopit](https://github.com/lontelis/cosmopit) (Cosmological python initial toolkit) is a package of cosmology built on basic numerical libraries of python such as Numpy/SciPy/Matplolib. Capabilities include common and not-so-common capabilities. 

Common capabilities include: various cosmological distances, and volumes calculations, transformation of cosmological coordinate systems, simple statistic quantities and also calculations of covariance matrices, fitting routines and bayesian inferences using corner plots, routines to read and write date in .fits and .txt files, example of calculating Monte Carlo Markov Chains for fittings and corner plot presentation tutorial, and a tutorial dedicated for transformation of cosmological coordinates.

Not-so common capabilities include: calculation of power spectra and correlation functions via the CLASS software, calculation routines of number counts, fractal dimension quantities, homogeneity scales, routines calculating statistics of simualated data, reading and writing routines for cute, astropy, pymangle, 
example of Baryon Acoustic Oscilation (BAO) detection and extraction feature from data, example of fitting the SuperNovae (SN) curve to extra simple magnitude parameters.

https://github.com/lontelis/cosmopit

### xAct and xPand

xAct is a suite of free packages for tensor computer algebra for Wolfram Mathematica. xAct implements state-of-the-art algorithms for fast manipulations of indices and has been modelled on the current geometric approach to General Relativity. It is highly programmable and configurable. Since its first public release in March 2004, xAct has been intensively tested and has solved a number of hard problems in GR.

xPand is a subpackage of the xAct distribution for efficient tensor manipulations a package for Mathematica. xPand provides tools to compute formally the cosmological perturbations for any tensor in any order around a homogeneous spacetime. xPand supports the most used gauges such as newtonian, comoving, synchronous, etc.

http://www.xact.es/
http://www2.iap.fr/users/pitrou/xpand.htm

### SageManifolds

The SageManifolds project aims at extending the modern Python-based computer algebra system SageMath towards differential geometry and tensor calculus.

SageManifolds deals with differentiable manifolds of arbitrary dimension. Various coordinate charts and vector frames can be introduced on the manifold, which does not need to be parallelizable. A given tensor field is then described by its sets of components in each vector frame, with automatic change-of-frame transformations for overlapping vector frames.

Generic pseudo-Riemannian manifolds can be considered, among which Riemannian manifolds and Lorentzian manifolds, with applications to General Relativity. In particular, the computation of the Riemann curvature tensor and associated tensors (Ricci, Weyl, Schouten and Cotton tensors) is implemented. SageManifolds can also deal with generic affine connections, not necessarily Levi-Civita ones. 

https://sagemanifolds.obspm.fr/

### HEALPix

HEALPix is a Hierarchical, Equal Area, and iso-Latitude Pixelation of the sphere designed to support efficiently (1) local operations on the pixel set, (2) a hierarchical tree structure for multi-resolution applications, and (3) the global Fast Spherical Harmonic transform. HEALPix based mathematical software meets the challenges presented by high resolution and large volume data sets, such as the WMAP and Planck CMB mission products. It's a must have for those who work with CMB data directly, allowing you to treat raw temperature data and easily extracting the Cls. It's implemented in C, C++, Fortran, Python, Java and IDL.

https://healpix.sourceforge.io/

Its Python implementation, Healpy, can be installed separately. It is very well documented and has tutorials to get started, which is great for those who want to start working with CMB data.

https://healpy.readthedocs.io/en/latest/index.html#

### EinsteinPy

EinsteinPy is an open source pure Python package dedicated to problems arising in General Relativity and gravitational physics, such as geodesics plotting for Schwarzschild, Kerr and Kerr Newman space-time model, calculation of Schwarzschild radius, calculation of Event Horizon and Ergosphere for Kerr space-time. Symbolic Manipulations of various tensors like Metric, Riemann, Ricci and Christoffel Symbols is also possible using the library. EinsteinPy also features Hypersurface Embedding of Schwarzschild space-time, which will soon lead to modelling of Gravitational Lensing.

https://einsteinpy.org/

### CosmoTransitions

The CosmoTransitions package is a set of python modules for calculating properties of effective potentials with one or more scalar fields. Most importantly, it can be used to find the instanton solutions which interpolate between different vacua in a given theory, allowing one to determine the probability for a vacuum transition.

https://github.com/clwainwright/CosmoTransitions

### CCL

The Core Cosmology Library (CCL) is a standardized library of routines to calculate basic observables used in cosmology. It will be the standard analysis package used by the LSST Dark Energy Science Collaboration (DESC). Please check the user policy!

https://github.com/LSSTDESC/CCL

### xmds

**Unsure if this belongs in cosmo really, will post an issue about it. I know some cosmo people use it for BEC etc**

This website provides the documentation for XMDS2 (an all-new version of XMDS), a software package that allows the fast and easy solution of sets of ordinary, partial and stochastic differential equations, using a variety of efficient numerical algorithms. If you publish work that has involved XMDS2, please cite it as Comput. Phys. Commun. 184, 201-208 (2013).

http://www.xmds.org

### Dark Emulator

A repository for a cosmology tool dark_emulator to emulate halo clustering statistics. The code is developed based on Dark Quest simulation suite (https://darkquestcosmology.github.io/). The current version supports the halo mass function and two point correlation function (both halo-halo and halo-matter cross). Dark Quest is a cosmological structure formation simulation campaign by Japanese cosmologists initiated in 2015. The primary goal of the project is to understand the complex parameter dependence of various large-scale structure probes, and provide a versatile tool to make predictions for parameter inference problems with observational datasets. The first series of simulations, Dark Quest. I. (DQ1), was completed in 2018 and we are now in the second phase (DQ2). A Gaussian-Process based emulation tool, Dark Emulator, was developed with the DQ1 database.

https://github.com/DarkQuestCosmology/dark_emulator_public

### Colossus

Colossus is a python toolkit for cosmology, large-scale structure, and dark matter halos. The main design goals are intuitive use and performance; the code is extensively documented. Colossus consists of three top-level modules. The cosmology module handles LCDM cosmologies with curvature, relativistic species, different dark energy equations of state, and so on. It includes densities, times, power spectra, variance, and correlation functions, among others. The large-scale structure module deals with peaks in Gaussian random fields and the statistical properties of halos such as peak height, halo bias, and the mass function. The halo module deals with masses, density profiles, concentration, and other halo properties. Colossus contains numerous fitting functions from the literature for convenience. 

https://bdiemer.bitbucket.io/colossus/

### CalPriorSNIa

[CalPriorSNIa](https://github.com/valerio-marra/CalPriorSNIa) computes the effective calibration prior on the absolute magnitude M_B of Type Ia supernovae that corresponds to a given determination of H_0. See Camarena & Marra [arXiv:1906.11814](https://arxiv.org/abs/1906.11814) and [arXiv:2101.08641](https://arxiv.org/abs/2101.08641) for more details. Example use case comments in two articles dated [19-Mar-2021](https://www.sunnyvagnozzi.com/blog/what-is-the-hubble-tension-really) and [22-Jan-2021](https://www.sunnyvagnozzi.com/blog/top-arxiv-week-3-2021) by Sunny Vagnozzi. 

[https://github.com/valerio-marra/CalPriorSNIa](https://github.com/valerio-marra/CalPriorSNIa)

***

## Inflationary Cosmology

### ASPIC

Aspic is a collection of fast modern fortran routines for computing various observable quantities used in Cosmology from definite single field inflationary models. It is distributed as a scientific library and aims at providing an efficient, extendable and accurate way of comparing theoretical inflationary predictions with cosmological data. Aspic currently supports 70 models of inflation, and more to come!
By observable quantities, we currently refer to as the Hubble flow functions, up to second order in the slow-roll approximation, which are in direct correspondence with the spectral index, the tensor-to-scalar ratio and the running of the primordial power spectrum. The aspic library also provides the field potential, its first and second derivatives, the energy density at the end of inflation, the energy density at the end of reheating, and the field value (or e-fold value) at which the pivot scale crossed the Hubble radius during inflation. All these quantities are computed in a way which is consistent with the existence of a reheating phase.

[ASPIC Library](http://cp3.irmp.ucl.ac.be/~ringeval/aspic.html)

[ASPIC Paper](https://arxiv.org/abs/1303.3787)

### NumDynInflation

This is a collection of codes for numerically simulating background dynamics as well as first-order perturbations for single field inflationary models. Efficient codes along with a little bit of manual work make for a very fast and accurate algorithm to compute various quantities of interest in single field inflation.

[Repository](https://github.com/bhattsiddharth/NumDynInflation)

[Paper](https://arxiv.org/abs/2212.00529)

***

## Big Bang Nucleosynthesis

### ACROPOLIS (A generiC fRamework fOr Photodisintegration Of LIght elementS)

ACROPOLIS is a generic framework to calculate the evolution of the light-element abundances due to photodisintegration reactions induced by different BSM particles. With ACROPOLIS, the widely discussed cases of decays as well as annihilations can be run without prior coding knowledge within example programs. However, its modular structure also makes it possible to easily implement other BSM physics scenarios.

https://acropolis.hepforge.org

### AlterBBN

AlterBBN is a C program which computes the abundances of the elements predicted by Big-Bang nucleosynthesis (BBN). Different cosmological scenarios are implemented in AlterBBN, which can alter the BBN predictions. Also, AlterBBN is included in the SuperIso Relic package so that the alternative models can be tested using BBN constraints.

https://alterbbn.hepforge.org

### PArthENoPE

PArthENoPE is a FORTRAN77 numerical code which computes the abundances of light nuclides produced during Big Bang Nucleosynthesis. Starting from nuclear statistical equilibrium conditions the program solves the set of corresponding coupled ordinary differential equations, follows the departure from chemical equilibrium of nuclear species, and determines their asymptotic abundances as function of several input cosmological parameters as the baryon density, the number of effective neutrino species, the value of cosmological constant and the neutrino chemical potential.

http://parthenope.na.infn.it

### PRIMAT (PRImordial MATter)

PRIMAT is a Mathematica code dedicated to the computation of light elements abundances at the end of the Big-Bang Nucleosynthesis (BBN). It computes the evolution of these abundances in the first minutes after the Big-Bang so as to obtain the frozen values when the temperature dropped below 108 degrees. It allows to explore the dependence of light elements abundances on cosmological parameters, such as baryon density or number of neutrino species, but also to estimate from a Monte-Carlo method the uncertainty in these predictions due to uncertainties in nuclear reaction rates.

http://www2.iap.fr/users/pitrou/primat.htm

## Einstein-Boltzmann Numerical Solvers

### CAMB

CAMB is a Python and Fortran code for computing CMB, CMB lensing, lensing, galaxy count and dark-age 21cm power spectra, transfer functions and matter power spectra, and background cosmological functions. Its object-oriented structure makes it easy to modify the code and add models such as exotic dark matter/dark energy and modified gravity to the analysis. 

https://camb.info (preferred citation method here: https://cosmologist.info/cosmomc/cosmomc.bib)

### CLASS

The purpose of CLASS is to simulate the evolution of linear perturbations in the universe and to compute CMB and large scale structure observables. Its name also comes from the fact that it is written in object-oriented style mimicking the notion of class. Classes are a wonderful programming feature available e.g. in C++ and Python, but these languages are known to be less vectorizable/parallelizable than plain C (or Fortran), and hence potentially slower. For CLASS we choose to use plain C for high performances, while organizing the code in a few modules that reproduce the architecture and philosophy of C++ classes, for optimal readability and modularity.

http://class-code.net

### Bolt

⚡ Bolt is the first differentiable cosmological Boltzmann-Solver. Bolt.jl is a pure-Julia integrator for the linearized Einstein-Boltzmann system. Bolt.jl predicts the CMB TT, TE, EE power spectra and accurately computes the linear matter power spectrum, which agrees to less than percent accuracy with CAMB/CLASS. Bolt computes the derivatives of cosmological perturbations and associated observables (currently) using forward-mode automatic differentiation. As Bolt
is written in Julia, the code makes use of unicode characters to minimize the code-equation translation barrier for those who want to build new cosmological physics into Bolt.

https://github.com/xzackli/Bolt.jl


### MGCAMB

Modified Growth with CAMB (MGCAMB) is a patch for the Einstein Boltzmann solver CAMB that intrdouces phenomenological Modifications of Growth (MG) along with dynamical Dark Energy (DE). It includes several phenomenological parametrizations.

https://github.com/sfu-cosmo/MGCAMB

### EFTCAMB

EFTCAMB is a patch of the public Einstein-Boltzmann solver CAMB, which implements the Effective Field Theory approach to cosmic acceleration. The code can be used to investigate the effect of different EFT operators on linear perturbations as well as to study perturbations in any specific DE/MG model that can be cast into EFT framework. To interface EFTCAMB with cosmological data sets, we equipped it with a modified version of CosmoMC, namely EFTCosmoMC, creating a bridge between the EFT parametrization of the dynamics of perturbations and observations.

eftcamb.org

### hi_class

The hi_class code is an extension of the Einstein-Boltzmann solver CLASS including gravitational theories in the Horndeski scalar-tensor category. hi_class can work with effective descriptions (expansion history + alpha functions) or covariant Lagrangians (+ field initial conditions) and includes consistent initial conditions and a flexible quasi-static approximation scheme. It can compute background quantities and linear observables (transfer functions, 3D matter power spectra, angular spectra for CMB, lensing potentials or galaxy number counts including relativistic effects) and be readily interfaced with samplers (Montepython, Cobaya, CosmoSIS).

www.hiclass-code.net

### CLASS-PT

This is a modification of the CLASS code that computes the non-linear power spectra of dark matter and biased tracers in one-loop cosmological perturbation theory.

https://github.com/Michalychforever/CLASS-PT

### CLASS delens

Implementation of delensing procedure discussed in https://arxiv.org/abs/2111.15036: by authors Selim C. Hotinli, Joel Meyers, Cynthia Trendafilova, Daniel Green, Alex van Engelen

Code produces delensed CMB spectra (TT, TE, EE and BB) and lensing-reconstruction noise for given CMB experiment specifications and cosmology.

Delensing reverses the effects of lensing on the observed CMB temperature and polarization maps. This provides various benefits. Delensed CMB spectra have sharper acoustic peaks and more prominent damping tails, allowing for improved inferences of cosmological parameters that impact those features. Delensing reduces B-mode power, aiding the search for primordial gravitational waves and allowing for lower variance reconstruction of lensing and other sources of secondary CMB anisotropies. Lensing-induced power spectrum covariances are reduced by delensing, simplifying analyses and improving constraints on primordial non-Gaussianities. Please refer to https://arxiv.org/abs/2111.15036 for a detailed demonstration of the benefits of CMB delensing.

This code can be used as a submodule for the Fisher forecasting tool software https://github.com/ctrendafilova/FisherLens.

### FisherLens

This code provides a wrapper for the CLASS_delens code to facilitate Fisher forecasting of cosmological parameter constraints from CMB spectra.
https://github.com/ctrendafilova/FisherLens

CLASS_delens https://github.com/selimhotinli/class_delens

Authors: Selim C. Hotinli, Joel Meyers, Cynthia Trendafilova, Daniel Green, Alexander van Engelen

### GDM_CLASS

GDM_CLASS is a modification of the CLASS code which adds the phenomenologically motivated Generalized Dark Matter (GDM) model from [Hu (1998)](https://arxiv.org/abs/astro-ph/9801234). A detailed investigation of the theoretical aspects related to this model can be found in the [following article](https://arxiv.org/abs/1605.00649).

https://github.com/s-ilic/gdm_class_public

Authors: Stéphane Ilic, Michael Kopp, Daniel B. Thomas, Constantinos Skordis

***

## Cosmological Parameter Estimation and Statistical Analysis

### CosmoMC

CosmoMC is a Fortran 2008 Markov-Chain Monte-Carlo (MCMC) engine for exploring cosmological parameter space, together with Fortran and python code for analysing Monte-Carlo samples and importance sampling (plus a suite of scripts for building grids of runs, plotting and presenting results). The code does brute force (but accurate) theoretical matter power spectrum and Cl calculations with CAMB. See the [original paper](https://arxiv.org/abs/astro-ph/0205436) for an introduction and descriptions, and [up-to-date sampling algorithm](https://arxiv.org/abs/1304.4473) details. It can also be compiled as a generic sampler without using any cosmology codes.

https://cosmologist.info/cosmomc

### Monte Python

Monte Python is a Monte Carlo code for Cosmological Parameter extraction. It contains likelihood codes of most recent experiments, and interfaces with the Boltzmann code class for computing the cosmological observables.

https://github.com/brinckmann/montepython_public

### Cobaya

Cobaya (code for bayesian analysis, and Spanish for Guinea Pig) is a framework for sampling and statistical modelling: it allows you to explore an arbitrary prior or posterior using a range of Monte Carlo samplers (including the advanced MCMC sampler from CosmoMC, and the advanced nested sampler PolyChord). The results of the sampling can be analysed with GetDist. It supports MPI parallelization (and very soon HPC containerization with Docker/Shifter and Singularity).

https://cobaya.readthedocs.io/en/latest/

### CosmoSIS

CosmoSIS is a cosmological parameter estimation code. It is now at version 1.6. It is a framework for structuring cosmological parameter estimation in a way that eases re-usability, debugging, verifiability, and code sharing in the form of calculation modules. It consolidates and connects together existing code for predicting cosmic observables, and makes mapping out experimental likelihoods with a range of different techniques much more accessible. CosmoSIS is described in Zuntz et al.: http://arxiv.org/abs/1409.3409. If you make use of it in your research, please cite that paper and include the URL of this repository in your acknowledgments. Thanks!

https://bitbucket.org/joezuntz/cosmosis/wiki/Home

### CosmoLike

CosmoLike is a collaborative software development project to analyze cosmological data sets and to forecast future missions.

https://github.com/CosmoLike

### ECLAIR

The ECLAIR suite of codes is meant to be used as a general inference tool, allowing to sample via MCMC techniques the posterior distribution of a set of parameters corresponding to a particular physical model, under the constraint of a number of datasets/likelihoods. It also contains a robust maximizer aimed at finding the point in parameter space corresponding to the best likelihood of any considered model.The suite also include a plotting script allowing to conveniently diagnose and check the convergence of a chain, as well as produce summary statistics on the parameters of interest.

https://github.com/s-ilic/ECLAIR

Author: Stéphane Ilić

***

## Correlation function and higher order statistics

### Corrfunc

Corrfunc is a set of high-performance routines to measure two-point clustering statistics, in particular computing angular correlation functions, projected correlation functions and 3D correlation functions. This is written in C and Python and fully parallelized.

https://corrfunc.readthedocs.io/en/master/

### encore

encore is a C++ code for estimating isotropic 3-, 4-, 5- and 6-point correlation functions of galaxy surveys and simulations. This is uses spherical harmonic decompositions to give an algorithm with quadratic complexity, as in [Philcox et al. 2021](https://arxiv.org/abs/2105.08722), and includes corrections for the survey geometry. The code also features a GPU implementation via CUDA, and is actively being developed.

https://github.com/oliverphilcox/encore

### Pynkowski

[Pynkowski](https://github.com/javicarron/pynkowski) is a fully documented Python package to compute Minkowki Functional and other higher order statistics of input fields, as well as computing the expectation values for different kinds of fields. It can compute these statistics on different kinds of data, such as healpix maps (scalar and spin 2) and two- and three-dimensional numpy arrays. It can also compute the theoretical expectation for different fields, including Gaussian and $\chi^2$ isotropic fields. It includes Minkowski Functionals (also called Lipschitz-Killing Curvatures) and maxima/minima distributions. The code has been designed to easily support more data formats, theoretical fields, and statistics.

https://github.com/javicarron/pynkowski (Documentation [here](https://javicarron.github.io/pynkowski/pynkowski.html))

***



## Simulations

### Quijote Simulations

The Quijote simulations is a suite of 44,100 full N-body simulations designed to:

Quantify the information content on cosmological observables
Provide enough statistics to train machine learning algorithms

[Quijote Simulations ReadTheDocs](https://quijote-simulations.readthedocs.io/en/latest/)
[Quijote Simulations GH](https://github.com/franciscovillaescusa/Quijote-simulations)
[Citation](https://quijote-simulations.readthedocs.io/en/latest/citation.html)

### nbodykit

nbodykit is an open source project written in Python that provides a set of state-of-the-art, large-scale structure algorithms useful in the analysis of cosmological datasets from N-body simulations and observational surveys. All algorithms are massively parallel and run using the Message Passing Interface (MPI).

[nbodykit Website](https://nbodykit.readthedocs.io/en/latest/index.html)

### CAMELS

CAMELS stands for Cosmology and Astrophysics with MachinE Learning Simulations, and is a suite 4,233 cosmological simulations: 2,049 N-body and 2,184 state-of-the-art (magneto-)hydrodynamic.

The CAMEL simulations are described in detail in the CAMELS presentation paper, [https://arxiv.org/abs/2010.00619](2010.00619).

This repository contains the codes and scripts we wrote to carry out the analysis of the CAMELS presentation paper.


[CAMELS ReadTheDocs](https://camels.readthedocs.io/en/latest/)

[CAMELS GH Repository](https://github.com/franciscovillaescusa/CAMELS)

[CAMELS Website](https://www.camel-simulations.org)

[Citation](https://camels.readthedocs.io/en/latest/citation.html)

### CAMELS Multifield Dataset (CMD)

CMD is a publicly available collection of hundreds of thousands 2D maps and 3D grids containing different properties of the gas, dark matter, and stars from more than 2,000 different universes. The data has been generated from thousands of state-of-the-art (magneto-)hydrodynamic and gravity-only N-body simulations from the CAMELS project.

[CMD GH Repository](https://github.com/franciscovillaescusa/CMD)

[CMD ReadTheDocs](https://camels-multifield-dataset.readthedocs.io/en/latest/)

[CMD Paper](https://arxiv.org/abs/2109.10915)

[Citation](https://camels-multifield-dataset.readthedocs.io/en/latest/citation.html)



### N-GenIC

This is the N-GenIC code used to generate the N-body initial conditions.
The codes uses the Zel'dovich approximation to generate the IC.
The code can create the ICs for pure CDM simulations, simulations containing
CDM & Baryons, CDM & neutrinos and CDM Baryons & neutrinos.

[N-GenIC GH Repository](https://github.com/franciscovillaescusa/N-GenIC_growth)

### monofonIC

High order LPT/PPT cosmological IC generator for single resolution cosmological simulations.

BEWARE: This software package is currently in a beta testing phase! Proceed with caution before using it for large-scale simulations.

The project is currently still under active development, so check the documentation on which features are currently implemented and how new modules/plugins can be added.

[Repository](https://bitbucket.org/ohahn/monofonic/src/master/)

### CosMomentum

As of now, this code can compute
* the probability distribution function (PDF) of the matter density field,
* the cumulant generating function (CGF) of the matter density field,
* individual cumulants of the matter density field.
It does this for
* Gaussian and 3 types of non-Gaussian initial conditions.
* both line-of-sight projections of the matter density field as well as the 3D density field
* all of this for PDFs of lensing convergence
* statistics of biased & stochastic tracers (e.g. galaxies)

Cite as Friedrich et al. ([2019, arXiv:1912.06621](https://arxiv.org/abs/1912.06621)) if you use this code.
[Repository](https://github.com/OliverFHD/CosMomentum)



### Pylians

Pylians stands for Python libraries for the analysis of numerical simulations. They are a set of python libraries, written in python, cython and C, whose purposes is to facilitate the analysis of numerical simulations (both N-body and hydro). Among other things, they can be used to:

Compute density fields
Compute power spectra
Compute bispectra
Compute correlation functions
Identify voids
Populate halos with galaxies using an HOD
Apply HI+H2 corrections to the output of hydrodynamic simulations
Make 21cm maps
Compute DLAs column density distribution functions
Plot density fields and make movies
Pylians were the native or inhabitant of the Homeric town of Pylos.

Notice that Pylians only works with python2. For a python3 version see Pylians3.

From 2020 only Pylians3 will be supported/updated.

[Pylians GH Repository](https://github.com/franciscovillaescusa/Pylians)



### Pylians3

Pylians stands for Python libraries for the analysis of numerical simulations. They are a set of python libraries, written in python, cython and C, whose purposes is to facilitate the analysis of numerical simulations (both N-body and hydrodynamic). Pylians3 evolved from Pylians to support python3. Among other things, they can be used to:

Compute density fields
Compute power spectra
Compute bispectra
Compute correlation functions
Identify voids
Populate halos with galaxies using an HOD
Apply HI+H2 corrections to the output of hydrodynamic simulations
Make 21cm maps
Compute DLAs column density distribution functions
Plot density fields and make movies
Pylians were the native or inhabitant of the Homeric town of Pylos

[Pylians3 ReadTheDocs](https://pylians3.readthedocs.io/en/master/)

[PYlians3 GH Repository](https://github.com/franciscovillaescusa/Pylians3)


### Abacus Cosmos

N-body simulation data products from the 2017 Abacus Cosmos project, including halos catalogs, particle subsamples, power spectra, and initial conditions.

[Abacus Cosmos Website](https://lgarrison.github.io/AbacusCosmos/)

[Abacus Cosmos Paper](https://lgarrison.github.io/AbacusCosmos/papers/)


### Abacus Summit

AbacusSummit is a suite of large, high-accuracy cosmological N-body simulations. These simulations were designed to meet (and exceed!) the Cosmological Simulation Requirements of the Dark Energy Spectroscopic Instrument (DESI) survey. AbacusSummit was run on the Summit supercomputer at the Oak Ridge Leadership Computing Facility under a time allocation from the DOE’s ALCC program.

[Abacus Summit ReadtheDocs](https://abacussummit.readthedocs.io/en/latest/)

[Abacus Summit GH Repository](https://github.com/abacusorg/AbacusSummit)

[Papers and Citation](https://abacussummit.readthedocs.io/en/latest/citation.html)


### abacusutils

abacusutils is a package for reading and manipulating data products from the Abacus N-body project. In particular, these utilities are intended for use with the AbacusSummit suite of simulations. We provide multiple interfaces: primarily Python 3, but also C/C++ [coming soon!] and language-agnostic interfaces like Unix pipes.

[abacusutils GH repository](https://github.com/abacusorg/abacusutils)

[abacusutils ReadtheDocs](https://abacusutils.readthedocs.io/en/latest/)


### zeldovich-PLT

This code generates Zel'dovich approximation (ZA) initial conditions (i.e. first-order Lagrangian perturbation theory) for cosmological N-body simulations, optionally applying particle linear theory (PLT) corrections. This code does not provide second-order ICs (2LPT), but one can use these ICs with the config-space 2LPT detailed in Garrison et al. (2016). This is the primary IC generator used by the Abacus N-body code.

If you do not intend to use the config-space 2LPT, then it's better to use a Fourier-space 2LPT code (e.g. 2LPTic) than to rely on ZA, even with PLT corrections.

[zeldovich-PLT GH Repository](https://github.com/abacusorg/zeldovich-PLT)

[Citation](https://github.com/abacusorg/zeldovich-PLT#citation)


### 2LPTic 


2LPT initial conditions parallel and serial code. The MPI-parallel version of 2LPTic is based on the Zel'dovich n-genic code written by Volker Springel.

The serial code is written in Fortran77 and requires the installation of FFTW version 2.1.5. You can read Appendix D2 in this paper to understand how the code works. 

[2LPTic Website](https://cosmo.nyu.edu/roman/2LPT/)


### Aemulus Project

* All simulations: 1050 Mpc/h, 14003 particles, resolving halos down to ~6x1012 Msol/h
* 40 Training Simulations, all different cosmologies and independent phases
* 7x5=35 Test Sample Simulations: 7 cosmologies, with 5 realizations per cosmology
* Additional “convergence” simulations testing box size, particle resolution, etc.

[Aemulus Website](https://aemulusproject.github.io/index.html)

[Aemulus Papers](https://aemulusproject.github.io/papers.html)


### BAHAMAS Project

BAHAMAS is a first attempt to do large-scale structure (LSS) cosmology using self-consistent full cosmological hydrodynamical simulations. Such simulations simultaneously follow all of the important matter components and allow us to accurately address the main uncertainty in the theoretical modelling of LSS, which is astrophysical effects (such as feedback) associated with the baryonic component.

A lay person's description of the project can be found in an [article](https://www.ljmu.ac.uk/about-us/news/articles/2017/12/6/funding-to-uncover-dark-side-of-the-universe) put together by LJMU's media team.

The initial BAHAMAS paper, which describes in detail our feedback modelling strategy, is McCarthy et al. (2017).

With a viable model for astrophysical effects, we are now gearing up to explore variations in cosmology. Both in the context of the standard model (LCDM), but also plausible extensions of LCDM, including massive neutrino cosmologies, time-variable dark energy models, and modified gravity models.

Recently, we posted our initial results on constraints on the summed mass of neutrinos in McCarthy et al. (2018). A popular science article about our results can be found in The Conversation.

[BAHAMAS Project Webpage](https://www.astro.ljmu.ac.uk/~igm/BAHAMAS/)

[BAHAMAS Paper 1](https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.2936M/abstract)

[BAHAMAS Paper 2](https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.2999M/abstract)


### CosmoSim

The CosmoSim database provides results from cosmological simulations performed within different projects: MultiDark and Bolshoi, CLUES, and Galaxies.

[CosmoSim database](https://www.cosmosim.org)

[MultiDark and Bolshoi Project](https://www.cosmosim.org/cms/documentation/projects/multidark-bolshoi-project/)

[CLUES Project](https://www.cosmosim.org/cms/documentation/projects/clues-project/)

[Galaxies](https://www.cosmosim.org/cms/documentation/projects/galaxies/)


### Massive Galaxy Project

The Massive Galaxy Project is a multi-pronged effort to characterize and understand the physical processes that drive the evolution of the most massive galaxies in the universe.
 
Our goals are to answer questions such as, are massive galaxies fundamentally multi-component systems and when did their components assemble?  What is their growth rate over cosmic time?  How do they reflect and influence their environment?
 
To tackle this broad range of questions, we make use of a variety of data sets from COSMOS, BOSS, the Stripe 82 co-adds, and IFU observations from the MaNGA survey.
 
The Massive Galaxy Project includes the following results and data-products:
 
* The Stripe 82 Massive Galaxy Catalog (S82-MGC), including a mass limited sample of 41,770 galaxies with log Mstar/Msun>11.2 to z=0.7
* Stellar mass completeness estimates for the BOSS CMASS and LOWZ samples
* Mock catalogs for the BOSS survey which reproduce the two-point correlation function and which account for the stellar mass completeness of the BOSS survey.

[The Massive Galaxy Project Website](https://www.ucolick.org/~kbundy/massivegalaxies/index.html)

[Publications](https://www.ucolick.org/~kbundy/massivegalaxies/publications.html)


### CurieHZ 

Large Scale simulations of Ly-α  and Ly-break galaxies in the high-z universe: Probing the epoch of reionization.

The  main goal of this project was to  generate a mass limited complete volume sample of resimulated objects at high-z (from z≈ 3 to 6).  To this end, we have simulated a cubic box  200 h-1 Mpc with different resolutions and physics, ranging from 5123 to 20483 particles.   These simulations of the whole box have been used to identify the objects of interest  in order to generate the  initial conditions for the resimulations.  We made use of a new Initial Conditions Generation code called GINNUNGAGAP  written by Steffen Knollmann at UAM.  This code has some new features such as  variable domain decompositions (slab or pencil) for the Fast Fourier Transform in 3D,  fully parallel random number generator and  the realizations of  random gaussian fluctuations in real space, rather than in Fourier space.  Increasing levels of mass refinements  can easily be  generated in this way to almost unlimited number, provided that there are sufficiente computational MPI nodes.   Thanks to the pencil domain decomposition, there is no limitation in the mesh sizes used to compute the density fields and to FFT them.

[CurieHZ Project Website](https://curiehz.ft.uam.es)


### EAGLE 

EAGLE (Evolution and Assembly of GaLaxies and their Environments) is a simulation aimed at understanding how galaxies form and evolve.  This computer calculation models the formation of structures in a cosmological volume, 100 Megaparsecs on a side (over 300 million light-years). This is large enough to contain 10,000 galaxies of the size of the Milky Way or bigger, enabling a comparison with the whole zoo of galaxies visible in the Hubble Deep field for example.

[EAGLE Webpage](http://icc.dur.ac.uk/Eagle/index.php)

[Publications](http://icc.dur.ac.uk/Eagle/publications.php)

### Aquarius

The Aquarius Project is a large-scale collaborative programme of the Virgo Consortium, similar in scope and scale to the Millennium Simulation project. At present, the principal set of Aquarius simulations contains six examples of an isolated halo similar in mass to that of the Milky Way. These are simulated in their full cosmological context (assuming the concordance LCDM cosmology) and at various resolutions up to about 200 million particles (counted within the radius where the enclosed density is 200 times the cosmic mean). One halo is also simulated at even higher resolution, resulting in almost 1.5 billion particles within this radius. These simulations are being used to understand the fine-scale structure predicted around the Milky Way by the standard structure formation model, and as the basis for simulation by various techniques of the growth of the stellar components of our Galaxy.

[Aquarius Project Website](https://wwwmpa.mpa-garching.mpg.de/aquarius/)

### Auriga


The Auriga Project is a large suite of high-resolution magneto-hydrodynamical simulations of Milky Way-sized galaxies, simulated in a fully cosmological environment by means of the 'zoom-in' technique. It is simulated with the state-of-the-art hydrodynamic moving mesh code AREPO, and includes a comprehensive galaxy formation model based on the successful cosmological simulation Illustris.

[Auriga Website](https://wwwmpa.mpa-garching.mpg.de/auriga/)

[Auriga Papers](https://wwwmpa.mpa-garching.mpg.de/auriga/publications.html)


### Illustris

The Illustris project is a large cosmological simulation of galaxy formation, completed in late 2013, using a state of the art numerical code and a comprehensive physical model. Building on several years of effort by members of the collaboration, the Illustris simulation represents an unprecedented combination of high resolution, total volume, and physical fidelity. 

[Illustris Website](https://www.illustris-project.org)

[Illustris Papers](https://www.illustris-project.org/results/)

### IllustrisTNG

The IllustrisTNG project is an ongoing series of large, cosmological magnetohydrodynamical simulations of galaxy formation. TNG aims to illuminate the physical processes that drive galaxy formation: to understand when and how galaxies evolve into the structures that are observed in the night sky, and to make predictions for current and future observational programs. The simulations use a state of the art numerical code which includes a comprehensive physical model and runs on some of the largest supercomputers in the world. TNG is a successor to the original Illustris simulation and builds on several years of effort by many people.

[IllustrisTNG Website](https://www.tng-project.org)

[TNG Papers](https://www.tng-project.org/results/)

### CosmoDC2

CosmoDC2 is a large synthetic galaxy catalog designed to support precision dark energy science with LSST, covering 440 sq. deg. of sky area to a redshift of z = 3, with a magnitude depth of 28 in the r band. A wide range of galaxy properties are available in cosmoDC2. To learn more about cosmoDC2, please see Korytov et al. (LSST DESC), ApJS, 245, 26 (2019).

[CosmoDC2 Website](https://data.lsstdesc.org/doc/cosmodc2)

### DC2 Simulated Sky Survey

The DC2 Simulated Sky Survey is a 300-sq-deq simulated survey in six optical bands with observations following a reference LSST observing cadence.

[DC2 SSS Website](https://data.lsstdesc.org/doc/dc2_sim_sky_survey)

### Magneticum simulations

The Magneticum simulations are aiming to follow the formation of cosmological structures over a large range of cosmological scales by performing a set of hydrodynamical simulations of different cosmological volumes, each of them sampled with a very large number of particles providing excellent spacial resolution of the different simulations. We take many physical processes into account to allow detailed comparisons to a variety of multi-wavelength observational data. See Simulations for a detailed description for astrophysicists.

[Magenticum Website](http://www.magneticum.org/index.html)

[Magneticum Publications](http://www.magneticum.org/publications.html)

### MICE

It includes extensive galaxy and lensing properties for ~200 million galaxies over 5000 sq.deg and up to a redshift z=1.4
This galaxy mock is built from the parent MICE-Grand Challenge (MICE-GC) lightcone N-body simulation and associated halo catalog, 
that is throughly validated in a series of forthcoming papers: Fosalba et al. 2013a,b, Crocce et al. 2013 
Further details about the galaxy mock method implemented will be given in Castander et al. 2014, Carretero et al. 2014.


[MICE Website](http://maia.ice.cat/mice/)

### CosmicProfiles

This Cython library provides shape and density profile analysis tools for point clouds such as halos and galaxies in cosmological simulations and beyond. Its features include overall halo shape determination (i.e. major, intermediate, minor axis vectors and shape quantities such as intermediate-to-major axis ratio or sphericity), halo shape profile determination, halo density profile estimation using direct binning and kernel-based approaches, density profile fitting functionalities assuming either an NFW, Hernquist, Einasto or alpha-beta-gamma (generalized NFW) profile model, a mock halo generator etc.


[CosmicProfiles ReadTheDocs](https://cosmic-profiles.readthedocs.io/en/latest/)

[CosmicProfiles GH Repository](https://github.com/tibordome/cosmic_profiles)


***

## Extended Cosmologies

### ReACT

ReACT is an extension of the software package Copter (0905:0479) and MG-Copter (1606.02520) which allows for the calculation of redshift and real space large scale structure observables for a wide class of gravity and dark energy models.

[ReACT GH Repository](https://github.com/nebblu/ReACT)

### pyLDT-cosmo

A Python package to generate matter PDF predictions in Large Deviation Theory for ΛCDM and alternative cosmologies.
Currently available cosmological models include:
* ΛCDM
* Hu-Sawicki f(R) gravity with ΛCDM background (0705.1158)
* nDGP gravity with ΛCDM background (0910.0235)
* w0waCDM (0009008, 0208512, 0808.3125)

pyLDT-cosmo is partly based on the Julia programming language, which uses a Just-In-Time (JIT) compiler to improve runtime performance. However, this feature comes at the expense of slow library loading when compared to pure or pre-compiled Python packages. On a modern computer pyLDT-cosmo takes about 80 seconds to load. After that computation is very efficient, taking only ~3 seconds to derive the matter PDF's simultaneuosly for all models, redshifts and smoothing radii.

[GH Repository](https://github.com/mcataneo/pyLDT-cosmo)

[Jupyter Notebook](https://github.com/mcataneo/pyLDT-cosmo/tree/main)

[Citation](https://github.com/mcataneo/pyLDT-cosmo#citation)

***

## Perturbation Theory

### Copter

The links below provide a tar-ball and README file for Copter, a C++ class library written by Jordan Carlson to compute predictions for different cosmological perturbation theory schemes. The name Copter loosely stands for cosmological perturbation theory. Specifically, Copter includes code for computing statistical observables in the large-scale structure of matter using various forms of perturbation theory, including linear theory, standard perturbation theory, renormalized perturbation theory, and many others.

[Copter Website](http://mwhite.berkeley.edu/Copter/)

[Copter Paper](https://arxiv.org/abs/0905.0479)

[Citation](https://arxiv.org/abs/0905.0479)

### nuCopter

Matter Power Spectrum and RSD Predictions including Massive Neutrinos and Dynamical Dark Energy.
Extended version of Copter.

[nuCopter Page](https://www.hep.anl.gov/cosmology/Pert/README.html)

Citation: [1](https://arxiv.org/abs/astro-ph/9911177), [2](https://arxiv.org/abs/0905.0479), [3](https://arxiv.org/abs/1309.5872)

### CLPT and GSRSD

Convolution Lagrangian Perturbation Theory (CLPT) and Gaussian Streaming Redshift-Space Distortions (GSRSD)

CLPT is a C++ code written by Lile Wang to compute predictions for the halo and matter clustering using Convolution Lagrangian Perturbation Theory; described in Carlson, Reid & White [MNRAS 429(2013)1674, arxiv:1209:0780]. It also has extensions for the Gaussian streaming model, described in described in Reid & White [MNRAS 417(2013)1913, arxiv:1105:4165] and Wang, Reid & White [MNRAS 437(2014)588, arxiv:1306:1804]

[CLPT GH Repository](https://github.com/wll745881210/CLPT_GSRSD)

[CLPT Page](http://mwhite.berkeley.edu/CLPT/)

[Citation](https://arxiv.org/abs/1306.1804)


### CLEFT-GSM

This code implements the Gaussian Streaming Model using components from Convolution Lagrangian Effective Field Theory.

The code is written (mostly) in C++. It can be run from the command line, or called from Python (wrappers provided).

[GH Repository](https://github.com/martinjameswhite/CLEFT_GSM)

[CLEFT-GSM Paper](https://arxiv.org/abs/1609.02908)

### velocileptors

Velocity-based perturbation theory (both Lagrangian (LPT) and Eulerian (EPT) formulations) expansions of redshift-space distortions and velocity statistics.

This code computes the real- and redshift-space power spectra and correlation functions of biased tracers using 1-loop perturbation theory (with effective field theory counter terms and up to cubic biasing) as well as the real-space pairwise velocity moments.

[Velocileptors GH Repository](https://github.com/sfschen/velocileptors)

Citation: [1](https://arxiv.org/abs/2005.00523), [2](https://arxiv.org/abs/2012.04636)

***

## Baryo and Leptogenesis

### ULYSSES

ULYSSES (Universal LeptogeneSiS Equation Solver) is a python package that calculates the baryon asymmetry produced from leptogenesis in the context of a type-I seesaw mechanism. The code provides: a selection of predefined Boltzmann equations with one, two and three decaying right-handed neutrinos in the single, two and three flavour regime; density matrix equations to fully take flavour effects into account, as well as oscillations of right-handed neutrinos for the ARS leptogenesis scenario; complete Boltzmann equations for thermal leptogenesis applying proper quantum statistics without assuming kinetic equilibrium of the right-handed neutrinos; primordial black hole-induced leptogenesis. Furthermore, the ULYSSES code provides tools for multi-dimensional parameter space exploration. The emphasis of the code is on user flexibility and rapid evaluation. 

[GH Repository](https://github.com/earlyuniverse/ulysses)

Citation: [1](https://arxiv.org/abs/2007.09150), [2](https://arxiv.org/abs/2301.05722)

### AMIQS

Fast code written in C++ that calculates the baryon asymmetry generated via quantum oscillations (ARS mechanism) of right-handed neutrinos in the early Universe.

[GH Repository] https://github.com/stefanmarinus/amiqs

Citation: [1](https://arxiv.org/abs/2207.01651)

# ASTRO

## Black hole perturbation theory

### Black hole perturbation toolkit

The Black Hole Perturbation Toolkit brings together software and data relating to black hole perturbation theory. These can then be used to model gravitational radiation from small mass-ratio binaries as well as from the ringdown of black holes. The former are key sources for the future space-based gravitational wave detector, LISA.

Our overall goal is for less researcher time to be spent writing code and more time spent doing physics. Currently there exist multiple scattered black hole perturbation theory codes developed by a wide array of individuals or groups over a number of decades. This project aims to bring together some of the core elements of these codes into a Toolkit that can be used by all.

The BHPToolkit is made up of many different tools which can be
individually installed by users depending on what they are interested
in.  Currently around ~20 packages.

* Web: http://bhptoolkit.org/
* Repos: https://github.com/BlackHolePerturbationToolkit
* Preferred citation method: https://bhptoolkit.org/BHPToolkit.bib

***

## Cosmic Rays

### Dragon

DRAGON adopts a second-order Cranck-Nicholson scheme with Operator Splitting and time overrelaxation to solve the diffusion equation. This provides fast a solution that is enough accurate for the average user. Occasionally, users may want to have very accurate solutions to their problem. To enable this feature, users may get close to the accurate solution by using the fast method, and then switch to a more accurate solution scheme, featuring the Alternating-Direction-Implicit (ADI) Cranck-Nicholson scheme.

Some parts of DRAGON are built following GALPROP, v50p. The first reason is that it is a waste of time to reimplement standard parts, like energy losses, in which nothing new has to be found. The second reason is that it is essential to be able to compare our predictions with that of the Galprop code, and this can be done only by following the details of its implementation. Therefore, we kept in the code some features and models used in Galprop, like nuclear cross-sections, the gas distribution, the convergence technique. However, each of these models is accompanied by other models, which can be selected by setting the appropriate switch. This is done very easily using the well known C++ structure of abstract/derived classes. The code is then very flexible and easy to manage and to modify or update.

https://github.com/cosmicrays

### USINE

A library with several semi-analytical Galactic cosmic-ray (GCR) propagation models. [Link to documentation](https://dmaurin.gitlab.io/USINE/_downloads/05d002129a26a2b732bd2a7d44e08ed4/usine.pdf)

https://dmaurin.gitlab.io/USINE/

### GALPROP

GALPROP is a numerical code for calculating the propagation of relativistic charged particles and the diffuse emissions produced during their propagation. The GALPROP code incorporates as much realistic astrophysical input as possible together with latest theoretical developments. The code calculates the propagation of cosmic-ray nuclei, antiprotons, electrons and positrons, and computes diffuse γ-rays and synchrotron emission in the same framework. Each run of the code is governed by a configuration file allowing the user to specify and control many details of the calculation. Thus, each run of the code corresponds to a potentially different 'model'.

https://galprop.stanford.edu/

### PICARD

Picard is a Galactic cosmic ray propagation code developed at Innsbruck University. The purpose of the code is the numerical solution of the cosmic ray transport equations with a focus on the observed cosmic ray spectra at Earth and the gamma-ray emission resulting from the interaction of the Galactic cosmic rays with the interstellar medium.

https://astro-staff.uibk.ac.at/~kissmrbu/Picard.html

### CORSIKA

CORSIKA (COsmic Ray SImulations for KAscade) is a program for detailed simulation of extensive air showers initiated by high energy cosmic ray particles. Protons, light nuclei up to iron, photons, and many other particles may be treated as primaries.
The particles are tracked through the atmosphere until they undergo reactions with the air nuclei or - in the case of instable secondaries - decay. The hadronic interactions at high energies may be described by several reaction models alternatively:The VENUS, QGSJET, and DPMJET models are based on the Gribov-Regge theory, while SIBYLL is a minijet model. The neXus model extends far above a simple combination of QGSJET and VENUS routines. The most recent EPOS model is based on the neXus framework but with important improvements concerning hard interactions and nuclear and high-density effect. HDPM is inspired by findings of the Dual Parton Model and tries to reproduce relevant kinematical distributions being measured at colliders.
Hadronic interactions at lower energies are described either by the GHEISHA interaction routines, by a link to FLUKA, or by the microscopic UrQMD model. In particle decays all decay branches down to the 1 % level are taken into account. For electromagnetic interactions a tailor made version of the shower program EGS4 or the analytical NKG formulas may be used. Options for the generation of Cherenkov radiation and neutrinos exist. The radio emission of showers may be treated by a link with the CoREAS (Corsika-based Radio Emission from Air Showers) code.

https://www.iap.kit.edu/corsika/

### CRPropa

CRPropa is a publicly available simulation framework to study the propagation of high- and ultra-high-energy particles (cosmic-ray nuclei, photons, electrons, and neutrinos) in the Galactic and extragalactic spaces. Photonuclear interactions between cosmic rays and background photons are implemented, including processes such as Bethe-Heitler pair production, photopion production, and photodisintegration of nuclei. The code also treats electromagnetic interactions involving photons and electrons, such as pair production and inverse Compton scattering as well as their higher-order counterparts, double and triplet pair production. In addition to the one-dimensional treatment of interactions, CRPropa can compute trajectories of charged relativistic particles propagating through various magnetic-field configurations, including turbulent magnetic fields, several models of the Galactic magnetic field, as well as custom magnetic-field grids. A low-energy extension which solves the transport equation using stochastic differential equations is available, being suitable for Galactic cosmic-ray propagation. CRPropa is written in C++ with shared-memory parallelisation and Python steering. It features a modular structure, which enables users to easily add their own modules for their specific types of studies.

https://crpropa.desy.de/

***

## Gamma-ray Astronomy

### Fermipy

Fermipy is a python package that facilitates analysis of data from the Large Area Telescope (LAT) with the Fermi Science Tools. For more information about the Fermi mission and the LAT instrument please refer to the Fermi Science Support Center (https://fermi.gsfc.nasa.gov/ssc/). The Fermipy package is built on the pyLikelihood interface of the Fermi Science Tools and provides a set of high-level tools for performing common analysis tasks:

Data and model preparation with the gt-tools (gtselect, gtmktime, etc.).
Extracting a spectral energy distribution (SED) of a source.
Generating TS and residual maps for a region of interest.
Finding new source candidates.
Localizing a source or fitting its spatial extension.

Fermipy uses a configuration-file driven workflow in which the analysis parameters (data selection, IRFs, and ROI model) are defined in a YAML configuration file. Analysis is executed through a python script that calls the methods of GTAnalysis to perform different analysis operations.

https://github.com/fermiPy/fermipy

### ctapipe

ctapipe is a python based low-level data processing pipeline software for CTA (the Cherenkov Telescope Array). It provides standard tools for Cherenkov Telescope image analysis and data processing.

https://github.com/cta-observatory/ctapipe

### gammapy

Gammapy is an open-source Python package for gamma-ray astronomy built on [Numpy](https://numpy.org/), [Scipy](https://scipy.org/) and [Astropy](https://www.astropy.org/).
It is used as core library for the Science Analysis tools of the [Cherenkov Telescope Array (CTA)](https://www.cta-observatory.org/), recommended by the [H.E.S.S.](https://www.mpi-hd.mpg.de/hfm/HESS/) collaboration to be used for Science publications, and is already widely used in the analysis of existing gamma-ray instruments, such as [MAGIC](https://magic.mpp.mpg.de/), [VERITAS](https://veritas.sao.arizona.edu/) and [HAWC](https://www.hawc-observatory.org/). 

Project Web site: https://gammapy.org/, GitHub repo: https://github.com/gammapy/gammapy

### gammaALPs

gammaALPS is a python package to calculate the conversion probability between photons and axions / axion-like particles in various astrophysical magnetic fields.

https://github.com/me-manu/gammaALPs

***

## Radiative Transport

### TARDIS

TARDIS is a tool that creates synthetic observations (spectra) for exploding stars (supernovae).

https://github.com/tardis-sn/tardis

***

## N-body Simulation

### AMUSE

A Python framework to combine existing astrophysical simulation codes in numerical experiments. With AMUSE you can simulate objects such as star clusters, proto-planetary disks and galaxies.

https://amusecode.github.io (preferred citation method here: https://amusecode.github.io/copyright)

### Arepo

Arepo is a massively parallel gravity and magnetohydrodynamics code for astrophysics, designed for problems of large dynamic range. It employs a finite-volume approach to discretize the equations of hydrodynamics on a moving Voronoi mesh, and a tree-particle-mesh method for gravitational interactions. Arepo is originally optimized for cosmological simulations of structure formation, but has also been used in many other applications in astrophysics.

https://arepo-code.org

### GADGET-2

GADGET is a freely available code for cosmological N-body/SPH simulations on massively parallel computers with distributed memory. GADGET uses an explicit communication model that is implemented with the standardized MPI communication interface. The code can be run on essentially all supercomputer systems presently in use, including clusters of workstations or individual PCs.

https://wwwmpa.mpa-garching.mpg.de/gadget/ (preferred citation method here: https://wwwmpa.mpa-garching.mpg.de/gadget/right.html#License)

### GADGET-4

A parallel cosmological N-body and SPH code meant for simulations of cosmic structure formation and calculations relevant for galaxy evolution and galactic dynamics.

http://gitlab.mpcdf.mpg.de/vrs/gadget4 (preferred citation method here: https://wwwmpa.mpa-garching.mpg.de/gadget4/)

### REBOUND

REBOUND is a multi-purpose N-body integrator written in C99 and comes with an easy-to-use python front-end. It is frequently used to simulate the orbital motion of stars, planets, and moons as well as the collisional dynamics of planetary rings. REBOUND offers several different integrators, including WHFast, IAS15, Mercurius, and SEI. 

https://github.com/hannorein/rebound

***

## MHD

Note: some of the entries listed can also be classified as N-body. Keep this in mind for the future or if you are browsing.-

### Athena++

Athena++ is a complete re-write of the Athena astrophysical magnetohydrodynamics (MHD) code in C++. Compared to earlier versions, the Athena++ code has (1) much more flexible coordinate and grid options including adaptive mesh refinement (AMR), (2) new physics including general relativity, (3) significantly improved performance and scalability, and (4) improved source code clarity and modularity.

https://www.athena-astro.app

### CASTRO

Castro is an adaptive-mesh compressible radiation / MHD / hydrodynamics code for astrophysical flows. Castro supports a general equation of state, full Poisson gravity, and reactive flows, and is parallelized with MPI + OpenMP for CPUs and MPI + CUDA for GPUs.

https://amrex-astro.github.io/Castro

### Dedalus

Dedalus solves differential equations using spectral methods. It's open-source, written in Python, and MPI-parallelized.
Dedalus developed and used to study fluid dynamics, but it is designed to solve initial-value, boundary-value, and eigenvalue problems involving nearly arbitrary equations sets.

https://dedalus-project.org

### Enzo

Enzo is a community-developed adaptive mesh refinement simulation code, designed for rich, multi-physics hydrodynamic astrophysical calculations.

https://enzo-project.org

### MAESTROeX

MAESTROeX solves the equations of low Mach number hydrodynamics for stratified atmospheres/full spherical stars with a general equation of state, and nuclear reaction networks in an adaptive-grid finite-volume framework. It includes reactions and thermal diffusion and can be used on anything from a single core to 100,000s of processor cores with MPI + OpenMP or 1,000s of GPUs.

https://amrex-astro.github.io/MAESTROeX

### MagIC

MagIC is a numerical code that can simulate fluid dynamics in a spherical shell. MagIC solves for the Navier-Stokes equation including Coriolis force, optionally coupled with an induction equation for Magneto-Hydro Dynamics (MHD), a temperature (or entropy) equation and an equation for chemical composition under both the anelastic and the Boussinesq approximations.

https://magic-sph.github.io

### Nyx

Nyx is an adaptive mesh, massively-parallel, cosmological simulation code that solves equations of compressible hydrodynamics flow describing the evolution of baryonic gas coupled with an N-body treatment of the dark matter in an expanding universe.

https://amrex-astro.github.io/Nyx

### Pencil

The Pencil Code is a high-order finite-difference code for compressible hydrodynamic flows with magnetic fields. It is highly modular and can easily be adapted to different types of problems. The code runs efficiently under MPI on massively parallel shared- or distributed-memory computers.

http://pencil-code.nordita.org

### Rayleigh

Rayleigh is a 3-D convection code designed for the study of dynamo behavior in spherical geometry. It evolves the incompressible and anelastic MHD equations in spherical geometry using a pseudo-spectral approach. Rayleigh employs spherical harmonics in the horizontal direction and Chebyshev polynomials in the radial direction.

https://github.com/geodynamics/Rayleigh

### PLUTO

PLUTO is a freely-distributed software for the numerical solution of mixed hyperbolic/parabolic systems of partial differential equations (conservation laws) targeting high Mach number flows in astrophysical fluid dynamics. The code is designed with a modular and flexible structure whereby different numerical algorithms can be separately combined to solve systems of conservation laws using the finite volume or finite difference approach based on Godunov-type schemes.

[PLUTO Website](http://plutocode.ph.unito.it)

[PLUTO Documentation](http://plutocode.ph.unito.it/documentation.html)

Papers: [1](https://arxiv.org/abs/1110.0740), [2](https://arxiv.org/abs/astro-ph/0701854)

***

## Numerical Relativity

### CosmoGRaPH

CosmoGRaPH is a c++ code designed to explore cosmological problems in a fully general relativistic setting. The code implements various methods for numerically solving the Einstein field equations, and includes as matter sources an N-body/collisionless phase-space sheet solver, analysis tools useful for 


such as raytracing, and full AMR capabilities via SAMRAI.

https://cwru-pat.github.io/cosmograph/

### Einstein Toolkit

The Einstein Toolkit is a community-driven software platform of core computational tools to advance and support research in relativistic astrophysics and gravitational physics. Currently, a large portion of the toolkit is made up by over 270 Cactus components ("thorns") for computational relativity along with associated tools for simulation management and visualization. This includes vacuum spacetime solvers and relativistic hydrodynamics solvers, along with components for initial data, analysis and computational infrastructure.

https://www.einsteintoolkit.org

### gevolution

gevolution is an open-source code for exploring relativistic effects in late-Universe cosmology. It is based on a weak field expansion of General Relativity and calculates all six metric degrees of freedom in Poisson gauge. The code is also capable of raytracing to compute relativistic observables.

https://github.com/gevolution-code

### GRChombo

GRChombo is a new open-source code for numerical relativity simulations. It is developed and maintained by a collaboration of numerical relativists with a wide range of research interests, from early universe cosmology to astrophysics and mathematical general relativity, and has been used in many papers since its first release in 2015.

https://www.grchombo.org/

### GR1D

GR1D is an open-source spherically-symmetric general-relativistic (GR) hydrodynamics code. It is based on the Eulerian formulation of GR hydrodynamics (GRHD) put forth by Romero-Ibanez-Gourgoulhon and employs radial-gauge, polar-slicing coordinates in which the 3+1 equations simplify substantially.

http://www.gr1dcode.org/

### NRPy+

The NRPy+/SENR software packages implement a generalized approach for solving Einstein's equations of general relativity in various coordinate systems. The project aims to reduce the cost in memory of numerical relativity black hole and neutron star binary simulations by ~100x, through adoption of numerical grids that fully exploit near-symmetries in these systems. With this cost savings, black hole binary merger simulations can be performed entirely on a consumer-grade desktop (or laptop) computer.

http://astro.phys.wvu.edu/bhathome/nrpy.html

### SpECTRE

SpECTRE is an open-source code for multi-scale, multi-physics problems in astrophysics and gravitational physics. In the future, we hope that it can be applied to problems across discipline boundaries in fluid dynamics, geoscience, plasma physics, nuclear physics, and engineering. It runs at petascale and is designed for future exascale computers.

SpECTRE is being developed in support of our collaborative Simulating eXtreme Spacetimes (SXS) research program into the multi-messenger astrophysics of neutron star mergers, core-collapse supernovae, and gamma-ray bursts.

* Web: https://spectre-code.org/
* Repo: https://github.com/sxs-collaboration/spectre
* Preferred citation method: https://doi.org/10.5281/zenodo.4290404

***

## Primordial Black Holes

### BlackHawk

BlackHawk is a public C program for calculating the Hawking evaporation spectra of any black hole distribution. This program enables the users to compute the primary and secondary spectra of stable or long-lived particles generated by Hawking radiation of the distribution of black holes, and to study their evolution in time.

https://blackhawk.hepforge.org/, https://arxiv.org/abs/1905.04268 

### PBHbounds

A collection of bounds on primordial black holes (PBHs) and code for plotting them.

https://github.com/bradkav/PBHbounds (archived at http://doi.org/10.5281/zenodo.3538998)

### SPriBHoS

Code using pseudo-spectral methods to perform numerical simulations of spherically symmetric black hole formation on a Friedman-Robertson-Walker universe.

https://sites.google.com/fqa.ub.edu/albertescriva/home, https://arxiv.org/abs/1907.13065 

***

## Stellar Modelling

### ESTER (Evolution STEllaire en Rotation)

The ambition of this project is to set out a two-dimensional stellar evolution code, which fully takes into account the effects of rotation, at any rate and in a self-consistent way.

* Web: http://ester-project.github.io/ester/
* Repo: https://github.com/ester-project/ester

### GYRE

GYRE is a stellar oscillation code. Given an input stellar model, GYRE calculates the eigenfrequencies and eigenfunctions for the normal oscillation modes of the model. These data can be put to a variety of uses; the most common is to compare them against observed oscillation frequencies of a star, allowing constraints on the star's fundamental parameters (mass, radius, etc.) to be established  the discipline of asteroseismology.

https://github.com/rhdtownsend/gyre

### MESA

Modules for Experiments in Stellar Astrophysics (MESA) is a suite of open source, robust, efficient, thread-safe libraries for a wide range of applications in computational stellar astrophysics. A one-dimensional stellar evolution module, MESAstar, combines many of the numerical and physics modules for simulations of a wide range of stellar evolution scenarios ranging from very low mass to massive stars, including advanced evolutionary phases.

https://docs.mesastar.org

## Modelling of Active Galactic Nuclei

### agnpy

agnpy is an astropy affiliated package designed to model the radiative processes of Active Galactic Nuclei with python, focusing on the numerical computation of the photon spectra produced by leptonic radiative processes in jetted Active Galactic Nuclei (AGN).

https://github.com/cosimoNigro/agnpy

***

## Multi-Messenger Analysis 

### 3ML

The Multi-Mission Maximum Likelihood framework (3ML) provides a common high-level interface and model definition, which allows for an easy, coherent and intuitive modeling of sources using all the available data, no matter their origin. At the same time, thanks to its architecture based on plug-ins, 3ML uses under the hood the official software of each instrument, the only one certified and maintained by the collaboration which built the instrument itself. This guarantees that 3ML is always using the best possible methodology to deal with the data of each instrument.

http://threeml.readthedocs.io/

***

## Population Synthesis

### popsynth

This framework provides an abstract way to generate populations from various luminosity functions and redshift distributions. Additionally, auxiliary quantities can be sampled and stored.

Populations can be saved and restored via an HDF5 files for later use.

Note that this is not Synth Pop. If you were expecting that… I suggest you check out Depeche Mode.

https://popsynth.readthedocs.io/

***
## Neutron Stars Equations of State

### CompOSE

The online service CompOSE (CompStar Online Supernovae Equations of State) provides data tables for different state of the art equations of state (EoS) ready for further usage in astrophysical applications, nuclear physics and beyond. 

https://compose.obspm.fr/

***

## Pulsar Timing

### Tempo

Tempo is an older pulsar timing package, however it is still considered the de-facto standard in many institutes. The code has been stable for some time and so provides a stability compared to more recent code that is more often being developed and improved.

The most recent version can be downloaded from:
http://tempo.sourceforge.net/

### PINT

PINT is a project to develop a pulsar timing solution based on python and modern libraries. It is still in active development, but it is in production use by the NANOGrav collaboration and it has been demonstrated produce residuals from most "normal" timing models that agree with Tempo and Tempo2 to within ~10 nanoseconds. It can be used within python scripts or notebooks, and there are several command line tools that come with it.

https://github.com/nanograv/PINT
