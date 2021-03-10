# HEP-ASTRO-COSMO

**This is a community effort to collect all HEP/ASTRO/COSMO open source packages/libraries/tools in one place. Everyone welcome to contribute!** :shipit:

*Original idea: [GF Bertone](https://twitter.com/gfbertone), [Twitter discussuion can be found here](https://twitter.com/malcfairbairn/status/1369178173884235776?s=20). First resource made by [Suchita Kulkarni](https://twitter.com/suchi_kulkarni), [this Gdoc](https://docs.google.com/document/d/1yDp4EfxR5ivlDhice2iQ3jeOJNN0GRKH_-ln5y3R7UY/edit).*

**Package/library/tool descriptions copied from respective pages.**


***

## Table of contents
- [HEP](#hep)
  * [Tools](#tools)
  * [Event Generators](#event-generators)
  * [Event Analysers](#event-analysers)
  * [Global Fitters](#global-fitters)
  * [Spectrum Generators](#spectrum-generators)
  * [Direct Detection](#direct-detection)
  * [Feyn Family](#feyn-family)
  * [Statistics](#statistics)
- [COSMO](#cosmo)
  * [Cosmic Microwave Background](#cosmic-microwave-background)
  * [General Cosmology](#general-cosmology)
- [ASTRO](#astro)
  * [Cosmic Rays](#cosmic-rays)
  * [N-body Simulations](#n-body-simulation)
  * [Primordial Black Holes](#primordial-black-holes)

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

### LHAPDF

LHAPDF is a general purpose C++ interpolator, used for evaluating PDFs from discretised data files. Previous versions of LHAPDF were written in Fortran 77/90 and are documented at http://lhapdf.hepforge.org/lhapdf5/. LHAPDF6 vastly reduces the memory overhead of the Fortran LHAPDF (from gigabytes to megabytes!), entirely removes restrictions on numbers of concurrent PDFs, allows access to single PDF members without needing to load whole sets, and separates a new standardised PDF data format from the code library so that new PDF sets may be created and released easier and faster. The C++ LHAPDF6 also permits arbitrary parton contents via the standard PDG ID code scheme, is computationally more efficient (particularly if only one or two flavours are required at each phase space point, as in PDF reweighting), and uses a flexible metadata system which fixes many fundamental metadata and concurrency bugs in LHAPDF5. Compatibility routines are provided as standard for existing C++ and Fortran codes using the LHAPDF5 and PDFLIB legacy interfaces, so you can keep using your existing codes. But the new interface is much more powerful and pleasant to work with, so we think you'll want to switch once you've used it! LHAPDF6 is documented in more detail in http://arxiv.org/abs/1412.7420

https://lhapdf.hepforge.org


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

### SARAH

SARAH is a Mathematica package for building and analyzing SUSY and non-SUSY models. It calculates all vertices, mass matrices, tadpoles equations, one-loop corrections for tadpoles and self-energies, and two-loop RGEs for a given model. SARAH writes model files for FeynArts, CalcHep/CompHep, which can also be used for dark matter studies using MicrOmegas, the UFO format which is supported by MadGraph 5 and for WHIZARD and OMEGA. 
SARAH was also the first available spectrum-generator-generator: based on derived analytical expressions it creates source code for SPheno. It is therefore possible to implement new models in SPheno without the need to write any Fortran code by hand. The output for Vevacious can be used to check for the global minimum for a given model and parameter point.
Running SARAH is fast, it already includes a long list of SUSY and non-SUSY models, and the implementation of new models is efficient and straightforward.

https://sarah.hepforge.org/

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

***

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


***
# COSMO

## Cosmic Microwave Background

### CAMB

CAMB is a Python and Fortran code for computing CMB, CMB lensing, lensing, galaxy count and dark-age 21cm power spectra, transfer functions and matter power spectra, and background cosmological functions.

https://camb.info (preferred citation method here: https://cosmologist.info/cosmomc/cosmomc.bib)

### CLASS

The purpose of CLASS is to simulate the evolution of linear perturbations in the universe and to compute CMB and large scale structure observables. Its name also comes from the fact that it is written in object-oriented style mimicking the notion of class. Classes are a wonderfull programming feature available e.g. in C++ and python, but these languages are known to be less vectorizable/parallelizable than plain C (or Fortran), and hence potentially slower. For CLASS we choose to use plain C for high performances, while organizing the code in a few modules that reproduce the architecture and philosophy of C++ classes, for optimal readability and modularity.

http://class-code.net

### CosmoMC

CosmoMC is a Fortran 2008 Markov-Chain Monte-Carlo (MCMC) engine for exploring cosmological parameter space, together with Fortran and python code for analysing Monte-Carlo samples and importance sampling (plus a suite of scripts for building grids of runs, plotting and presenting results). The code does brute force (but accurate) theoretical matter power spectrum and Cl calculations with CAMB. See the [original paper](https://arxiv.org/abs/astro-ph/0205436) for an introduction and descriptions, and [up-to-date sampling algorithm](https://arxiv.org/abs/1304.4473) details. It can also be compiled as a generic sampler without using any cosmology codes.

https://cosmologist.info/cosmomc

### Monte Python

Monte Python is a Monte Carlo code for Cosmological Parameter extraction. It contains likelihood codes of most recent experiments, and interfaces with the Boltzmann code class for computing the cosmological observables.

http://baudren.github.io/montepython.html

***

## General Cosmology

### CosmoloPy


A cosmology package for Python.

http://roban.github.io/CosmoloPy/

***
# ASTRO

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
