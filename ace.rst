:title: Ace Exam
:author: Bas Rustenburg
:description: Admission to candidacy exam
:keywords: exam, ace, phd
:css: ace.css

.. |lt_theta| image:: images/colored_theta.png
.. |lt_prior| image:: images/colored_prior.png
.. |lt_posterior| image:: images/colored_posterior.png
.. |lt_likelihood| image:: images/colored_likelihood.png
.. |lt_model| image:: images/colored_model.png
.. |lt_Bayes| image:: images/colored_bayes_rule.png
.. |lt_dG| image:: images/colored_dG.png
.. |lt_dH| image:: images/colored_dH.png
.. |lt_H0| image:: images/colored_H0.png
.. |lt_Xs| image:: images/colored_Xs.png
.. |lt_Mc| image:: images/colored_Mc.png
.. |lt_sigma| image:: images/colored_sigma.png
.. |lt_norm| image:: images/colored_norm_n.png
.. |lt_variance| image:: images/colored_variance.png

:data-transition-duration: 500

----

Admission to candidacy exam
===========================

Bas Rustenburg
--------------

Date : tbd
..........

.. image:: images/background.png
  :align: center
  :width: 800px

----

:hovercraft-path: M1,5L20,20L40,10L60,40L80,5L100,60

::
    :hovercraft-path: M600,350 l 50,-25,a25,25 -30 0,1 50,-25 l 50,-25, a25,50 -30 0,1 50,-25 l 50,-25,a25,75 -30 0,1 50,-25 l 50,-25, a25,100 -30 0,1 50,-25 l 50,-25

Alchemical free energy calculations
===================================

Alchemical free energy calculations are a powerful computational tool for computing binding free energies, as they allow for efficient sampling of the relevant states of protein-ligand complexes.

.. image:: images/alchem.png
  :width: 700px

www.alchemistry.org

----



Specific aims
=============

In this proposal, we address three of the most significant open challenges
in the quantitative modeling of small molecule recognition by alchemical free energy calculations.

Aims
----

1. Establish a correct quantitative treatment of alchemical free energy calculations for binding of charged ligands.

2. Quantify the magnitude of protonation state effects on binding

3. Develop a framework for alchemical free energy calculations to describe weak association and cooperative ligand binding.

----

Establish a correct quantitative treatment of alchemical free energy calculations for binding of charged ligands
================================================================================================================
Aim 1.
--------


+----------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| .. image:: images/reif_oostenbrink.png | In order to apply alchemical free energy calculations to charged ligands, one needs to eliminate artifacts introduced into the calculation arising from the modeling of bulk solvent behavior using a small periodic system. |
+----------------------------------------+                                                                                                                                                                                                                              |
| Image Source: [#]_                     | Ligand interactions with:                                                                                                                                                                                                    |
|                                        |                                                                                                                                                                                                                              |
|                                        | * solvent (Blue)                                                                                                                                                                                                             |
|                                        | * receptor (Red)                                                                                                                                                                                                             |
|                                        | * self-interaction (Green)                                                                                                                                                                                                   |
+----------------------------------------+------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. [#] MM Reif and C Oostenbrink. J Comput Chem 35.3 (Nov. 2013), pp. 227–243


----



Establish a correct quantitative treatment of alchemical free energy calculations for binding of charged ligands
================================================================================================================
Aim 1.
--------

+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| .. image:: images/reif_oostenbrink.png | Bulk liquids are approximated in simulation, either by using periodic boundary conditions, or an implicit solvent.                                                                    |
|                                        | Often, to further reduce computation cost, we introduce truncated,potentials and non-Coulombic electrostatics (such as particle mesh Ewald,[PME],and reaction field [RF] potentials). |
+----------------------------------------+                                                                                                                                                                                       |
| Image Source: [#]_                     |                                                                                                                                                                                       |
|                                        |                                                                                                                                                                                       |
|                                        |                                                                                                                                                                                       |
|                                        |                                                                                                                                                                                       |
|                                        |                                                                                                                                                                                       |
+----------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

.. [#] MM Reif and C Oostenbrink. J Comput Chem 35.3 (Nov. 2013), pp. 227–243


----



Establish a correct quantitative treatment of alchemical free energy calculations for binding of charged ligands
================================================================================================================
Aim 1.
------

A number of corrections have been proposed but:
 * They have not been compared to each other
 * Quantitative correctness of these methods has not been established

Sources:
 - MM Reif and C Oostenbrink. J Comput Chem 35.3 (Nov. 2013), pp. 227–243
 - GJ Rocklin et al. J Chem Phys 139.18 (2013), p. 184103.
 - YL Lin et al.  J Chem Theory Comput 10.7 (July 2014), pp. 2690–2709.

----



Establish a correct quantitative treatment of alchemical free energy calculations for binding of charged ligands
================================================================================================================
Aim 1.
------

Subaim 1.1:  Develop an accurate approach to quantifying experimental uncertainty in ITC using Bayesian inference.
..................................................................................................................

Because we need a reliable experimental dataset in order to make a quantitative comparison

Subaim 1.2: Perform a quantitative comparison of suggested correction models to experiments to establish a correct treatment of charged ligands in alchemical free energy calculations.
.......................................................................................................................................................................................

Evaluating the charge corrections, testing an alternative (counter ions), comparing to each other and experiment

----

The host-guest model system
===========================

Aim 1
-----

We will use cucurbit-\[7\]-uril as a model system

+-----------------------------------+------------------------------------+----------------------------------------------------------------------------------+
| .. image:: images/guest11_top.png | .. image:: images/guest11_side.png | The system is useful because:                                                    |
|   :width: 200px                   |   :width: 200px                    |                                                                                  |
|                                   |                                    | * Both guest and hosts are very soluble                                          |
+-----------------------------------+------------------------------------+ * They are small, with few degrees of freedom                                    +
| .. image:: images/Kd_guest2.png                                        | * The affinities are in the range of typical protein-small molecule interactions |
|   :width: 410px                                                        |                                                                                  |
+------------------------------------------------------------------------+----------------------------------------------------------------------------------+

----



Develop an accurate approach to quantifying experimental uncertainty in ITC using Bayesian inference.
=====================================================================================================

Subaim 1.1
----------
The experimental parameters, |lt_theta| , can be estimated using Bayes rule:
|lt_Bayes| , where

  - |lt_posterior| is the posterior distribution. The probability of the parameters given the observed data. *This is what we want to know!*
  - |lt_likelihood| is the likelihood. The probability of the observed data, given a single set of parameters.
  - |lt_prior| are distributions containing prior information. We can use this to propagate errors.


We can sample from the posterior distribution by using a technique called *Markov chain Monte Carlo*.

----

Bayes rule in effect
====================

We apply |lt_Bayes| for a value centered around zero,
with prior information that it is between -1 and 1, uniformly distributed.


.. figure:: images/bayes_dist.png


----




Sampling from a posterior distribution using MCMC
=================================================


Flipping an weighted coin


.. figure:: images/distributions.png


  http://bayesianbiologist.com

----





Develop an accurate approach to quantifying experimental uncertainty in ITC using Bayesian inference.
=====================================================================================================

Subaim 1.1
----------

The ITC model structure
.......................

.. image:: images/colored_parameters.png
  
Thermodynamic parameters include
  
  - binding affinity, |lt_dG|
  - enthalpy of binding, |lt_dH|
  - mechanical heats offset, |lt_H0|
  - concentration of syringe component, |lt_Xs|
  - concentration of cell component, |lt_Mc|
  - noise parameter, |lt_sigma|

We can use prior distributions |lt_prior| to propagate error estimates in concentrations, and include previous measurements.


----

Develop an accurate approach to quantifying experimental uncertainty in ITC using Bayesian inference.
=====================================================================================================

Subaim 1.1
----------


The ITC model structure
.......................

+--------------------------------------------------------------+-----------------------------------------------+
| The likelihood model, |lt_likelihood|, is defined as         | .. image:: images/normal.png                  |
|                                                              |   :height: 350px                              |
| .. image:: images/colored_model.png                          |                                               |
+--------------------------------------------------------------+-----------------------------------------------+
| Where the observed heats are sampled from a normal distribution |lt_norm|, with a variance of |lt_variance|. |
+--------------------------------------------------------------------------------------------------------------+

----

Compare the different charge correction models
==============================================

Subaim 1.2
----------

We will consider these approaches:

* Reif and Oostenbrink use thermodynamic cycles to eliminate individual components
* Rocklin et al. use Poisson-Boltzmann calculations with exact either numerical solutions to quantify the erroneous contributions.
* Lin et al. use potential of mean force (PMF) calculations in a large simulation system, pulling the ligand away from the protein non-alchemically.
* Eliminating a pair of ions, with a net charge of 0.

We will first check if the methods produce the same quantitative estimate.
Next, we will compare to experiment, to see if they produce a quantitatively correct answer.

This is the first comparison of any of these methods on the same system!

----


Quantify the magnitude of protonation state effects on binding
==============================================================
  
Aim 2.
------

Proteins and many small-molecule drugs contain titratable moieties that can change protonation state upon binding or sample mixtures of protonation states, often in a conformation-dependent manner. Evidence exists that for the binding of imatinib to Abl kinase, pH dependent effects may contribute to the binding affinity, and preliminary data indicates that it is the same for *many other kinase inhibitors*.


+---------------------------------------+--------------------------------------------+
| .. image:: images/inhibitor-pKas.png  | .. image:: images/imatinib_image_curve.png |
|   :width: 400px                       |   :width: 400px                            |
+---------------------------------------+--------------------------------------------+


----

Quantify the magnitude of protonation state effects on binding
==============================================================
  
Aim 2.
------

Subaim 2.1: Benchmark small molecule pKa prediction tools against experimental data for kinase inhibitors.
..........................................................................................................
We need reliable pKa estimates of small molecule kinase inhibitors. We will benchmark available tools and compare to experimental data.


Subaim 2.2: Survey the kinase:inhibitor cocrystal structures for possible protonation state effects in inhibitor binding.
.........................................................................................................................
We will identify kinase-inhibitor systems that show changes in the populations of protonation states from MCCE calculations.

Subaim 2.3: Dissect the determinants and impact of protonation state effects on binding affinity through free energy calculations and ITC experiments.
......................................................................................................................................................
The systems identified will be simulated using alchemical free energy calculations, and we will perform ITC experiments on them.
 
----

Benchmark small molecule pKa prediction tools against experimental data for kinase inhibitors.
==============================================================================================

Subaim 2.1
----------

We will consider a variety of tools that are capable of predicting small molecule pKa data.

* **MoKa** generates pKa s based on atomistic descriptors, defined by the surrounding atoms. The descriptors are based on molecular interaction fields calculated using GRID for a library of 3D fragments, but can successfully be applied on 2D structures.
 
* Schrodinger’s **Jaguar** provides means of estimating pKa values using quantum mechanical methods.
 
* **Epik** uses Hammett Taft linear free energy approaches [86] for predicting pKa values.

The results provided by these tools will be compared against available pKa data from a Sirius T3 instrument. We have already measured the pKa of several of these inhibitors. For the sake of having a completely computational framework to perform these calculations, we would like to find a reliable predictor.


----

Survey the kinase:inhibitor cocrystal structures for possible protonation state effects in inhibitor binding.
=============================================================================================================

Subaim 2.2
----------

We will investigate complex structures from the protein databank, using a framework called MCCE.

.. image:: images/imatinib_sites.png
  :width: 800px
  

----

Survey the kinase:inhibitor cocrystal structures for possible protonation state effects in inhibitor binding.
=============================================================================================================

MCCE samples multiple conformations of protein side-chains and estimates the most probably protonation state.
The framework has been extended to incorporate sampling of ligands. We will keep ligand conformations fixed to those found in crystal structures.

Subaim 2.2
----------
  
.. image:: images/mcce2_sharp.png
  :width: 400px

----



That's all folks!
=================
