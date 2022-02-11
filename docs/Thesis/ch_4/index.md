---
layout: default
title: Quasi Rigid Domain Decompostion Of Viral Capsids
nav_order: 4
parent: Thesis
output: pdf_document
bibliography: test.bib
---


# Chapter 4: Quasi Rigid Domain Decomposition Of Viral Capsids

![myimg](img.png)

[@a]

## 4.1: Augmenting The Classification Scheme
The ambiguities of the classification scheme described in chapter 2 gives us the opportunity to include additional
assumptions that encode structural properties of the capsid. Here we propose that the classification scheme be 
augmented with the following assumption.

*Assumption: Capsomers correspond to quasi rigid domains of a viral capsid*

Here we use the same definition of a quasi rigid domain as described in [ref]. A rigid structures is a structure in
which the distances between elements of the structure are fixed over time and under transformations in space. A quasi-rigid structure
is thus a structure where the fluctuations between elements of the structure are minimized. We calculate the pairwise
distance fluctuation of a structure in the following manner.

$$
\begin{equation}
    f^{2}_{ij} = Var(d^{2}_{ij}) = \langle d^{2}_{ij} \rangle - \langle d_{ij} \rangle ^{2}
\end{equation}
$$

A quasi rigid domain of a protein structure is a domain of the protein which satisfied our definition of a quasi rigid 
structure.

## 4.2 Normal Mode Analysis
Calculating the distance fluctuations between elements of a protein structure is not a simple task. The most direct, and
computationally expensive method would be using molecular dynamics, averaging the fluctuations between elements over the
trajectories. While accurate and intuitive this method is computationally infeasible with even small capsids containing
far too many elements. We can examine our requirements to select an alternate method. We are interested in the large
scale dynamics of the capsid, we are interested in the dynamics of the capsid near equilibrium, and we want to use
techniques that are scalable to large capsids.
The technique we will make use of is Normal Mode Analysis
![Alt Text](1a34.gif)

Normal Mode Analysis is a technique aimed towards describing the equilibrium dynamics of a physical system. It aims to
approximate the way the system fluctuates around the equilibrium by assuming oscillatory behavior and considering only
a subset of the normal modes of the system. The assumptions necessary to allow this technique are that the system has a 
specific equilibrium configuration and that all particles in the system interact under a simple harmonic potential. This
assumption is taken to be accurate only locally. The further the system strays from the equilibrium the less accurate this
technique will be. 

The method disregards any specific interactions and constraints in the system. As a result it describes only macroscopic
motions of the system and will fail to represent complex bonds within a protein. This means it is best paired with a model
that also doesn't concern itself with microscopic forces and constraints. It also suggests that the technique is best 
applied to systems large enough that collective motions are dominant.

The mathematical formulation of NMA begins by examining a taylor series of the potential about the equilibrium.

$$
\begin{equation}
    V(\vec{x}) = V(\vec{x_0}) + \sum_{i}(\frac{\partial V}{\partial \vec{x_0}}\delta \vec{x})
\end{equation}
$$

## 4.3 Elastic Network Models
Normal Mode Analysis is a model independent technique, and the choice of model is a separate and important decision.
Among coarse grained models aimed at describing large scale dynamics the most popular options are Elastic Network 
Models (ENM). 

## 4.4 Spectral Clustering

## 4.5 Results
# References

