---
layout: default
title: Quasi Rigid Domain Decompostion Of Viral Capsids
nav_order: 4
parent: Thesis

---


# Chapter 4: Quasi Rigid Domain Decompostion Of Viral Capsids

![myimg](img.png)

## 4.1: Augmenting The Classification Scheme
The ambiguities of the classification scheme described in chapter 2 gives us the opportunity to include additional
assumptions that encode structural information about the capsid. Here we propose that the classification scheme be 
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

![Alt Text](https://github.com/luquelab/cbrown_thesis_athena/blob/gh-pages/docs/Thesis/ch_4/1a34.gif)


## 4.3 Elastic Network Models
Normal Mode Analysis is a model independent technique, and the choice of model is a separate and important decision.
Among coarse grained models aimed at describing large scale dynamics the most popular options are Elastic Network 
Models (ENM). 

## 4.4 Spectral Clustering

## 4.5 Results

