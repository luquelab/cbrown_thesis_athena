---
layout: default
title: Approach
nav_order: 4
---

# 1. Quasi Rigid Domain Composition Approach

## 1.1 Data Acquisition

## 1.2 Coarse Graining

## 1.3 Anisotropic Network Model

## Normal Mode Analysis

## Spectral Clustering

### Spectral Decomposition

### K-Means Clustering

## Interpretation Of Results

## Symmetry Reduction Of Matrices

Knowing the icosahedral rotation matrices ahead of time can allow us to simplify the memory requirements using icosahedral symmetry. If we have our points labeled by which asymmetric unit they belong to and their position within that asymmetric unit we can calculate the distance vector between points using the following equation.

$$
\begin{equation}
    \textbf{r}_{ij,kl} = I_k \textbf{x}_{i} - I_l \textbf{x}_{j}
\end{equation}
$$

A property of symmetry groups is that products of group actions yield another member of the group. We can use this fact to express distance vectors between members of the asymmetric unit copies in terms of the distance vectors between the base unit and all others because of the following relations.

$$
\begin{equation}
    I_k^{-1} \textbf{r}_{ij,kl} = I_k^{-1} I_k \textbf{x}_{i} - I_k^{-1} I_l \textbf{x}_{j} = \textbf{x}_{i} - I_m \textbf{x}_{j} = \textbf{r}_{ij,0m}
\end{equation}
$$

$$
\begin{equation}
    \textbf{r}_{ij,kl} = I_k \textbf{r}_{ij,0m}
\end{equation}
$$


Here $I_m$ is a solution to the following equation:

$$
\begin{equation}
    \textbf{I}_m =  \textbf{I}_k^{-1} \textbf{I}_l
\end{equation}
$$
We calculate the solutions to eq. 4 ahead of time using the group table calculated from the rotation matrices.

The Distance matrix can be reduced without having to apply a rotation. This also applies to the connectivity matrix since
it depends solely upon distance.

$$
\begin{equation}
    \textbf{D}_{ij,kl} =  \textbf{D}_{ij,0m}
\end{equation}
$$

The Hessian Matrix is calculated using the outer product of the distance vectors.



$$
\begin{equation}
    H_{ij,kl} = \frac{\textbf{$\Gamma$}_{ij,kl}}{\textbf{D}_{ij,kl}} \textbf{r}_{ij,kl} \otimes \textbf{r}_{ij,kl}
\end{equation}
$$

The outer product between two vectors expressed using equation (3) yields:

$$
\begin{equation}
    \textbf{r}_{ij,kl} \otimes \textbf{r}_{ij,kl} = I_k \textbf{r}_{ij,0m} \otimes \textbf{r}_{ij,0m} I_k^{-1}
\end{equation}
$$

Since all of the matrices in equation (6) can be reduced this way the final reduction for the hessian is:

$$
\begin{equation}
    \textbf{H}_{ij,kl} =  I_k \textbf{H}_{ij,0m} I_k^{-1}
\end{equation}
$$