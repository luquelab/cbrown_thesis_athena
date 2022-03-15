---
layout: default
title: Methods
nav_order: 3
math: mathjax3
---

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>

# 2. Methods

## 2.1 Data Acquisition

All atomic models are acquired in PDB format and are, with few exceptions, from the RCSB Protein Databank. The
majority of .pdb files for viral capsids contain an asymmetric unit and a set of icosahedral rotations that build the
full capsid. {%cite PDB101 %} These models are built by fitting folded proteins into density distributions measured using X-Ray Crystallography
or Cryo-electron Microscopy. 



| ![myimg](2e0z_pdb.png) |
|:--:| 
| *Figure 1: The PDB file 2e0z visualized in ChimeraX* |


## 2.2 The Anisotropic Network Model

Elastic Network Models (ENMs) are among the most popular models for describing large scale protein dynamics. They represent
proteins as a network of masses and springs in an equilibrium state. They require very few parameters to fully describe the system, and are
also able to be coarse grained to any level depending on computational needs. We select the Anisotropic Network Model (ANM),
the most commonly used ENM, for its ability to describe protein conformational changes in three dimensions.{% cite Bahar2010 %} 
We construct our model by coarse-graining to the level of protein residues, selecting the carbon alpha atoms as the representative
coordinates of each residue, and connect them to other residues within a cutoff distance.

| ![](2e0z_enm.png) |
|:--:| 
| *Figure 2: A representation of an Elastic Network Model using the example pdb 2e0z.* |

The overall potential of the system is thus the sum of harmonic potentials between each residue. The summation is performed
only over connected residues determined by a spring connectivity matrix.

$$
\begin{equation}
    V(\vec{x}) =  \frac{1}{2 \sum_{i|i \neq j} \Gamma_{ij} (||\vec{x}_i - \vec{x}_j|| - ||\vec{x}^0_i - \vec{x}^0_j||) }
\end{equation}
$$

Where $$\vec{x}_i$$ is the coordinate vector of residue i and $$\vec{x}_i^0$$ is the equilibrium coordinate vector for that residue.
$$\Gamma$$ is our inter-residue connectivity matrix, with each entry determined using our cutoff distance and choice of spring constant.  {% cite Eyal2006 %}
The connectivity matrix on its own represents the topology of our system, similar to a connected graph.

$$
\begin{equation}
    \Gamma_{ij} = \biggr \{
    \begin{array}{ll}
      \gamma, & R_{ij} \leq r_c \\
      0, & R_{ij} > r_c
    \end{array} 
\end{equation}
$$

To simplify the model we set the spring constant to 1 for all residues. Our cutoff distance is set to $$18Å$$, which yields
the best agreement between residue square fluctuations and experimental b-factors.




## 2.3 Normal Mode Analysis

We are interested in the large scale dynamics of the capsid near equilibrium. This prompts us to make use of a technique
called Normal Mode Analysis.

|![Alt Text](Test.gif)|
|:--:| 
| *Figure 3: An animation showing vibration along one of the normal modes* |

Normal Mode Analysis is a technique for describing the near equilibrium dynamics of a physical system. It aims to
approximate vibrations around the equilibrium by assuming harmonic potentials and considering only
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
    V(\vec{x}) = V(\vec{x^0}) + \sum_{i}\Delta x_i \frac{\partial V}{\partial x_i }\biggr|_{x=x^0}  + \sum_{i,j}\Delta x_i \Delta x_j \frac{\partial^2 V}{\partial x_i \partial x_j }\biggr|_{x=x^0} + \dots
\end{equation}
$$

The first and second terms of this expansion are zero in any equilibrium conformation. Truncating the remaining terms
gives us our second order expansion of our potential about the equilibrium.

The matrix of second derivatives of our potential around the equilibrium is called the Hessian Matrix.

$$
\begin{equation}
    H_{ij} = (\frac{\partial^2 V}{\partial x_i \partial x_j})^0
\end{equation}
$$

Our equation of motion may be written using the Hessian as follows:

$$
\begin{equation}
    \boldsymbol{M} \frac{d^2 \Delta \vec{x}}{dt^2} + \boldsymbol{H} \Delta \vec{x} = 0
\end{equation}
$$

Where the matrix M is a mass matrix, which in our case is the identity matrix and can be ignored. The normal modes of
the system are thus solutions to the following eigenvalue problem.

$$
\begin{equation}
    \boldsymbol{H} \vec{v_k} = \omega^2 \vec{v_k}
\end{equation}
$$

The Hessian of our ANM can be derived from our potential in Eq. (1). Because ANM uses three dimensional coordinates the
Hessian is a $$3Nx3N$$ blovk matrix that consists of $$NxN$$ blocks for each residue interaction. The off-diagonal blocks
have the following form.

$$
\begin{equation}
    \mathbf{H}_{ij} = \frac{\textbf{$\Gamma$}_{ij}}{R_{ij}^2} \vec{r}_{ij} \otimes \vec{r}_{ij}
\end{equation}
$$

Where $$\vec{r}_{ij}$$ is the distance vector between residues, $$R_{ij}^2$$ is the distance between residues, and
$$\otimes$$ denotes the outer product of two vectors yielding a $$3 x 3$$ block.
The diagonal blocks of our Hessian Matrix are the sum of all other blocks in that row.

$$
\begin{equation}
    \mathbf{H}_{ii} = - \sum_{i|i \neq j} \mathbf{H}_{ij}
\end{equation}
$$

From the results of NMA we can determine the pairwise fluctuations in distance between residues by first constructing the
cross correlation between the fluctuation of residues. The correlation matrix is related to the inverse of the Hessian by taking
the trace of each 3x3 submatrix.

$$
\begin{equation}
    \mathbf{C}_{ij} = \langle \Delta \mathbf{R}_i \Delta \mathbf{R}_j \rangle = tr(\mathbf{H}^{-1}_{ij})
\end{equation}
$$

The Hessian matrix is, however, singular and cannot be exactly inverted, having 6 zero eigenvalues. We instead construct
a pseudo-inverse from the eigenvectors/normal modes we calculated.

$$
\begin{equation}
    \mathbf{H}^{-1} = \sum_{k=1}^{3N - 6} \frac{1}{\omega_k^2} \vec{v_k} \otimes \vec{v_k}
\end{equation}
$$


| ![myimg](distflucts.png) |
|:--:| 
| *Figure 4: A matrix of pairwise distance fluctuations* |



## 2.4 Spectral Clustering

$$
\begin{equation}
    f^{2}_{ij} = Var(d^{2}_{ij}) = \langle d^{2}_{ij} \rangle - \langle d_{ij} \rangle ^{2}
\end{equation}
$$

*Assumption: Capsomers correspond to quasi rigid domains of a viral capsid*

Here we use the same definition of a quasi rigid domain as described in {% cite Ponzoni2015 %}. A rigid structures is a structure in
which the distances between elements of the structure are fixed over time and under transformations in space. A quasi-rigid structure
is thus a structure where the fluctuations between elements of the structure are minimized. We calculate the pairwise
distance fluctuation of a structure in the following manner.



A quasi rigid domain of a protein structure is a domain of the protein which satisfied our definition of a quasi rigid 
structure.


Now that we have determined the pairwise distance fluctuations between the residues of the capsid we need to determine
an optimal subdivision, or clustering, of the system. There exist many algorithms to identify optimal clusterings of
data. One of the most effective algorithms used when dealing with large, sparsely connected systems is Spectral Clustering.
This method requires us to first transform our measure of dissimilarity, distance fluctuations, into a measure of similarity.

$$
\begin{equation}
    S_{i,j} = e^{-D_{i,j}^2 / 2 \bar{D}^2}
\end{equation}
$$

We can use the nature of connectivity in our model to simplify our similarity matrix by setting the similarity of unconnected
residues to zero. 

### Spectral Graph Embedding

Spectral embedding is a technique based on graph theory, and requires as an input a Laplacian Matrix representing a graph.
We can transform a similarity matrix into a Laplacian matrix, specifically the Symmetric Normalized Laplacian, with the 
following identity.

$$
\begin{equation}
    \mathbf{L} = \mathbf{I} - \mathbf{D}^{-1/2} \mathbf{S} \mathbf{D}^{-1/2}
\end{equation}
$$


### Clustering Embedded Points

The eigenvectors of this graph now represent a set of points in a higher dimensional space that can be clustered
using one of many methods. We choose a 

## 2.5 Scoring & Selection

Since our methods take the number of clusters as input, we need to compare results across different numbers of clusters
and select the optimal clustering. We

| ![myimg](2e0z_32_domains.png) |
|:--:| 
| *Figure 5: Comparative plots for Pyrococcus Furiosus VLP* |

## 2.6 Classification & Visualization

The labels assigned to each residue allow us to visualize the results of the clustering in ChimeraX. We color each residue
based on its cluster label and then overlay a 3d structure with a lattice that fits the clustering most accurately.

| ![myimg](2e0z_subdivision.png) |
|:--:| 
| *Figure 6: The results of visualizing Pyrococcus Furiosus VLP* |


# References

{% bibliography --cited %}