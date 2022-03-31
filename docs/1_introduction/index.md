---
layout: default
title: Introduction
nav_order: 2
math: mathjax3
---

<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>

# Introduction

Tailed phages are among the most abundant and diverse organisms in the ecosystem. They regulate bacterial communities and
thus play an essential role in systems from the human gut biome {%cite gut2019 %} to the global carbon cycle {%cite carbon2017%}. The genome lengths
observed in tailed phages, ranging from 5 to 550 kilobase pairs (kbp){%cite kb2020%}{%cite diana2022%}, are key to enabling their diverse roles across
ecosystems {cite}. This genetic diversity necessitates diversity in the size of the protein shells, called capsids, that 
contain these genomes{%cite diana2022 %}. Indeed, capsids of tailed phages range from 30 to 180 nanometers (nm) in diameter to contain 
these genomes{%cite hk2019 %}. Despite the diversity in size, tailed phages use only variants of a single type of protein fold, the 
HK97 fold{%cite hk2019 %} {%cite hk2015 %}. The ubiquity of this fold suggests that its versatility helped tailed phages evolve to fill
many different roles{%cite hk2019 %}.

The majority (80-90%) of tailed phages form icosahedrally symmetric capsids {%cite phages2020%}. A highly symmetric architecture is favorable 
for viruses because it reduces the genome space required to encode capsid proteins{%cite Twarock2019%}. Icosahedral symmetry is the most common 
type of symmetry seen in viral capsids. However, tailed phages stand out for how they build capsids with sizes and architectures 
not typically seen in most families of icosahedral virus{%cite phages2020 %}. <span style="color:red">Explain MCPs</span>

Generalized theories of icosahedral capsid geometry describe several possible lattices capsids could adopt. Each 2D lattice 
represents a different way that Major Capsid Proteins can form into multimeric groups called capsomers, with some tilings 
incorporating additional minor Capsid Proteins. The triangulation number, or T-number, determines the total number of proteins, 
and thus the total number of capsomers, used to build a capsid. The allowable T-numbers depend on the distance between vertices 
of the icosahedron represented in hexagonal lattice coordinates h and k. A capsid requires $$60T_0$$ total MCPs to build a 
full capsid.

<span style="color:red">I'm forgoing discussing the T-numbers differing by alpha with each lattice because the *T-number* as a measure of size/volume isn't important for my topic. I just need the base T# as a way to predict the number of capsomers</span>.

$$
\begin{equation}
    T_0 = h^2 + hk + k^2
\end{equation}
$$

The number of total capsomers in a capsid depends on the lattice as well as the T-number.

$$
\begin{equation}
    N_H = n_p + n_h = 12 + 10(T-1) = 10T + 2
\end{equation}
$$

$$
\begin{equation}
    N_T = 20T
\end{equation}
$$

$$
\begin{equation}
    N_D = 30T
\end{equation}
$$

$$
\begin{equation}
    N_{trihex} = n_p + n_h + n_t = 10T + 2 + 20T = 30T + 2
\end{equation}
$$

### <span style="color:red">[I will place a figure here demonstrating the geometry of the architectures]</span>


{% bibliography --cited %}