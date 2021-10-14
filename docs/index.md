---
layout: default
title: Home
nav_order: 1
has_children: True
---

# Home

My name is Colin Brown, an MS student in Physics, and this repository contains a record of my contributions to the lab and
the development of my thesis. I primarily analyze the structural properties of viral capsids with coarse grained models.
The majority of the content here will be a part of my Masters Thesis in Physics.

## 1. Thesis Outline

### 1.1 Introduction To Icosahedral Capsids

### 1.2 Geometry Of Icosahedral Capsids

### 1.3 Coarse Grained Models Of Capsid Geometries

### 1.4 Quasi Rigid Domain Decompostion Of Viral Capsids

#### 1.4.1 Geometry Of Quasi Rigid Subunits

#### 1.4.2 Normal Mode Analysis

#### 1.4.3 Elastic Network Models

#### 1.4.4 Spectral Clustering

#### 1.4.5 Implementation On HPC Cluster

#### 1.4.6 Results

### 1.5 Discussion & Conclusion

#### 1.5.1-3 Results Of Each Chapter

#### 1.5.4 Synthesis Of Results

#### 1.5.5 Computational Improvements

#### 1.5.6 Future Research

The associated file `index.md` contains a YAML front matter to indicate the layout, title, and navigation options. The repo's website is based on Jekyll's theme [Just-the-Docs](https://pmarsceill.github.io/just-the-docs/). Explore their [documentation]([Just-the-Docs](https://pmarsceill.github.io/just-the-docs/)) and associated [GitHub repo](https://github.com/pmarsceill/just-the-docs) to adapt your project's website to your needs.
---

The `head.html` file from the original template has been modified in `/docs/_includes` to include MathJax, so you can write math using latex format. Here are some examples:

Math equation using MathJax: $$5+5$$

$$
\begin{align*}
  & \phi(x,y) = \phi \left(\sum_{i=1}^n x_ie_i, \sum_{j=1}^n y_je_j \right)
  = \sum_{i=1}^n \sum_{j=1}^n x_i y_j \phi(e_i, e_j) = \\
  & (x_1, \ldots, x_n) \left( \begin{array}{ccc}
      \phi(e_1, e_1) & \cdots & \phi(e_1, e_n) \\
      \vdots & \ddots & \vdots \\
      \phi(e_n, e_1) & \cdots & \phi(e_n, e_n)
    \end{array} \right)
  \left( \begin{array}{c}
      y_1 \\
      \vdots \\
      y_n
    \end{array} \right)
\end{align*}
$$

\( 5+5 \)

When $$a \ne 0$$, there are two solutions to $$ax^2 + bx + c = 0$$ and they are
\begin{equation}
  x = \frac{-b \pm \sqrt{b^2-4ac}}{2a}
\end{equation}

May Athena be with you.
