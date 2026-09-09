# Code for PhD Thesis 

This repository contains code accompanying the PhD thesis
**Mixed-precision preconditioned Jacobi algorithms for accurate eigenvalue and singular value decompositions**.

This repository combines code from three repositories maintained by the author:
- [`Code_htwz25`](https://github.com/zhengbo0503/Code_htwz25):
   computing accurate eigenvalues.
- [`Code_twz26a`](https://github.com/zhengbo0503/Code_twz26a):
   computing accurate singular values.
- [`Code_twz26b`](https://github.com/zhengbo0503/Code_twz26b):
   computing accurate eigenvectors and singular vectors.

## Structure 
The three main folders correspond to 
Chapters 3, 4 and 5 of the thesis.

## Manuals
For each chapter folder, 
the `README.md` file contains a brief description 
of how the thesis code relates to the code accompanying the original papers
- [chapter 3 readme](ch3_preconditioning/README.md)
- [chapter 4 readme](ch4_eigenvalue_eigenvector/README.md)
- [chapter 5 readme](ch5_singular_value_singular_vector/README.md)

## Tested environment 
The code was tested in the following environment

|                                                                                     |                                   |
| ----------------------------------------------------------------------------------- | --------------------------------- |
| OS                                                                                  | macOS 26.6.2 (build 25G83), arm64 |
| CPU                                                                                 | Apple M3 Pro                      |
| GPU                                                                                 | Apple M3 Pro                      |
| Memory                                                                              | 36 GiB                            |
| MATLAB Version                                                                      | MATLAB R2026a                     |
| [Advanpix](https://www.advanpix.com/)                                               | Version 5.4.8 Build 16251         |
| [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS/releases#release-v0.3.29)        | 0.3.29 dev                        |
| [Julia](https://julialang.org/)                                                     | 1.11.6                            |
| [Julia Quadmath](https://github.com/JuliaMath/Quadmath.jl/releases#release-v0.5.13) | 0.5.13                            |

## Generative AI Disclosure 
The original research code was written by the author without using any generative AI tools.
The author has used Codex with [GPT-6 Astra](https://openai.com/index/gpt-6-astra/) as an assistant to:
- extend the original real-arithmetic code to complex arithmetic,
- check for bugs and inconsistencies between the code and the thesis,
- write supplementary scripts for running the experiments and performing smoke tests,
- merge the author’s three repositories into one, and
- write README files for the chapter folders.

The author takes full responsibility for the code and the results reported in the thesis.
