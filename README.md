# MR Thermometry algorithm code: 
# harmonic initialized model-based multi-echo (HIMM)
***
## Multi-echo MR Thermometry in the upper leg at 7T using near-harmonic 2D reconstruction for initialization

#### M.W.I. Kikken, B.R. Steensma, C.A.T. van den Berg, and A.J.E. Raaijmakers
#### Computational Imaging Group for MRI Therapy & Diagnostics, University Medical Center Utrecht
***
This code repository implements an algorithm for real-time-compatible water/fat separated MR thermometry in aqueous and mixed water/fat tissues. An article detailing the work has been submitted for publication to a peer reviewed journal. Use of this algorithm is allowed given proper citation of the original work and in accordance with the license file.
The algorithm provided in this repository is based on the Multi-echo MR thermometry algorithm using iterative separation of baseline water and fat images with l1 regularization. This method was proposed by Poorman et al. (https://github.com/poormanme/waterFatSeparated_MRThermometry) and published in Magnetic Resonance in Medicine [1].
***
The current algorithm mainly focuses on the measurements of RF heating. RF heating distributions are more diffisive, so the principle of l1 regularization does not work anymore (as that requires a localized heating hotspot, with very different spatial distribution from the drift field). The l1 regularization was removed from the code. In order to still separate the drift fields from the temprature increases, we propose to initialize the drift fields using near-harmonic 2D reconstruction (proposed by Salomir et al.). Near-harmonic 2D reconstuction used the fat layer (not susceptible to temperature-induced phase changes) that surrounds the anatomy to estimate the drift field in the whole region of interest. This drift field is decomposed into spherical harmonics and the corresponding coefficients are iteratively optimized together with the spatial temperature map and the coefficients of the library.
***
* __RunHIMM.m__: Example script to run algorithm on simulated data
* __thermo_hybrid_waterfat__: Main function script
***
[1] Poorman, M. E., Braškutė, I., Bartels, L. W., & Grissom, W. A. (2019). Multi‐echo MR thermometry using iterative separation of baseline water and fat images. Magnetic resonance in medicine, 81(4), 2385-2398.
