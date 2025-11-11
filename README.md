# Online Conformal Inference with Retrospective Adjustment

This repository contains the implementation of the paper
> **"Online Conformal Inference with Retrospective Adjustment for Faster Adaptation to Distribution Shift"**  
> Jungbin Jun and Ilsang Ohn, 2025.

Preprint available on arXiv:
<https://arxiv.org/abs/2511.04275>

The code implements the proposed **RetroAdj** method with efficient leave-one-out computation, as well as several competing online conformal prediction algorithms.
This repository provides scripts to reproduce the synthetic experiments and real-data applications reported in the paper.

## Repository structure

```text
.
├── R/                   # Core functions for Conformal Methods
├── aci/                 # Miscoverage Update Algorithms
├── simulations/         # Scripts for experiments and figure reproduction
├── call_dependencies.R  # Installs and loads all required R packages
└── .gitignore
```

## Citation
If you find this repository useful, in addition to the relevant methods, please cite:
```text
@article{jun2025online,
  title={Online Conformal Inference with Retrospective Adjustment for Faster Adaptation to Distribution Shift},
  author={Jun, Jungbin and Ohn, Ilsang},
  journal={arXiv preprint arXiv:2511.04275},
  year={2025}
}
```

## Contact

For questions about the implementation or replication of the paper, please contact:  
- Jungbin Jun - jungbini03@inha.edu

## Acknowledgements

The implementations of **DtACI** and **AgACI** algorithms included in this repository are adapted from the implemented code released by **Issac Gibbs**. [https://github.com/isgibbs/DtACI](https://github.com/isgibbs/DtACI)

We gratefully acknowledge the authors for making their work publicly available,  
which served as a foundation for reproducing and extending the ACI-family algorithms in this study.
