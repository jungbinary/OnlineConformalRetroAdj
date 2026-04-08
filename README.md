# Online Conformal Inference with Retrospective Adjustment

> **Online Conformal Inference with Retrospective Adjustment for Faster Adaptation to Distribution Shift**  
> Jungbin Jun and Ilsang Ohn, 2025

Preprint available on arXiv: <https://arxiv.org/abs/2511.04275>

This repository provides code for the proposed RetroAdj method, including an efficient leave-one-out implementation, along with the baseline methods used in the paper. It also includes scripts for reproducing the synthetic experiments and real-data analyses reported in the manuscript.

## Repository Structure

```text
.
├── R/                         # Core R implementations
├── python/                    # Python helpers for River-based baselines
├── aci/                       # Adaptive conformal update rules
├── experiments/               # All runnable experiments
│   ├── main/
│   │   ├── simulation_study/  # Main synthetic experiments
│   │   └── real_data/         # Main real-data experiments
│   └── appendix/              # Appendix-only experiments
│       ├── window_sensitivity/
│       └── runtime_benchmark/
├── setup.R                    # Loads required packages and sources the main code
├── reproducibility_notes.txt  # Experiment map in manuscript order
├── requirements.txt
├── README.md
└── .gitignore
```

## Reproducibility

To reproduce the experiments in the paper, start from `experiments/`.  
The main results are organized under `experiments/main/`, and additional appendix analyses are provided under `experiments/appendix/`.

## Citation

If you use this repository in your research, please cite:

```bibtex
@article{jun2025online,
  title={Online Conformal Inference with Retrospective Adjustment for Faster Adaptation to Distribution Shift},
  author={Jun, Jungbin and Ohn, Ilsang},
  journal={arXiv preprint arXiv:2511.04275},
  year={2025}
}
```

## Contact

For questions regarding the code or reproducibility of the results, please contact:

- Jungbin Jun — jungbini03@inha.edu

## Acknowledgements

The implementations of **DtACI** and **AgACI** included in this repository are adapted from the publicly available code released by **Isaac Gibbs**:

<https://github.com/isgibbs/DtACI>
