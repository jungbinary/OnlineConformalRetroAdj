# Online Conformal Inference with Retrospective Adjustment

> **Online Conformal Inference with Retrospective Adjustment for Faster Adaptation to Distribution Shift**  
> Jungbin Jun and Ilsang Ohn  
> *Pattern Recognition*, 2026 [[Link](https://www.sciencedirect.com/science/article/abs/pii/S0031320326013713)] [[arXiv](https://arxiv.org/abs/2511.04275)]

This repository provides code for the proposed *RetroAdj* method, including an efficient leave-one-out implementation, along with the baseline methods used in the paper. It also includes scripts for reproducing the synthetic experiments and real-data analyses reported in the manuscript.

## Repository Structure

```text
.
├── R/                         # Core R implementations
├── python/                    # Python helpers for River-based baselines
├── aci/                       # ACI implementations
├── experiments/               # All runnable experiments
│   ├── main/
│   │   ├── simulation_study/  # Main synthetic experiments
│   │   └── real_data/         # Main real-data experiments
│   └── appendix/              # Appendix experiments
│       ├── window_sensitivity/
│       └── runtime_benchmark/
├── setup.R                    # Loads required packages and sources the main code
├── reproducibility_notes.txt 
├── requirements.txt
├── README.md
└── .gitignore
```

## Reproducibility

To reproduce the experiments in the paper, start from `experiments/`.  
The main results are organized under `experiments/main/`, and additional analyses are provided under `experiments/appendix/`.

## Citation

If you find this repository useful for your research, please cite:
```bibtex
@article{jun2026online,
  title   = {Online conformal inference with retrospective adjustment for faster adaptation to distribution shift},
  author  = {Jun, Jungbin and Ohn, Ilsang},
  journal = {Pattern Recognition},
  volume  = {180},
  pages   = {114406},
  year    = {2026},
  doi     = {10.1016/j.patcog.2026.114406}
}
```

## Contact

For questions regarding the code or reproducibility of the results, please contact:

- Jungbin Jun — jungbini03@inha.edu

## Acknowledgements

The implementations of **DtACI** and **AgACI** included in this repository are adapted from the publicly available code released by **Isaac Gibbs**:

<https://github.com/isgibbs/DtACI>
