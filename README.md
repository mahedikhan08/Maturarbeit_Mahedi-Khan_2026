# Optimierung einer Protein-Design-Pipeline zur Erzeugung thermostabiler Proteine

> **High-Performance Computing & Deep Learning Workflows for Thermostable Protein Design**  
> *Graduation Project / Maturaarbeit (2026) – Developed in collaboration with the Schwede Group (Biozentrum, University of Basel)*

---

## 📌 Project Overview

Proteins are fundamental biological tools in medicine and industry, but their natural stability is usually limited to physiological conditions . This project focuses on the systematic calibration of a thermodynamic stability weight ($w_{\text{stab}}$) within the Monte-Carlo-based protein design pipeline **TEA-Leaves** . 

By balancing three fundamental forces—**biological plausibility** (ESM2 language model), **structural wildtype similarity** (TEA-Alphabet $T\%$), and **physical stability** (abs_dG-predictor)—we aim to generate thermostable protein sequences without compromising structural integrity or biological function .

### Key Scientific Findings
- **Scale:** Evaluated over **82 representative CATH protein families** with 100 independent designs per family across 4 weight groups ($w_{\text{stab}} = 0, 0.001, 0.01, 0.1$), totaling **32,800 simulated trajectories** .
- **Sweet Spot Identified:** A weight of **$w_{\text{stab}} = 0.01$** provides the optimal balance between thermodynamic stability gain ($\Delta G$) and high structural success rate ($\text{pTM} > 0.7$, $\text{TM-score} > 0.7$) .
- **Trade-off Boundary:** Setting $w_{\text{stab}} = 0.1$ over-prioritizes pure physics, leading to a $>50\%$ drop in structurally viable designs .

---

## 🗂 Repository Structure

The repository is organized into distinct modules for reproducibility and clear data provenance:

* cath_database/
  * cath_sequences_selected.fasta (82 target wildtype CATH sequences)
* Sequences/
  * all_designs/ (All 32,800 simulated trajectory FASTA files)
  * successful_designs/ (Designs filtered by global topology: pTM > 0.7)
  * lab/ (Top candidates filtered for lab validation: pTM > 0.7 & TM > 0.7)
* results/
  * Plots/ (Categorized analytical plots: scatter, violin, tendency)
  * dG_averages_summary.txt (Aggregated raw stability scores)
  * analysis_results.txt (Wilcoxon signed-rank tests & compromise score evaluation)
* utils/
  * generate_scripts/ (Python scripts for automated SLURM job generation)
  * plot_scripts/ (Data visualization and statistical plotting scripts)
  * submission_files/ (SLURM array submission bash scripts)

---

## 📊 Key Results Summary

| Weight Group ($w_{\text{stab}}$) | Total Designs | Successful Designs ($\text{TM} \ge 0.7$) | Success Rate (%) | Mean abs_dG (kcal/mol) | Compromise Score | Status |
| :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| **`dG_0`** (Control) | 8,200 | 363 | 4.43% | 2.743 | -0.038 | Baseline |
| **`dG_0.001`** | 8,200 | 354 | 4.32% | 2.665 | -0.198 | Ineffective |
| **`dG_0.01`** | 8,200 | **350** | **4.27%** | **3.366** | **+0.268** | **Optimal Sweet Spot** |
| **`dG_0.1`** | 8,200 | 153 | 1.87% | 6.029 | -0.032 | Severe quality drop |

*Statistical significance confirmed via paired Wilcoxon signed-rank tests (p < 0.001$ for w_stab = 0.01 vs dG_0).* 

---

## 💻 Infrastructure & Execution

All simulations were executed on the **sciCORE High-Performance Computing (HPC) cluster** at the University of Basel using GPU-accelerated nodes managed via **SLURM** .

> ⚠️ **Note on Pipeline Dependencies:**  
> The core `TEA-Leaves` pipeline package relies on internal software developed by the Schwede Group (Biozentrum, University of Basel) which is currently not publicly available . Therefore, the scripts in `utils/` are provided for **reproducibility and documentation purposes** of this research project . 

### Environment Setup
- **Python 3.10+**
- **PyTorch** & **Transformers** (Hugging Face)
- **ESM2** (Evolutionary Scale Modeling) & **ESMFold** 
- **Absolute Stability Predictor** (`ESM3AG`) 

To re-run job submission script generation:
python utils/generate_scripts/cath_generate_abs_stability.py

To automatically re-submit incomplete or failed array jobs:
python utils/generate_scripts/make_missing_jobs.py

---

## 📜 License & Citation

This repository contains code and data developed for a high school graduation project (Maturaarbeit) at Gymnasium Liestal in collaboration with the Schwede Group .

* **Author:** Mahedi Khan 
* **Supervision:** Dr. Gabriel Studer, Dr. Janani Durairaj (Biozentrum, University of Basel) 
