# Example notebooks

Four Jupyter notebooks that use `tepsim`, the Python package built from this
repository, to run the Tennessee Eastman Process and to reproduce results from
its literature. They are committed with their outputs, so they can be read
without running anything.

| Notebook | What it covers |
|---|---|
| `01-getting-started.ipynb` | Running the plant, the 41 measurements and 12 manipulated variables, the twenty faults and what each one really does, injecting a fault, reading the ground-truth labels, reproducibility from a seed, and what a trip looks like. |
| `02-fault-detection-pca.ipynb` | The standard PCA monitoring scheme: T-squared and SPE, their control limits, false alarm rate, detection rate and detection delay. Checked digit for digit against the Rust reference implementation in `book/src/tutorials/a-detector.md`. |
| `03-hard-faults.ipynb` | The best known empirical result about this benchmark, that `IDV(3)`, `IDV(9)` and `IDV(15)` are effectively undetectable by these methods, measured three ways: on the published `d00`-`d21` files, on simulated runs at the seeds those files were generated with, and over a ten-seed ensemble. |
| `04-custom-scenarios.ipynb` | What this port can express that the original cannot: scheduled faults that arrive and clear, composed faults, continuous fault magnitudes, a choice of integrator, and the scenario text and digest that make a run reproducible from its description alone. |

`pcamon.py` beside them is the PCA monitoring implementation the second and
third notebooks share. It is NumPy and the Python standard library only: the
eigendecomposition is `numpy.linalg.eigh`, the normal quantile is
`statistics.NormalDist`, and the F quantile the T-squared limit needs is an
incomplete beta with a bisection on top.

## Running them

The notebooks need `tepsim`, `numpy`, `matplotlib` and `jupyter`. Nothing else:
no SciPy, no scikit-learn, no pandas.

If you have installed `tepsim` from a wheel, add the notebook tools and go:

```bash
pip install tepsim jupyter matplotlib
jupyter lab notebooks/
```

To run them against the working tree rather than a released wheel, build the
wheel and its virtualenv from the repository root, then add the notebook tools
to that virtualenv:

```bash
cargo xtask python                                  # builds .xtask-python/venv
.xtask-python/venv/bin/pip install jupyter matplotlib
.xtask-python/venv/bin/jupyter lab notebooks/
```

Note that `cargo xtask python` deletes and rebuilds `.xtask-python/venv`
every time it runs, so `jupyter` and `matplotlib` have to be reinstalled after
each rebuild.

To re-execute them non-interactively, which is how the committed outputs were
produced:

```bash
cd notebooks
../.xtask-python/venv/bin/jupyter nbconvert \
    --to notebook --execute --inplace *.ipynb
```

All four together take a few minutes, most of it in notebook 3, which runs a
ten-seed ensemble of 210 simulations of 48 hours each.

## Data

Notebook 3 reads the published `d00` through `d21` datasets from
`reference/data/`, which are vendored in this repository and are never edited.
The other three notebooks generate everything they use. Notebook 3 locates the
repository root by walking up from the working directory, so it works whether
it is run from `notebooks/` or from the root.

## References

The notebooks cite only these, and each one lists the ones it uses at the
bottom. They are the entries of `tep-rust.bib` in the repository root.

- J. J. Downs and E. F. Vogel, "A plant-wide industrial process control
  problem", *Computers & Chemical Engineering* **17**(3), 245-255 (1993).
  doi:10.1016/0098-1354(93)80018-I
- N. L. Ricker, "Decentralized control of the Tennessee Eastman challenge
  process", *Journal of Process Control* **6**(4), 205-221 (1996).
  doi:10.1016/0959-1524(96)00031-5
- E. L. Russell, L. H. Chiang and R. D. Braatz, "Fault detection in industrial
  processes using canonical variate analysis and dynamic principal component
  analysis", *Chemometrics and Intelligent Laboratory Systems* **51**(1), 81-93
  (2000). doi:10.1016/S0169-7439(00)00058-7
- L. H. Chiang, E. L. Russell and R. D. Braatz, *Fault Detection and Diagnosis
  in Industrial Systems*, Springer London (2001).
  doi:10.1007/978-1-4471-0347-9
- A. Bathelt, N. L. Ricker and M. Jelali, "Revision of the Tennessee Eastman
  process model", *IFAC-PapersOnLine* **48**(8), 309-314 (2015).
  doi:10.1016/j.ifacol.2015.08.199
- C. A. Rieth, B. D. Amsel, R. Tran and M. B. Cook, "Additional Tennessee
  Eastman process simulation data for anomaly detection evaluation", Harvard
  Dataverse (2017). doi:10.7910/DVN/6C3JR1
