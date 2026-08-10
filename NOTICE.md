# Notices and attribution

This repository contains original work by John Kitchin (Carnegie Mellon
University) under the BSD 3-Clause License in [`LICENSE`](LICENSE), together
with work derived from the Tennessee Eastman Process Fortran code.

## Derived work

The Rust process model in `crates/tepsim-core` and the control layer in
`crates/tepsim-control` are ports of Fortran originally written by James J.
Downs and Ernest F. Vogel of the Tennessee Eastman Company, as later modified
and distributed by the Large Scale Systems Research Laboratory at the
University of Illinois under Professor Richard D. Braatz.

That original code is licensed under the **University of Illinois/NCSA Open
Source License**, reproduced verbatim in [`LICENSE-NCSA`](LICENSE-NCSA).

Two notes on that file. Upstream distributions of this code have sometimes been
described as BSD 3-Clause; the actual license text is NCSA, which is an
MIT-style grant carrying three BSD-style conditions. The file is also
reproduced exactly as distributed, including its typographical artifacts (stray
`%` characters and irregular whitespace in the disclaimer), because a license
should be carried verbatim rather than tidied.

The NCSA conditions require that the copyright notice, the list of conditions,
and the disclaimers be retained in redistributions in **both source and binary
form**. For this project that means `LICENSE-NCSA` and this file are included
in the published Python wheels, in the WebAssembly bundle, and in any binary
release, not only in the source repository.

Unmodified copies of the original Fortran and of the published `d00` through
`d21` datasets are vendored under `reference/`, with provenance and checksums
recorded in `reference/README.org`. Those files are never modified.

## Citing this work

Per the upstream license, users should cite the original code:

- J. J. Downs and E. F. Vogel, *A plant-wide industrial process control
  problem*, Presented at the AIChE 1990 Annual Meeting, Session on Industrial
  Challenge Problems in Process Control, Paper #24a, Chicago, Illinois,
  November 14, 1990.
- J. J. Downs and E. F. Vogel, *A plant-wide industrial process control
  problem*, Computers and Chemical Engineering, 17:245-255 (1993).
  <https://doi.org/10.1016/0098-1354(93)80018-I>

and the modified code:

- E. L. Russell, L. H. Chiang, and R. D. Braatz, *Data-driven Techniques for
  Fault Detection and Diagnosis in Chemical Processes*, Springer-Verlag,
  London, 2000. <https://doi.org/10.1007/978-1-4471-0409-4>
- L. H. Chiang, E. L. Russell, and R. D. Braatz, *Fault Detection and Diagnosis
  in Industrial Systems*, Springer-Verlag, London, 2001.
  <https://doi.org/10.1007/978-1-4471-0347-9>
- L. H. Chiang, E. L. Russell, and R. D. Braatz, *Fault diagnosis in chemical
  processes using Fisher discriminant analysis, discriminant partial least
  squares, and principal component analysis*, Chemometrics and Intelligent
  Laboratory Systems, 50:243-252, 2000.
  <https://doi.org/10.1016/S0169-7439(99)00061-1>
- E. L. Russell, L. H. Chiang, and R. D. Braatz, *Fault detection in industrial
  processes using canonical variate analysis and dynamic principal component
  analysis*, Chemometrics and Intelligent Laboratory Systems, 51:81-93, 2000.
  <https://doi.org/10.1016/S0169-7439(00)00058-7>

A citation for this Rust implementation will be added when it is published.
