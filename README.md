# QCxMS2
Program package for the quantum mechanical calculation of EI mass spectra using automated reaction network exploration.

[![License](https://img.shields.io/github/license/grimme-lab/QCxMS2)](https://github.com/grimme-lab/QCxMS2/blob/main/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/grimme-lab/QCxMS2)](https://github.com/grimme-lab/QCxMS2/releases/latest)


This is the download repository for the QCxMS2 program. 

<div align="center">
<img src="./assets/logo/qcxms2.svg" alt="Mass spectra calculation" width="420">
</div>

**Installation**

### Binary 

Statically linked binaries (Intel Compiler 21.3.0) can be found at the [latest release page](https://github.com/grimme-lab/QCxMS2/releases/latest).


Untar the zipped archive:

```bash
tar -xvzf QCxMS2.vX.X.tar.xz

The following files are being extracted: `qcxms2`
```
Place the executables into your ``$HOME/bin/`` directory or path.


> [!WARNING]
> For any installation make sure that you have correctly installed and sourced the following external programs before attempting any calculations with QCxMS2:

**required external programs:**

**xtb** (version >= 6.7.1)
[`xtb`](https://github.com/grimme-lab/xtb)


**CREST** (version >= 3.0.2)
[`crest`](https://github.com/crest-lab/crest)


**molbar** (version >= 1.1.3)
[`molbar`](https://git.rwth-aachen.de/bannwarthlab/molbar)


**orca** (version >= 6.0.0)
[`orca`](https://orcaforum.kofo.mpg.de)


**geodesic_interpolate** (version 1.0.0, optional)
[`geodesic_interpolate`](https://github.com/virtualzx-nad/geodesic-interpolate)





## Compilers 

ifort(<=2021.10.0), icc(<=2021.10.0)

### Meson

Using [meson](https://mesonbuild.com/) as build system requires you to install a fairly new version like 0.062 or newer.
To use the default backend of meson you have to install [ninja](https://ninja-build.org/) version 1.7 or newer.

```bash
export FC=ifort CC=icc
meson setup build -Dfortran_link_args="-lifcoremt -static" 
ninja -C build 
```

This will build a static linked binary in the ``build`` folder. Copy the binary from ``build/qcxms2`` file into a directory in your path, e.g. ``~/bin/``.


### Usage

To calculate an electron ionization mass spectrum with the default settings, run qcxms2 from the command line for an input molecule in an xyz coordinate file

```bash
qcxms2 in.xyz 
```

To see the various optional settings, run

```bash
qcxms2 --help
```


**Documentation**

A more detailed documentation on topics like input settings can be found at [read-the-docs](https://xtb-docs.readthedocs.io/en/latest/qcxms2_doc/qcxms2.html). 

Coming soon ...

**Citations**
1. J.Gorges, S. Grimme, *ChemRxiv*, **2025**, "QCxMS2 - a program for the calculation of electron ionization mass spectra via automated reaction network discovery".
  DOI: [10.26434/chemrxiv-2025-277zm](https://dx.doi.org/10.26434/chemrxiv-2025-277zm)


**License**

QCxMS2 is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

QCxMS2 is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.
