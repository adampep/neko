# Neko
![CI](https://github.com/ExtremeFLOW/neko/workflows/CI/badge.svg)
## Cloning the project

```bash
git clone https://github.com/ExtremeFLOW/neko
```

## Building the project
To build the project you will need: A Fortran compiler supporting the Fortran-08 standard, a working MPI installation and BLAS/lapack.
We use automake to build the project. These instructions should work in general, but as the project is quickly developing, things might change.

```bash
cd neko
./regen.sh
./configure --prefix=/path/to/neko_install --with-pfunit=/path/to/pFUnit/installed/PFUNIT-VERSION
make install
```
## Running examples
After the project has been built

```bash
cd examples/tgv
/path/to/neko_install/bin/makeneko tgv.f90
mpirun -np 4 ./neko tgv.case
```

## Testing the Code
Assuming you configured with pFUnit you should be able to test the code with
```bash
make check
```

## Documentation
To generate the documentation, you need to have both doxygen and dot installed (they will be picked up by configure). Once installed, you should be able to generate the documentation with
```bash
make html
```