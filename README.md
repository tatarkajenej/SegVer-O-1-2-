These are computations to accompany the paper "Secant varieties of Segre-Veronese varieties $\mathbb P^m\times\mathbb P^n$ embedded by $\mathcal O(1,2)$ are non-defective for $n\gg m^3$, $m\geq3$" ([arXiv:2503.21972](https://arxiv.org/abs/2503.21972)), written as part of the [TENORS network](https://github.com/tenors-network).

This work has been supported by European Union’s HORIZON–MSCA-2023-DN-JD programme under under 
the Horizon Europe (HORIZON) Marie Skłodowska-Curie Actions, grant agreement 101120296 (TENORS).


Please find the C++ source codes in this directory as well as (instances of) all relevant certificates in the subdirectory `certificates`.

The code is undoubtedly suboptimal in some aspects, so suggestions and comments are welcome.


# Dependencies

The code depends on:

* Eigen 3.4.0; see https://eigen.tuxfamily.org
* FFLAS-FFPACK 2.5.0; see https://linbox-team.github.io/fflas-ffpack/
* Givaro 4.2.1; see https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro/
* OpenMP; included in the GCC compiler
* a BLAS implementation, we used OpenBLAS 0.3.26; see http://www.openmathlib.org/OpenBLAS/

(The versions given are the ones we used, other recent versions might also work.)

# Compiling

Two executables should be compiled, one for the nice case and one for the ugly case. This should be achievable on an up-to-date Linux distribution with the following command:
> `g++ <name>.cpp -o"build/<name>" -fopenmp -lblas -lgivaro`

where `<name>` is one of 'nice' or 'ugly'.
This is automated in the attached `Makefile`, so it is also possible to just use
> `make build/<name>`

or:
> `make all`

If some of the library dependencies are located elsewhere than `/usr/include/` or `/usr/local/include/`, add appropriate paths to the library locations to the compilation command (e.g. by modifying your `Makefile`).


# Executing

Both executables take four integer arguments and a fifth optional string argument like so:
> `build/<name> w1 w2 m n [certpath]`

The arguments `w1`, `w2` determine which family of coordinate configurations will be used (see below for details), `m`, `n` specify the ambient $\mathbb P^m \times \mathbb P^n$ in which the coordinate configuration will be constructed and the optional `certpath` argument is a path specifying where the certificate of the computation should be stored. If omitted, `certpath` defaults to `cert/<name>-w1-w2-m-n.run` (hence they do not mix with or overwrite the presupplied certificates in `certificates`).

Note that the input is not currated in any sophisticated way, therefore, bad inputs (e.g., `w1` or `w2` out of range, negative `m` or `n` or an ugly (`m`,`n`) passed to the executable for the nice case or vice versa) may variously result in segmentation faults or errors on failed assert statements.

For ease of use, we provide both cases with a python script that runs the computation on all of the induction base cases used in the article. These can be executed by
> `python3 <name>.py`

and produce a small log `<name>.log` that contain information about which of the attempted computations succeeded (found non-defectivity) or failed; by default, up to 3 attempts are made for each computation if non-defectivity is not found immediately. Additionally, these scripts replace the `<name>-w1-w2` portion of the default certificate paths with a more legible designation in line with the paper, e.g. `B0` or `C1hat`. Please beware that especially the nice cases may be quite memory-intensive.

Both scripts may also be run through make by:
> `make run`


### Indexing of families of coordinate configurations

In the nice case, `w1` corresponds to the 'letter' used to designate the family in the paper (i.e. 0 means $A$, 1 means $B$, etc. up to 6 means $G$), whereas `w2` corresponds to the index. For example, `build/nice 3 0 7 111` executes the computation for $D_0(7,111)$.

In the ugly case, `w1` still corresponds to the letter, but starts with 1, i.e. 1 means $\hat A$ up to 4 means $\hat D$. `w2` no longer corresponds to the index. This is because the parameter functions for most of the families involve quasipolynomial functions, so we implement them internally as two separate families. As a result, the index of the family is $\left\lceil\frac{w_2+1}{2}\right\rceil$.

For both cases, the indexing can also be explicitly seen from the dictionary `configuration_names` in the respective python scripts.

