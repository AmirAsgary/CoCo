# CoCo

CoCo is a software suite for different **Co**nsensus **Co**rrection applications using spaced k-mer count profiles of short reads or contigs. CoCo is open source GPL-licensed software implemented in C++.

### Compile from source
To compile CoCo `git`, a `C++/11` capable compiler (e.g. `gcc 4.7+`, `clang 3.5+`) and `cmake` (3.10 or higher) are required. Afterwards, the CoCo binary will be located in the `build/` directory.

      git clone https://github.com/soedinglab/CoCo.git
      cd CoCo
      git submodule update --init
      mkdir build && cd build
      cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
      make -j 4

### Rust port

[`Coco-rust/`](Coco-rust/) contains a Rust port of this tool. It reproduces the
C++ binary's output byte for byte on every command, runs 2.7× faster on one core
and 61× faster on 36, and uses about a tenth of the memory. Same command-line
interface, plus a `--threads` option.

      cd Coco-rust/coco && cargo build --release

See [Coco-rust/README.md](Coco-rust/README.md) for installation and usage,
[Coco-rust/ALGORITHM.md](Coco-rust/ALGORITHM.md) for how the algorithm works, and
[Coco-rust/MODIFICATIONS.md](Coco-rust/MODIFICATIONS.md) for how the port was
verified against this implementation.
