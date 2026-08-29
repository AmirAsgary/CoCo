//! CoCo in Rust: consensus correction from spaced k-mer count profiles.
//!
//! A port of <https://github.com/soedinglab/CoCo> (Annika Jochheim,
//! MPI for Multidisciplinary Sciences). The module layout mirrors the C++ one
//! closely enough that a reader can hold both open side by side:
//!
//! | Rust | C++ |
//! | --- | --- |
//! | [`types`] | `src/types.{h,cpp}` |
//! | [`kmer`] | `src/kmer.{h,cpp}` |
//! | [`translator`] | `src/KmerTranslator.{h,cpp}` |
//! | [`lookuptable`] | `src/Lookuptable.{h,cpp}`, `src/HashTable.h` |
//! | [`dsk`] | GATB-core's `Storage`/`Partition<Count>` |
//! | [`seq`] | `lib/kseq/kseq.h`, `src/KSeqWrapper.cpp`, `src/SequenceInfo.h` |
//! | [`countprofile`] | `src/CountProfile.{h,cpp}` |
//! | [`runner`] | `src/runner.{h,cpp}` |
//! | [`filehandling`] | `src/filehandling.{h,cpp}` |
//! | [`options`] | `src/Options.{h,cpp}` |
//! | [`preprocessing`] | `src/preprocessing.{h,cpp}` |
//! | [`commands`] | `src/{correction,filter,abundanceEstimator,profile,counts2flat}.cpp` |
//!
//! Deliberate departures are documented where they occur; the significant ones are
//! the parallel read loop in [`runner`], the open addressing count table in
//! [`lookuptable`], the `PEXT` k-mer gather in [`translator`], and bounds checks
//! where the C++ relies on undefined behaviour.

pub mod commands;
pub mod countprofile;
pub mod dsk;
pub mod filehandling;
pub mod info;
pub mod kmer;
pub mod lookuptable;
pub mod options;
pub mod preprocessing;
pub mod runner;
pub mod seq;
pub mod sliding;
pub mod translator;
pub mod types;

pub const TOOL_NAME: &str = "CoCo";
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
pub const MAIN_AUTHOR: &str = "Annika Jochheim (annika.jochheim@mpinat.mpg.de)";
pub const TOOL_INTRODUCTION: &str = "CoCo is an open-source software suite for \
different COnsensus COrrection applications using spaced k-mer count profiles of \
short reads or contigs";
