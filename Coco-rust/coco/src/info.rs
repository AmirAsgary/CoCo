//! Verbosity-gated output, ported from `src/Info.h`.
//!
//! The levels and the stdout/stderr split are the C++'s: `INFO` goes to stdout,
//! everything else to stderr, and nothing is emitted above the configured level.

#[derive(Debug, Clone, Copy)]
pub struct Info {
    level: i32,
}

impl Info {
    pub const ERROR: i32 = 1;
    pub const WARNING: i32 = 2;
    pub const INFO: i32 = 3;
    pub const DEBUG: i32 = 4;

    pub fn new(verbose_level: i32) -> Self {
        Info { level: verbose_level }
    }

    pub fn info(&self, s: &str) {
        if Self::INFO <= self.level {
            use std::io::Write;
            let mut o = std::io::stdout();
            let _ = o.write_all(s.as_bytes());
        }
    }

    pub fn warning(&self, s: &str) {
        if Self::WARNING <= self.level {
            eprint!("{s}");
        }
    }

    pub fn error(&self, s: &str) {
        if Self::ERROR <= self.level {
            eprint!("{s}");
        }
    }

    pub fn debug(&self, s: &str) {
        if Self::DEBUG <= self.level {
            eprint!("{s}");
        }
    }

    pub fn level(&self) -> i32 {
        self.level
    }
}
