//! Command-line parsing, ported from `src/Options.cpp`.
//!
//! Option names, defaults, validation order and the parameter banner are kept
//! identical to the C++ so that scripts driving CoCo do not have to change and the
//! two binaries' stdout can be diffed directly.

use crate::translator::DEFAULT_PATTERN;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CommandId {
    Correction,
    Filter,
    Abundance,
    Profile,
    Counts2Flat,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OptType {
    Str,
    Int,
    Dbl,
    Bool,
}

struct OptSpec {
    name: &'static str,
    display: &'static str,
    description: &'static str,
    ty: OptType,
    is_file: bool,
}

macro_rules! opts {
    ($($id:ident => ($name:expr, $display:expr, $descr:expr, $ty:expr, $file:expr)),* $(,)?) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
        pub enum OptId { $($id),* }
        const SPECS: &[(OptId, OptSpec)] = &[
            $((OptId::$id, OptSpec { name: $name, display: $display, description: $descr,
                                     ty: $ty, is_file: $file })),*
        ];
    };
}

opts! {
    Reads => ("--reads", "Unpaired Reads",
        "file with unpaired/single/merged reads (fasta/fastq format)", OptType::Str, true),
    ForwardReads => ("-1", "Forward Reads",
        "file with forward paired-end reads (fasta/fastq format)", OptType::Str, true),
    ReverseReads => ("-2", "Reverse Reads",
        "file with reverse paired-end reads (fasta/fastq format)", OptType::Str, true),
    CountFile => ("--counts", "Count File",
        "pre computed solid kmer count file in hdf5 format (dsk output)", OptType::Str, true),
    Outdir => ("--outdir", "Outdir", "output directory", OptType::Str, false),
    Outprefix => ("--outprefix", "Outprefix", "prefix to use for resultfile(s)", OptType::Str, false),
    SpacedKmerPattern => ("--spaced-pattern", "Spaced pattern",
        "user-specified spaced k-mer pattern\n(span must be <=64, 12 <= weight <=32 and symmetric)",
        OptType::Str, false),
    Skip => ("--skip", "skip", "skip sequences with less than this many k-mers", OptType::Int, false),
    Threshold => ("--threshold", "Threshold",
        "percent drop threshold relative to neighborhood counts", OptType::Dbl, false),
    // `filter` swaps in its own wording for --threshold; see FILTER_THRESHOLD_DESCR.
    Pseudocount => ("--pseudocount", "Pseudocount",
        "untrusted count added to the pseudocount parameter", OptType::Int, false),
    LowerBound => ("--lowerbound", "Lower bound",
        "lower bound for neighborhood counts to be considered", OptType::Int, false),
    MaxCorrNum => ("--max-corr-num", "Max number of correction per read",
        "maximal number of corrections performed per read, changes are discarded otherwise",
        OptType::Int, false),
    MaxTrimLen => ("--max-trim-len", "Max number of trimmed nucleotides",
        "maximal number of nucleotides trimmed from the beginning or end of a read if error could not be corrected",
        OptType::Int, false),
    UpdateLookup => ("--update-lookup", "Update lookup table",
        "update counts in lookuptable after a sequence is corrected\n(slow down and creates read order dependency but might help in low coverage regions)",
        OptType::Bool, false),
    DropLevel1 => ("--drop-level1", "Drop-level1", "local drop criterion (range 0-0.33)", OptType::Dbl, false),
    DropLevel2 => ("--drop-level2", "Drop-level2", "global drop criterion (range 0-0.33)", OptType::Dbl, false),
    SoftFilter => ("--soft", "Soft Filtering",
        "less strict filtering mode due to more strict masking strategy ", OptType::Bool, false),
    Aligned => ("--aligned", "Aligned",
        "optimize abundance estimation for reads that span the same region (amplicon sequence data)",
        OptType::Bool, false),
    Threads => ("--threads", "Number of threads", "number of threads", OptType::Int, false),
    Verbose => ("--verbose", "Verbosity level",
        "verbosity level, 0: quiet 1: Errors, 2: +Warnings, 3: +Info, 4: +Debug", OptType::Int, false),
    CountMode => ("--count-mode", "Count mode",
        "way to store counts for concurrent spaced kmers (expert option)\n 0: sum \n 1: maximize",
        OptType::Int, false),
}

fn spec(id: OptId) -> &'static OptSpec {
    &SPECS.iter().find(|(i, _)| *i == id).unwrap().1
}

/// `filter()` replaces the --threshold help text with this one.
pub const FILTER_THRESHOLD_DESCR: &str = "percent drop threshold between two successive windows";

pub const COUNT_MODE_SUM: i32 = 0;
pub const COUNT_MODE_MAX: i32 = 1;

#[derive(Debug, Clone)]
pub struct Options {
    pub reads: String,
    pub forward_reads: String,
    pub reverse_reads: String,
    pub count_file: String,
    pub outdir: String,
    pub outprefix: String,
    pub spaced_kmer_pattern: String,

    pub count_mode: i32,
    pub skip: i32,
    pub threshold: f64,
    pub pseudocount: i32,
    pub lower_bound: i32,
    pub max_trim_len: i32,
    pub max_corr_num: i32,
    pub update_lookup: bool,
    pub drop_level1: f64,
    pub drop_level2: f64,
    pub threads: i32,
    pub verbose: i32,
    pub aligned: bool,
    pub soft_filter: bool,

    set: std::collections::HashSet<OptId>,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            reads: String::new(),
            forward_reads: String::new(),
            reverse_reads: String::new(),
            count_file: String::new(),
            outprefix: String::new(),
            outdir: "coco_out/".into(),
            spaced_kmer_pattern: DEFAULT_PATTERN.into(),
            skip: 10,
            pseudocount: 1,
            threshold: 0.01,
            lower_bound: 5,
            max_trim_len: 0,
            max_corr_num: 10,
            update_lookup: false,
            drop_level1: 0.33,
            drop_level2: 0.33,
            aligned: false,
            soft_filter: false,
            count_mode: COUNT_MODE_SUM,
            // The C++ default is 1 because its thread support was never finished.
            // Here the work really is parallel, so default to the machine.
            threads: 0,
            verbose: 3,
            set: std::collections::HashSet::new(),
        }
    }
}

/// Option order per command, matching the `*Workflow` vectors in `Options.cpp`.
pub fn workflow(cmd: CommandId) -> &'static [OptId] {
    use OptId::*;
    match cmd {
        CommandId::Correction => &[
            ForwardReads, ReverseReads, Reads, CountFile, CountMode, Outdir, Outprefix,
            SpacedKmerPattern, Skip, Threshold, Pseudocount, LowerBound, MaxCorrNum, MaxTrimLen,
            UpdateLookup, Threads, Verbose,
        ],
        CommandId::Filter => &[
            ForwardReads, ReverseReads, Reads, CountFile, CountMode, Outdir, Outprefix,
            SpacedKmerPattern, Skip, Threshold, Threads, Verbose,
        ],
        CommandId::Abundance => &[
            ForwardReads, ReverseReads, Reads, CountFile, CountMode, Outdir, Outprefix,
            SpacedKmerPattern, Skip, Threads, Verbose,
        ],
        CommandId::Profile => &[
            ForwardReads, ReverseReads, Reads, CountFile, CountMode, Outdir, Outprefix,
            SpacedKmerPattern, Skip, Threads, Verbose,
        ],
        CommandId::Counts2Flat => &[CountFile, CountMode, Outdir, Outprefix, SpacedKmerPattern],
    }
}

#[derive(Debug)]
pub enum ParseOutcome {
    Ok(Box<Options>),
    /// `-h`/`--help` was given; the usage text has already been produced.
    HelpRequested(String),
    Error(String),
}

impl Options {
    pub fn is_set(&self, id: OptId) -> bool {
        self.set.contains(&id)
    }

    fn value_string(&self, id: OptId) -> String {
        use OptId::*;
        match id {
            Reads => self.reads.clone(),
            ForwardReads => self.forward_reads.clone(),
            ReverseReads => self.reverse_reads.clone(),
            CountFile => self.count_file.clone(),
            Outdir => self.outdir.clone(),
            Outprefix => self.outprefix.clone(),
            SpacedKmerPattern => self.spaced_kmer_pattern.clone(),
            Skip => self.skip.to_string(),
            Threshold => fmt_double(self.threshold),
            Pseudocount => self.pseudocount.to_string(),
            LowerBound => self.lower_bound.to_string(),
            MaxCorrNum => self.max_corr_num.to_string(),
            MaxTrimLen => self.max_trim_len.to_string(),
            UpdateLookup => (self.update_lookup as i32).to_string(),
            DropLevel1 => fmt_double(self.drop_level1),
            DropLevel2 => fmt_double(self.drop_level2),
            SoftFilter => (self.soft_filter as i32).to_string(),
            Aligned => (self.aligned as i32).to_string(),
            Threads => self.threads.to_string(),
            Verbose => self.verbose.to_string(),
            CountMode => self.count_mode.to_string(),
        }
    }

    fn set_value(&mut self, id: OptId, raw: &str) -> Result<(), String> {
        use OptId::*;
        match id {
            Reads => self.reads = raw.into(),
            ForwardReads => self.forward_reads = raw.into(),
            ReverseReads => self.reverse_reads = raw.into(),
            CountFile => self.count_file = raw.into(),
            Outdir => self.outdir = raw.into(),
            Outprefix => self.outprefix = raw.into(),
            SpacedKmerPattern => self.spaced_kmer_pattern = raw.into(),
            // `atoi` on a non-number yields 0 rather than failing; kept so a typo
            // behaves the same in both binaries.
            Skip => self.skip = atoi(raw),
            Pseudocount => self.pseudocount = atoi(raw),
            LowerBound => self.lower_bound = atoi(raw),
            MaxCorrNum => self.max_corr_num = atoi(raw),
            MaxTrimLen => self.max_trim_len = atoi(raw),
            Threads => self.threads = atoi(raw),
            Verbose => self.verbose = atoi(raw),
            CountMode => self.count_mode = atoi(raw),
            Threshold => self.threshold = parse_double(raw)?,
            DropLevel1 => self.drop_level1 = parse_double(raw)?,
            DropLevel2 => self.drop_level2 = parse_double(raw)?,
            UpdateLookup | SoftFilter | Aligned => {
                let v = parse_bool(raw, spec(id).name)?;
                match id {
                    UpdateLookup => self.update_lookup = v,
                    SoftFilter => self.soft_filter = v,
                    _ => self.aligned = v,
                }
            }
        }
        self.set.insert(id);
        Ok(())
    }

    /// `Options::printParameterSettings`.
    pub fn parameter_settings(&self, cmd: CommandId) -> String {
        let ids = workflow(cmd);
        let max_param_width = ids.iter().map(|&i| spec(i).display.len()).max().unwrap_or(0) + 6;
        let mut s = String::new();
        for &id in ids {
            let display = spec(id).display;
            let pad = if max_param_width < display.len() { 1 } else { max_param_width - display.len() };
            s.push_str("  ");
            s.push_str(display);
            s.push_str(&" ".repeat(pad));
            let v = self.value_string(id);
            if !v.is_empty() {
                s.push_str(" [");
                s.push_str(&v);
                s.push(']');
            }
            s.push('\n');
        }
        s
    }

    /// `printToolUsage(command, EXTENDED)`.
    pub fn tool_usage(&self, cmd: CommandId, author: &str, usage: &str, cmd_name: &str) -> String {
        let ids = workflow(cmd);
        let mut s = format!("© {author}\n\nUsage: coco {cmd_name} {usage}\n\n");
        let max_param_width = ids.iter().map(|&i| spec(i).name.len()).max().unwrap_or(0) + 6;
        for &id in ids {
            let sp = spec(id);
            let pad = if max_param_width < sp.name.len() { 1 } else { max_param_width - sp.name.len() };
            let mut line = format!("  {}{}", sp.name, " ".repeat(pad));
            let description_start = line.chars().count();
            let descr = if cmd == CommandId::Filter && id == OptId::Threshold {
                FILTER_THRESHOLD_DESCR
            } else {
                sp.description
            };
            let mut parts = descr.split('\n');
            line.push_str(parts.next().unwrap_or(""));
            let v = self.value_string(id);
            if !v.is_empty() {
                line.push_str(&format!(" [{}]", fmt_usage_value(id, self)));
            }
            for p in parts {
                line.push('\n');
                line.push_str(&" ".repeat(description_start));
                line.push_str(p);
            }
            s.push_str(&line);
            s.push('\n');
        }
        s.push('\n');
        s
    }
}

/// `std::to_string(double)`: fixed six decimal places.
fn fmt_double(v: f64) -> String {
    format!("{v:.6}")
}

/// The usage banner uses `sprintf("%.3lf")` instead.
fn fmt_usage_value(id: OptId, o: &Options) -> String {
    match spec(id).ty {
        OptType::Dbl => {
            let v = match id {
                OptId::Threshold => o.threshold,
                OptId::DropLevel1 => o.drop_level1,
                _ => o.drop_level2,
            };
            format!("{v:.3}")
        }
        _ => o.value_string(id),
    }
}

/// C `atoi`: leading whitespace, optional sign, digits, stop at the first
/// non-digit, 0 if nothing parsed.
fn atoi(s: &str) -> i32 {
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    let neg = i < b.len() && b[i] == b'-';
    if i < b.len() && (b[i] == b'-' || b[i] == b'+') {
        i += 1;
    }
    let mut v: i64 = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        v = v.saturating_mul(10).saturating_add((b[i] - b'0') as i64);
        if v > i32::MAX as i64 + 1 {
            v = i32::MAX as i64 + 1;
        }
        i += 1;
    }
    let v = if neg { -v } else { v };
    v.clamp(i32::MIN as i64, i32::MAX as i64) as i32
}

fn parse_double(s: &str) -> Result<f64, String> {
    s.trim()
        .parse::<f64>()
        .map_err(|_| format!("ERROR: Invalid numeric argument '{s}'"))
}

fn parse_bool(s: &str, name: &str) -> Result<bool, String> {
    match s {
        "true" | "TRUE" | "1" => Ok(true),
        "false" | "FALSE" | "0" => Ok(false),
        _ => Err(format!("ERROR: Invalid boolean argument for option {name}\n")),
    }
}

/// `Options::parseOptions`.
pub fn parse_options(
    argv: &[String],
    cmd: CommandId,
    author: &str,
    usage: &str,
    cmd_name: &str,
) -> ParseOutcome {
    let mut o = Options::default();
    // `filter()` in the C++ overwrites the threshold default *before* calling
    // parseOptions, so the filter command means something different by
    // --threshold: a percent drop between two successive windows rather than a
    // drop relative to neighbourhood counts. A user-supplied value still wins,
    // because this happens before parsing.
    if cmd == CommandId::Filter {
        o.threshold = 0.1;
    }
    let ids = workflow(cmd);

    if argv.is_empty() {
        return ParseOutcome::HelpRequested(o.tool_usage(cmd, author, usage, cmd_name));
    }
    if argv.iter().any(|a| a == "-h" || a == "--help") {
        return ParseOutcome::HelpRequested(o.tool_usage(cmd, author, usage, cmd_name));
    }

    let mut i = 0usize;
    while i < argv.len() {
        let optname = &argv[i];
        let found = ids.iter().copied().find(|&id| spec(id).name == optname);
        let Some(id) = found else {
            return ParseOutcome::Error(format!("ERROR: Unrecognized parameter {optname}\n"));
        };
        if o.set.contains(&id) {
            return ParseOutcome::Error(format!("ERROR: Duplicate option {}\n\n", spec(id).name));
        }
        let sp = spec(id);
        let next_is_flag = i + 1 >= argv.len() || argv[i + 1].starts_with('-');
        if sp.ty != OptType::Bool && next_is_flag {
            return ParseOutcome::Error(format!(
                "ERROR: No argument for option {} is given \n",
                sp.name
            ));
        }
        if sp.ty == OptType::Bool {
            // A bare boolean flag means true; an explicit value is consumed.
            if next_is_flag {
                if let Err(e) = o.set_value(id, "true") {
                    return ParseOutcome::Error(e);
                }
            } else {
                if let Err(e) = o.set_value(id, &argv[i + 1]) {
                    return ParseOutcome::Error(e);
                }
                i += 1;
            }
        } else {
            if let Err(e) = o.set_value(id, &argv[i + 1]) {
                return ParseOutcome::Error(e);
            }
            i += 1;
        }
        i += 1;
    }

    // Required-argument checks, in the C++'s order.
    if cmd != CommandId::Counts2Flat
        && !(o.is_set(OptId::Reads)
            || (o.is_set(OptId::ForwardReads) && o.is_set(OptId::ReverseReads)))
    {
        return ParseOutcome::Error(
            "ERROR: Either --reads or -1 and -2 must be set\n".to_string(),
        );
    }
    if cmd == CommandId::Counts2Flat && !o.is_set(OptId::CountFile) {
        return ParseOutcome::Error("ERROR: --counts must be set\n".to_string());
    }

    for &id in ids {
        let sp = spec(id);
        if o.is_set(id) && sp.is_file {
            let v = o.value_string(id);
            if !std::path::Path::new(&v).exists() {
                return ParseOutcome::Error(format!(
                    "ERROR: In Option {} file '{}' does not exist\n",
                    sp.name, v
                ));
            }
        }
    }

    if !o.outdir.is_empty() {
        if !o.outdir.ends_with('/') {
            o.outdir.push('/');
        }
        if !std::path::Path::new(&o.outdir).is_dir() {
            if let Err(e) = std::fs::create_dir_all(&o.outdir) {
                return ParseOutcome::Error(format!(
                    "ERROR: Failed to create output directory {} ({e})\n",
                    o.outdir
                ));
            }
        }
    }

    if o.is_set(OptId::Outprefix) && o.outprefix.contains('/') {
        return ParseOutcome::Error(
            "ERROR: Value for option --outprefix is a path, please use --outdir to set a directory\n"
                .to_string(),
        );
    }

    if o.threads <= 0 {
        o.threads = std::thread::available_parallelism().map(|n| n.get() as i32).unwrap_or(1);
    }
    if o.update_lookup && o.threads > 1 {
        // Not an error: the option is inherently sequential, so honour it and say so.
        o.threads = 1;
    }
    if o.verbose >= 4 {
        // The CDEBUG correction trace is emitted per read; keeping it in input
        // order means keeping the reads in input order.
        o.threads = 1;
    }

    ParseOutcome::Ok(Box::new(o))
}

pub use crate::filehandling::get_filename;

#[cfg(test)]
mod tests {
    use super::*;

    fn args(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn atoi_matches_c_semantics() {
        assert_eq!(atoi("42"), 42);
        assert_eq!(atoi("-7"), -7);
        assert_eq!(atoi("12abc"), 12);
        assert_eq!(atoi("abc"), 0);
        assert_eq!(atoi(""), 0);
        assert_eq!(atoi("  8 "), 8);
    }

    #[test]
    fn get_filename_strips_extension_and_directory() {
        assert_eq!(get_filename("/a/b/reads.fq"), "reads");
        assert_eq!(get_filename("reads.err0.1pct.fq"), "reads.err0.1pct");
        assert_eq!(get_filename("noext"), "noext");
    }

    #[test]
    fn defaults_match_the_cpp() {
        let o = Options::default();
        assert_eq!(o.spaced_kmer_pattern, DEFAULT_PATTERN);
        assert_eq!(o.skip, 10);
        assert_eq!(o.pseudocount, 1);
        assert_eq!(o.threshold, 0.01);
        assert_eq!(o.lower_bound, 5);
        assert_eq!(o.max_trim_len, 0);
        assert_eq!(o.max_corr_num, 10);
        assert!(!o.update_lookup);
        assert_eq!(o.drop_level1, 0.33);
        assert_eq!(o.drop_level2, 0.33);
        assert_eq!(o.count_mode, COUNT_MODE_SUM);
        assert_eq!(o.outdir, "coco_out/");
        assert_eq!(o.verbose, 3);
    }

    #[test]
    fn rejects_duplicate_and_unknown_options() {
        let a = args(&["--reads", "/dev/null", "--reads", "/dev/null"]);
        assert!(matches!(
            parse_options(&a, CommandId::Correction, "", "", "correction"),
            ParseOutcome::Error(e) if e.contains("Duplicate option")
        ));
        let a = args(&["--nope", "1"]);
        assert!(matches!(
            parse_options(&a, CommandId::Correction, "", "", "correction"),
            ParseOutcome::Error(e) if e.contains("Unrecognized parameter")
        ));
    }

    #[test]
    fn requires_reads_or_a_pair() {
        let a = args(&["--threshold", "0.02"]);
        assert!(matches!(
            parse_options(&a, CommandId::Correction, "", "", "correction"),
            ParseOutcome::Error(e) if e.contains("must be set")
        ));
    }

    #[test]
    fn trailing_boolean_flag_does_not_crash() {
        // The C++ guards `argv[argIdx+1]` against running off the end only for
        // non-boolean options, so a bare boolean flag in final position reads
        // argv[argc] -- NULL -- and segfaults. Every bool option is affected:
        // --update-lookup, --soft, --aligned.
        for a in [
            vec!["--reads", "/dev/null", "--update-lookup"],
            vec!["--reads", "/dev/null", "--max-corr-num", "3", "--update-lookup"],
        ] {
            let a: Vec<String> = a.iter().map(|s| s.to_string()).collect();
            match parse_options(&a, CommandId::Correction, "", "", "correction") {
                ParseOutcome::Ok(o) => assert!(o.update_lookup),
                other => panic!("unexpected {other:?}"),
            }
        }
    }

    #[test]
    fn bare_boolean_flag_is_true() {
        let a = args(&["--reads", "/dev/null", "--update-lookup", "--threshold", "0.5"]);
        match parse_options(&a, CommandId::Correction, "", "", "correction") {
            ParseOutcome::Ok(o) => {
                assert!(o.update_lookup);
                assert_eq!(o.threshold, 0.5);
                // update-lookup forces the sequential path.
                assert_eq!(o.threads, 1);
            }
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn parameter_banner_layout_matches_the_cpp() {
        let o = Options::default();
        let s = o.parameter_settings(CommandId::Counts2Flat);
        // Longest display among the counts2flat options is "Spaced pattern" (14),
        // so the column is 14 + 6 = 20 wide after the two-space indent.
        let first = s.lines().next().unwrap();
        assert!(first.starts_with("  Count File"));
        assert_eq!(&first[2..2 + 20], "Count File          ");
        // Empty values print no bracket at all.
        assert_eq!(first.trim_end(), "  Count File");
        assert!(s.contains(&format!(" [{DEFAULT_PATTERN}]")));
    }

    #[test]
    fn filter_gets_its_own_threshold_default() {
        // The C++ filter() sets opt.threshold = 0.1 before parsing, so the filter
        // command's default differs from every other command's.
        let a = args(&["--reads", "/dev/null"]);
        match parse_options(&a, CommandId::Filter, "", "", "filter") {
            ParseOutcome::Ok(o) => assert_eq!(o.threshold, 0.1),
            other => panic!("unexpected {other:?}"),
        }
        match parse_options(&a, CommandId::Correction, "", "", "correction") {
            ParseOutcome::Ok(o) => assert_eq!(o.threshold, 0.01),
            other => panic!("unexpected {other:?}"),
        }
        // An explicit value still wins.
        let a = args(&["--reads", "/dev/null", "--threshold", "0.42"]);
        match parse_options(&a, CommandId::Filter, "", "", "filter") {
            ParseOutcome::Ok(o) => assert_eq!(o.threshold, 0.42),
            other => panic!("unexpected {other:?}"),
        }
    }

    #[test]
    fn doubles_print_with_six_decimals() {
        let o = Options::default();
        let s = o.parameter_settings(CommandId::Correction);
        assert!(s.contains("[0.010000]"), "got:\n{s}");
    }
}
