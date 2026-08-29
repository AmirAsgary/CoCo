//! Path and file helpers, ported from `src/filehandling.cpp`.
//!
//! `get_filename` is the one the commands use, to derive an output prefix from an
//! input path. The other two have no caller in the released C++; they are ported
//! so the translation is complete.
//!
//! `openFileOrDie`, `fileExists`, `directoryExists` and `_mkdir` have no separate
//! home here -- `std::fs` covers them directly at their call sites.

/// `getFilename`: drop the last extension, then take the basename.
///
/// Note the order: the extension goes first, so a directory containing a dot does
/// not confuse it, but a filename with several dots keeps all but the last.
pub fn get_filename(whole_file_path: &str) -> String {
    let stripped = match whole_file_path.rfind('.') {
        Some(i) => &whole_file_path[..i],
        None => whole_file_path,
    };
    match stripped.rfind('/') {
        Some(i) => stripped[i + 1..].to_string(),
        None => stripped.to_string(),
    }
}

/// `getFileExtension`: the last dot and everything after it, or empty.
pub fn get_file_extension(whole_file_path: &str) -> String {
    match whole_file_path.rfind('.') {
        Some(i) => whole_file_path[i..].to_string(),
        None => String::new(),
    }
}

/// `endsWith(suffix, str)`. Note the C++ takes the suffix first.
pub fn ends_with(suffix: &str, s: &str) -> bool {
    s.len() >= suffix.len() && s.ends_with(suffix)
}

/// `getFileList`: one path per line.
///
/// Blank lines are kept, as in the C++, which pushes whatever `getline` returns.
pub fn get_file_list(file_list_filename: &str) -> std::io::Result<Vec<String>> {
    let text = std::fs::read_to_string(file_list_filename)?;
    // `while (getline(fp, line))` yields no entry for a trailing newline, so a
    // file ending in "\n" does not produce a final empty string.
    let mut out: Vec<String> = text.split('\n').map(|s| s.to_string()).collect();
    if out.last().map(|s| s.is_empty()).unwrap_or(false) {
        out.pop();
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn get_filename_matches_the_cpp() {
        assert_eq!(get_filename("/a/b/reads.fq"), "reads");
        assert_eq!(get_filename("reads.err0.1pct.fq"), "reads.err0.1pct");
        assert_eq!(get_filename("noext"), "noext");
        assert_eq!(get_filename("/a/b/noext"), "noext");
        // The extension is stripped *before* the basename, so a dot in a directory
        // name truncates the path. `path_helpers_match_the_cpp` in tests/oracle.rs
        // pins the corner cases against the C++ rather than against a guess.
        assert_eq!(get_filename("/a.b/c"), "a");
    }

    #[test]
    fn get_file_extension_matches_the_cpp() {
        assert_eq!(get_file_extension("reads.fq"), ".fq");
        assert_eq!(get_file_extension("a.b.c"), ".c");
        assert_eq!(get_file_extension("noext"), "");
    }

    #[test]
    fn ends_with_matches_the_cpp() {
        assert!(ends_with(".gz", "reads.fq.gz"));
        assert!(!ends_with(".gz", "reads.fq"));
        // Shorter than the suffix is false, not a panic.
        assert!(!ends_with(".fastq", "a"));
        assert!(ends_with("", "anything"));
    }

    #[test]
    fn get_file_list_reads_one_path_per_line() {
        let p = std::env::temp_dir().join(format!("coco_filelist_{}", std::process::id()));
        std::fs::write(&p, "a.fq\nb.fq\n").unwrap();
        assert_eq!(get_file_list(p.to_str().unwrap()).unwrap(), vec!["a.fq", "b.fq"]);
        std::fs::write(&p, "a.fq\nb.fq").unwrap();
        assert_eq!(get_file_list(p.to_str().unwrap()).unwrap(), vec!["a.fq", "b.fq"]);
        let _ = std::fs::remove_file(&p);
    }
}
