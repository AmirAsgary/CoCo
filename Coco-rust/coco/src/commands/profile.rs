//! `coco profile`, ported from `src/profile.cpp`. A developer tool that writes the
//! raw count profile of every read.

use crate::countprofile::CountProfile;
use crate::info::Info;
use crate::options::Options;
use crate::preprocessing::TableKind;
use crate::runner::{process_reads, ReadTask};
use crate::seq::SequenceInfo;
use crate::translator::KmerTranslator;

struct ProfileTask;

impl ReadTask for ProfileTask {
    type Stats = ();

    fn process(
        &self,
        cp: &mut CountProfile,
        seq: &mut SequenceInfo,
        skip: bool,
        out: &mut Vec<u8>,
        _stats: &mut (),
    ) {
        if skip {
            return;
        }
        out.push(b'#');
        out.extend_from_slice(&seq.name);
        if !seq.comment.is_empty() {
            out.push(b' ');
            out.extend_from_slice(&seq.comment);
        }
        out.push(b'\n');
        let mut num = itoa_buf();
        for (idx, &c) in cp.counts().iter().enumerate() {
            out.extend_from_slice(write_u64(&mut num, idx as u64));
            out.push(b'\t');
            out.extend_from_slice(write_u64(&mut num, c as u64));
            out.push(b'\n');
        }
    }

    fn merge(_lhs: &mut (), _rhs: &()) {}
}

/// Small decimal formatter; `format!` per line dominates this command otherwise.
fn itoa_buf() -> [u8; 20] {
    [0u8; 20]
}

fn write_u64(buf: &mut [u8; 20], mut v: u64) -> &[u8] {
    if v == 0 {
        buf[19] = b'0';
        return &buf[19..];
    }
    let mut i = 20;
    while v > 0 {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
    }
    &buf[i..]
}

pub fn run(opt: &Options, info: &Info) -> Result<i32, String> {
    let translator = KmerTranslator::new(&opt.spaced_kmer_pattern).map_err(|e| format!("{e}\n"))?;
    let table = super::make_lookuptable(opt, &translator, info, TableKind::Compact)?;

    info.info("Step 2: Print profiles...\n");
    for (src, suffix) in [
        (&opt.reads, ".profile.reads.txt"),
        (&opt.forward_reads, ".profile.1.txt"),
        (&opt.reverse_reads, ".profile.2.txt"),
    ] {
        if src.is_empty() {
            continue;
        }
        let path = super::out_path(opt, src, suffix);
        let mut w = std::io::BufWriter::with_capacity(8 << 20, std::fs::File::create(&path).map_err(|e| e.to_string())?);
        info.info(&format!("...process reads from file {src}\n"));
        process_reads(src, table.as_ref(), &translator, &ProfileTask, opt.skip as i64,
                      opt.threads as usize, &mut w).map_err(|e| format!("{e}\n"))?;
        use std::io::Write;
        w.flush().map_err(|e| e.to_string())?;
        info.info("...completed\n");
    }
    Ok(0)
}
