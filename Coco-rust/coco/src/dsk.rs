//! Reading DSK k-mer counts out of a GATB HDF5 file.
//!
//! In the C++ this is GATB-core's job (`Storage`, `Partition<Count>`,
//! `LargeInt<2>`), which is the only reason CoCo depends on GATB at all. The file
//! itself is simple: a `/dsk` group carrying a `kmer_size` attribute, and a
//! `/dsk/solid` group holding one dataset per DSK partition. Each dataset is a
//! compound of a 128-bit little-endian integer `value` and a `u32` `abundance`.
//!
//! `hdf5-metno` maps HDF5 types onto Rust types and has nothing to map a 128-bit
//! integer to, so the records are read through the raw C API using the file's own
//! datatype as the memory type. HDF5 treats that as an identity conversion and
//! hands back the bytes exactly as stored, which is all this needs.

use hdf5_sys::h5::hsize_t;
use hdf5_sys::h5d::{H5Dget_space, H5Dget_type, H5Dread};
use hdf5_sys::h5i::H5I_INVALID_HID;
use hdf5_sys::h5p::H5P_DEFAULT;
use hdf5_sys::h5s::{H5Sget_simple_extent_npoints, H5S_ALL};
use hdf5_sys::h5t::{H5Tclose, H5Tget_member_name, H5Tget_member_offset, H5Tget_member_type,
                    H5Tget_nmembers, H5Tget_size};
use std::ffi::CStr;

/// One `(canonical k-mer, abundance)` record, `Kmer<>::Count` in GATB.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DskCount {
    pub value: u128,
    pub abundance: u32,
}

#[derive(Debug)]
pub struct DskError(pub String);
impl std::fmt::Display for DskError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for DskError {}

impl From<hdf5_metno::Error> for DskError {
    fn from(e: hdf5_metno::Error) -> Self {
        DskError(e.to_string())
    }
}

/// Serialises every call into libhdf5.
///
/// The vendored HDF5 is built without `--enable-threadsafe`, so concurrent calls
/// from different threads corrupt its internal state. In the binary this only ever
/// runs once on the main thread, but the library is also driven from tests that
/// run in parallel, and a library should not silently require single-threaded use.
static HDF5_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn hdf5_guard() -> std::sync::MutexGuard<'static, ()> {
    HDF5_LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// A DSK counts file, opened and ready to stream its solid k-mers.
pub struct DskCountsFile {
    file: hdf5_metno::File,
    pub kmer_size: u32,
    partitions: Vec<String>,
}

/// Where the two compound members live inside a record, discovered from the file
/// rather than assumed. GATB pads `Count` out to 32 bytes on x86-64; relying on
/// that number would break the moment a file came from a differently padded build.
#[derive(Debug, Clone, Copy)]
struct RecordLayout {
    size: usize,
    value_off: usize,
    value_size: usize,
    abundance_off: usize,
    abundance_size: usize,
}

/// Read a string-valued attribute regardless of how it was encoded.
fn read_string_attr(attr: &hdf5_metno::Attribute) -> Option<String> {
    use hdf5_metno::types::{FixedAscii, FixedUnicode, VarLenAscii, VarLenUnicode};
    if let Ok(v) = attr.read_raw::<VarLenAscii>() {
        return v.first().map(|s| s.as_str().to_string());
    }
    if let Ok(v) = attr.read_raw::<VarLenUnicode>() {
        return v.first().map(|s| s.as_str().to_string());
    }
    if let Ok(v) = attr.read_raw::<FixedAscii<64>>() {
        return v.first().map(|s| s.as_str().to_string());
    }
    if let Ok(v) = attr.read_raw::<FixedUnicode<64>>() {
        return v.first().map(|s| s.as_str().to_string());
    }
    None
}

impl Drop for DskCountsFile {
    fn drop(&mut self) {
        // Closing the file is an HDF5 call too.
        let _g = hdf5_guard();
    }
}

impl DskCountsFile {
    pub fn open(path: &str) -> Result<Self, DskError> {
        let _g = hdf5_guard();
        let file = hdf5_metno::File::open(path)
            .map_err(|e| DskError(format!("cannot open count file {path}: {e}")))?;
        let dsk = file
            .group("dsk")
            .map_err(|e| DskError(format!("{path} has no /dsk group: {e}")))?;
        let attr = dsk
            .attr("kmer_size")
            .map_err(|e| DskError(format!("{path} has no dsk/kmer_size attribute: {e}")))?;
        // GATB writes this as a variable-length ASCII string. Which spelling HDF5
        // reports depends on the writing library, so try each in turn.
        let text = read_string_attr(&attr)
            .ok_or_else(|| DskError(format!("cannot read dsk/kmer_size from {path}")))?;
        let kmer_size: u32 = text
            .trim()
            .trim_end_matches(char::from(0))
            .parse()
            .map_err(|_| DskError(format!("kmer_size attribute {text:?} is not a number")))?;

        let solid = file
            .group("dsk/solid")
            .map_err(|e| DskError(format!("{path} has no dsk/solid group: {e}")))?;
        let mut partitions = solid.member_names()?;
        // DSK names partitions "0", "1", ...; the C++ iterator walks them in that
        // numeric order, and the lookup table's construction order follows.
        partitions.sort_by_key(|n| n.parse::<u64>().unwrap_or(u64::MAX));

        Ok(DskCountsFile { file, kmer_size, partitions })
    }

    pub fn num_partitions(&self) -> usize {
        self.partitions.len()
    }

    /// Total number of solid k-mers across all partitions.
    pub fn num_items(&self) -> Result<usize, DskError> {
        let _g = hdf5_guard();
        let mut n = 0usize;
        for p in &self.partitions {
            let ds = self.file.dataset(&format!("dsk/solid/{p}"))?;
            n += ds.shape().first().copied().unwrap_or(0);
        }
        Ok(n)
    }

    /// Call `f` once per DSK partition, with that partition's records decoded.
    ///
    /// Partitions arrive in the C++ iterator's order. Handing over a whole
    /// partition at a time lets the caller translate and insert it in parallel
    /// while peak memory stays at one partition's worth of records.
    pub fn for_each_partition<F: FnMut(&[DskCount])>(&self, mut f: F) -> Result<usize, DskError> {
        let _g = hdf5_guard();
        let mut total = 0usize;
        let mut decoded: Vec<DskCount> = Vec::new();
        for p in &self.partitions {
            let ds = self.file.dataset(&format!("dsk/solid/{p}"))?;
            let (layout, bytes) = read_raw_records(ds.id())?;
            let n = if layout.size == 0 { 0 } else { bytes.len() / layout.size };
            decoded.clear();
            decoded.reserve(n);
            for i in 0..n {
                decoded.push(decode(&bytes[i * layout.size..(i + 1) * layout.size], &layout));
            }
            f(&decoded);
            total += n;
        }
        Ok(total)
    }

    /// Call `f` once per solid k-mer, in the same order the C++ iterator yields.
    ///
    /// Streams a partition at a time so peak memory stays at one partition's worth
    /// of records rather than the whole file's.
    pub fn for_each<F: FnMut(DskCount)>(&self, mut f: F) -> Result<usize, DskError> {
        let _g = hdf5_guard();
        let mut total = 0usize;
        for p in &self.partitions {
            let ds = self.file.dataset(&format!("dsk/solid/{p}"))?;
            let recs = read_raw_records(ds.id())?;
            let (layout, bytes) = recs;
            let n = if layout.size == 0 { 0 } else { bytes.len() / layout.size };
            for i in 0..n {
                let rec = &bytes[i * layout.size..(i + 1) * layout.size];
                f(decode(rec, &layout));
            }
            total += n;
        }
        Ok(total)
    }
}

fn decode(rec: &[u8], l: &RecordLayout) -> DskCount {
    let mut v = [0u8; 16];
    let vs = l.value_size.min(16);
    v[..vs].copy_from_slice(&rec[l.value_off..l.value_off + vs]);
    let mut a = [0u8; 4];
    let as_ = l.abundance_size.min(4);
    a[..as_].copy_from_slice(&rec[l.abundance_off..l.abundance_off + as_]);
    DskCount { value: u128::from_le_bytes(v), abundance: u32::from_le_bytes(a) }
}

/// Read a whole dataset as raw bytes in its on-disk record layout.
fn read_raw_records(dset_id: i64) -> Result<(RecordLayout, Vec<u8>), DskError> {
    unsafe {
        let ftype = H5Dget_type(dset_id);
        if ftype == H5I_INVALID_HID {
            return Err(DskError("H5Dget_type failed".into()));
        }
        let layout = describe(ftype)?;

        let space = H5Dget_space(dset_id);
        if space == H5I_INVALID_HID {
            H5Tclose(ftype);
            return Err(DskError("H5Dget_space failed".into()));
        }
        let npoints = H5Sget_simple_extent_npoints(space);
        if npoints < 0 {
            H5Tclose(ftype);
            return Err(DskError("H5Sget_simple_extent_npoints failed".into()));
        }
        let mut buf = vec![0u8; npoints as usize * layout.size];
        // Passing the file datatype as the memory datatype makes this an identity
        // conversion, so HDF5 copies the records verbatim.
        let status = H5Dread(
            dset_id,
            ftype,
            H5S_ALL,
            H5S_ALL,
            H5P_DEFAULT,
            buf.as_mut_ptr() as *mut std::ffi::c_void,
        );
        H5Tclose(ftype);
        if status < 0 {
            return Err(DskError("H5Dread failed on a dsk/solid partition".into()));
        }
        let _ = npoints as hsize_t;
        Ok((layout, buf))
    }
}

unsafe fn describe(ftype: i64) -> Result<RecordLayout, DskError> {
    let size = H5Tget_size(ftype);
    let n = H5Tget_nmembers(ftype);
    if n < 2 {
        return Err(DskError(format!(
            "dsk/solid records have {n} members, expected a (value, abundance) compound"
        )));
    }
    let mut value_off = usize::MAX;
    let mut value_size = 0usize;
    let mut abundance_off = usize::MAX;
    let mut abundance_size = 0usize;
    for i in 0..n {
        let raw = H5Tget_member_name(ftype, i as u32);
        if raw.is_null() {
            continue;
        }
        let name = CStr::from_ptr(raw).to_string_lossy().into_owned();
        let off = H5Tget_member_offset(ftype, i as u32);
        let mt = H5Tget_member_type(ftype, i as u32);
        let msz = H5Tget_size(mt);
        H5Tclose(mt);
        hdf5_sys::h5::H5free_memory(raw as *mut std::ffi::c_void);
        match name.as_str() {
            "value" => {
                value_off = off;
                value_size = msz;
            }
            "abundance" => {
                abundance_off = off;
                abundance_size = msz;
            }
            _ => {}
        }
    }
    if value_off == usize::MAX || abundance_off == usize::MAX {
        return Err(DskError(
            "dsk/solid compound is missing a 'value' or 'abundance' member".into(),
        ));
    }
    if value_size > 16 {
        return Err(DskError(format!(
            "dsk/solid 'value' is {value_size} bytes; only k-mers up to 128 bits are supported"
        )));
    }
    Ok(RecordLayout { size, value_off, value_size, abundance_off, abundance_size })
}
