"""Dump a DSK HDF5 counts file to TSV and report its physical layout.

h5py cannot map HDF5's 128-bit integer to a numpy dtype, so the compound is read
through the low-level API as an opaque blob and decoded by hand. That decoding is
exactly what the Rust reader has to do, so this doubles as the spec for it.
"""
import sys, h5py, numpy as np
from h5py import h5d, h5t, h5s

path = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else None

f = h5py.File(path, "r")
print("kmer_size attr:", f["/dsk"].attrs["kmer_size"])
g = f["/dsk/solid"]
names = sorted(g.keys(), key=int)
print("partitions:", names)

did = g[names[0]].id
ftype = did.get_type()
print("compound size:", ftype.get_size(), "nmembers:", ftype.get_nmembers())
for i in range(ftype.get_nmembers()):
    mt = ftype.get_member_type(i)
    print(f"  [{i}] name={ftype.get_member_name(i).decode()} offset={ftype.get_member_offset(i)}"
          f" size={mt.get_size()} class={mt.get_class()}")
dcpl = did.get_create_plist()
print("layout:", dcpl.get_layout(), "nfilters:", dcpl.get_nfilters())

RECSIZE = ftype.get_size()
opaque = h5t.create(h5t.OPAQUE, RECSIZE)

total, rows = 0, []
for n in names:
    did = g[n].id
    npts = did.shape[0]
    buf = np.zeros(npts * RECSIZE, dtype=np.uint8)
    sp = h5s.create_simple((npts,))
    did.read(sp, did.get_space(), buf.view(np.dtype((np.void, RECSIZE))), opaque)
    raw = buf.tobytes()
    total += npts
    for i in range(npts):
        rec = raw[i*RECSIZE:(i+1)*RECSIZE]
        val = int.from_bytes(rec[0:16], "little")
        ab = int.from_bytes(rec[16:20], "little")
        rows.append((n, val, ab))
print("total solid kmers:", total)
print("first 4:", rows[:4])
if out:
    with open(out, "w") as fh:
        fh.write("partition\tkmer_u128\tabundance\n")
        for p, v, a in rows:
            fh.write(f"{p}\t{v}\t{a}\n")
    print("wrote", out, len(rows), "rows")
