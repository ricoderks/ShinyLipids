#!/usr/bin/env python3
"""
One-off converter: MS-DIAL .lbm2 lipid library -> SQLite.

The .lbm2 file is a sequence of LZ4-compressed MessagePack chunks, each holding
a MessagePack array of MspFormatCompoundInformationBean records serialised as
positional arrays. See the lbm2r package README for the field map.

MS/MS peaks are stored as a BLOB of little-endian float32 (mz, intensity) pairs
rather than as rows in a peak table: the source values are already float32, so
this is lossless, and it keeps the whole spectrum a single-row fetch.

Usage:
    python3 scripts/lbm2_to_sqlite.py DB/<file>.lbm2 DB/lipids.sqlite
"""

import os
import sqlite3
import struct
import sys
import time

import lz4.block
import msgpack

# Positional indices within a record array.
I_PRECURSOR_MZ = 1
I_CHROM_XS = 2
I_ION_MODE = 3
I_PEAKS = 4
I_NAME = 5
I_FORMULA_BEAN = 6
I_SMILES = 8
I_INCHIKEY = 9
I_ADDUCT = 10
I_CCS = 11
I_CLASS = 14

ION_MODE = {0: "Positive", 1: "Negative"}

# chrom_xs holds four axes as [key, [value, type, unit]]: 1 = retention time,
# 2 = retention index, 3 = drift time, 4 = m/z. In this library only the
# retention time is populated (100% of records, 0.325-18.697 min); the other
# three are -1 throughout, so only RT is carried across.
CHROM_RT = 1


def retention_time(rec):
    """Retention time in minutes, or None where the library leaves it unset."""
    for item in rec[I_CHROM_XS]:
        if (isinstance(item, list) and len(item) == 2
                and item[0] == CHROM_RT and isinstance(item[1], list)):
            v = item[1][0]
            return v if v is not None and v > -1 else None
    return None

SCHEMA = """
PRAGMA journal_mode = OFF;
PRAGMA synchronous  = OFF;
PRAGMA temp_store   = MEMORY;
PRAGMA cache_size   = -2000000;

CREATE TABLE lipid (
    id           INTEGER PRIMARY KEY,
    name         TEXT    NOT NULL COLLATE NOCASE,
    lipid_class  TEXT    NOT NULL COLLATE NOCASE,
    precursor_mz REAL    NOT NULL,
    adduct       TEXT,
    ion_mode     TEXT,
    ccs          REAL,
    retention_time REAL,
    formula      TEXT,
    exact_mass   REAL,
    inchikey     TEXT,
    smiles       TEXT,
    n_peaks      INTEGER NOT NULL,
    peaks        BLOB    NOT NULL
);
"""

# Built after the bulk load: creating indexes up front slows inserts badly.
#
# name_fts is a trigram index over *distinct* names (2.6M of the 5.3M rows),
# which turns a substring search from a 200-400 ms full scan into ~1 ms.
# Note: the trigram tokeniser cannot match patterns shorter than 3 characters,
# so the app must fall back to LIKE for 1- and 2-character queries.
INDEXES = """
CREATE INDEX ix_lipid_mz        ON lipid (precursor_mz);
CREATE INDEX ix_lipid_class_mz  ON lipid (lipid_class, precursor_mz);
CREATE INDEX ix_lipid_name      ON lipid (name COLLATE NOCASE);
CREATE INDEX ix_lipid_inchikey  ON lipid (inchikey);
CREATE INDEX ix_lipid_rt        ON lipid (retention_time);

CREATE VIRTUAL TABLE name_fts USING fts5(name, tokenize='trigram');
INSERT INTO name_fts(name) SELECT DISTINCT name FROM lipid;

-- Precomputed so the class dropdown does not GROUP BY 5.3M rows at startup.
CREATE TABLE lipid_class_summary (
    lipid_class TEXT PRIMARY KEY,
    n           INTEGER NOT NULL,
    min_mz      REAL,
    max_mz      REAL,
    min_rt      REAL,
    max_rt      REAL
);
INSERT INTO lipid_class_summary
SELECT lipid_class, COUNT(*), MIN(precursor_mz), MAX(precursor_mz),
       MIN(retention_time), MAX(retention_time)
  FROM lipid GROUP BY lipid_class;

CREATE TABLE adduct_summary (
    adduct   TEXT PRIMARY KEY,
    ion_mode TEXT,
    n        INTEGER NOT NULL
);
INSERT INTO adduct_summary
SELECT adduct, MIN(ion_mode), COUNT(*) FROM lipid GROUP BY adduct;
"""


def iter_chunks(path):
    """Yield each decompressed MessagePack chunk payload from an .lbm2 file."""
    size = os.path.getsize(path)
    pos = 0
    with open(path, "rb") as fh:
        while pos < size:
            fh.seek(pos)
            header = fh.read(11)
            if len(header) < 11:
                break
            if header[0] != 0xC9 or header[6] != 0xD2:
                raise ValueError(f"unexpected chunk header at offset {pos}")
            comp_len = struct.unpack(">I", header[1:5])[0]
            unc_len = struct.unpack(">i", header[7:11])[0]
            fh.seek(pos + 11)
            payload = fh.read(comp_len - 5)
            yield lz4.block.decompress(payload, uncompressed_size=unc_len)
            pos += 6 + comp_len


def rows(records, start_id):
    """Turn raw positional records into tuples ready for executemany()."""
    rid = start_id
    for rec in records:
        peaks = rec[I_PEAKS]
        blob = b"".join(
            struct.pack("<ff", pk[0], pk[1]) for pk in peaks
        )
        bean = rec[I_FORMULA_BEAN]
        adduct = rec[I_ADDUCT]
        yield (
            rid,
            rec[I_NAME],
            rec[I_CLASS],
            rec[I_PRECURSOR_MZ],
            adduct[2] if adduct else None,
            ION_MODE.get(rec[I_ION_MODE]),
            rec[I_CCS],
            retention_time(rec),
            bean[0] if bean else None,
            bean[1] if bean else None,
            rec[I_INCHIKEY],
            rec[I_SMILES],
            len(peaks),
            blob,
        )
        rid += 1


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]

    if os.path.exists(dst):
        os.remove(dst)

    con = sqlite3.connect(dst)
    con.executescript(SCHEMA)

    insert = (
        "INSERT INTO lipid (id, name, lipid_class, precursor_mz, adduct, "
        "ion_mode, ccs, retention_time, formula, exact_mass, inchikey, smiles, "
        "n_peaks, peaks) "
        "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)"
    )

    t0 = time.time()
    total = 0
    for n, raw in enumerate(iter_chunks(src), start=1):
        records = msgpack.unpackb(raw, raw=False, strict_map_key=False)
        del raw
        con.executemany(insert, rows(records, total))
        con.commit()
        total += len(records)
        del records
        print(f"  chunk {n}: {total:,} records loaded ({time.time() - t0:.0f}s)",
              flush=True)

    print("building indexes ...", flush=True)
    con.executescript(INDEXES)
    con.executescript("PRAGMA optimize;")
    con.commit()
    con.close()

    print(f"done: {total:,} records in {time.time() - t0:.0f}s -> {dst} "
          f"({os.path.getsize(dst) / 1e9:.2f} GB)")


if __name__ == "__main__":
    main()
