"""Minimal reader for netCDF-3 classic / 64-bit-offset files (no dependencies beyond numpy).

Only supports fixed-size (non-record) variables, which is all ERDDAP griddap .nc output needs.
"""
import struct
import numpy as np

_TYPES = {1: ('i1', 1), 2: ('S1', 1), 3: ('>i2', 2), 4: ('>i4', 4), 5: ('>f4', 4), 6: ('>f8', 8)}


class _Reader:
    def __init__(self, buf):
        self.buf = buf
        self.pos = 0

    def int(self):
        v = struct.unpack_from('>i', self.buf, self.pos)[0]
        self.pos += 4
        return v

    def int64(self):
        v = struct.unpack_from('>q', self.buf, self.pos)[0]
        self.pos += 8
        return v

    def name(self):
        n = self.int()
        s = self.buf[self.pos:self.pos + n].decode('utf-8')
        self.pos += (n + 3) & ~3
        return s

    def values(self, nc_type, n):
        dtype, size = _TYPES[nc_type]
        raw = self.buf[self.pos:self.pos + n * size]
        self.pos += (n * size + 3) & ~3
        if nc_type == 2:
            return raw.decode('utf-8', 'replace')
        return np.frombuffer(raw, dtype=dtype)

    def attrs(self):
        tag, n = self.int(), self.int()
        out = {}
        if tag == 0:
            return out
        assert tag == 0x0C, tag
        for _ in range(n):
            key = self.name()
            nc_type = self.int()
            count = self.int()
            out[key] = self.values(nc_type, count)
        return out


def read(path):
    """Return (dims, global_attrs, variables) where variables maps name -> (dim_names, attrs, ndarray)."""
    with open(path, 'rb') as f:
        buf = f.read()
    assert buf[:3] == b'CDF', 'not a netCDF-3 file'
    version = buf[3]
    r = _Reader(buf)
    r.pos = 4
    r.int()  # numrecs (unused: no record variables expected)

    dims = []
    tag, n = r.int(), r.int()
    if tag:
        assert tag == 0x0A
        for _ in range(n):
            dims.append((r.name(), r.int()))

    gattrs = r.attrs()

    variables = {}
    tag, n = r.int(), r.int()
    if tag:
        assert tag == 0x0B
        for _ in range(n):
            vname = r.name()
            ndims = r.int()
            dimids = [r.int() for _ in range(ndims)]
            vattrs = r.attrs()
            nc_type = r.int()
            r.int()  # vsize
            begin = r.int64() if version == 2 else r.int()
            shape = tuple(dims[d][1] for d in dimids)
            if any(s == 0 for s in shape):
                raise NotImplementedError('record variables are not supported')
            dtype, size = _TYPES[nc_type]
            count = int(np.prod(shape)) if shape else 1
            data = np.frombuffer(buf, dtype=dtype, count=count, offset=begin).reshape(shape)
            variables[vname] = ([dims[d][0] for d in dimids], vattrs, data)
    return dims, gattrs, variables


if __name__ == '__main__':
    import sys
    d, g, v = read(sys.argv[1])
    print('dims', d)
    for k, (dn, a, arr) in v.items():
        print(k, dn, arr.dtype, arr.shape, float(arr.min()), float(arr.max()), {kk: (vv if isinstance(vv, str) else vv.tolist()) for kk, vv in a.items() if kk in ('units', '_FillValue', 'actual_range')})
