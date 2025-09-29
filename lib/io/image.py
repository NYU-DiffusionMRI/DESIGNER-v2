'''
Load and save images in MRtrix format.

Copyright (c) 2017 - Daan Christiaens (daan.christiaens@gmail.com)
'''

import numpy as np
import sys
import copy
import gzip


_dtdict = {'Int8': '|i1', 'UInt8': '|u1', 'Int16': '=i2', 'UInt16': '=u2', 'Int16LE': '<i2', 'UInt16LE': '<u2', 'Int16BE': '>i2', 'UInt16BE': '>u2', 'Int32': '=i4', 'UInt32': '=u4', 'Int32LE': '<i4', 'UInt32LE': '<u4', 'Int32BE': '>i4', 'UInt32BE': '>u4', 'Float32': '=f4', 'Float32LE': '<f4', 'Float32BE': '>f4', 'Float64': '=f8', 'Float64LE': '<f8', 'Float64BE': '>f8', 'CFloat32': '=c8', 'CFloat32LE': '<c8', 'CFloat32BE': '>c8', 'CFloat64': '=c16', 'CFloat64LE': '<c16', 'CFloat64BE': '>c16'}
_dtdict_inv = {v: k for k, v in _dtdict.items()}


class Image (object):
    '''
    Lightweight wrapper class that stores MRtrix images in numpy ndarray objects.

    Class attributes:
      data:        np.ndarray that stores the image data with its datatype and shape
      vox:         image voxel size
      transform:   image transformation matrix
      grad:        image gradient table
      comments:    header comments

    The class also exposes these data attributes:
      shape, ndim, dtype, size, strides, nbytes
    '''

    def __init__(self, data=None, vox=(),
                       transform=np.eye(4),
                       grad=None, comments=None):
        self.data = data
        self.vox = vox
        self.transform = transform
        self.grad = grad
        self.comments = comments if comments is not None else []


    _array_attr = ['shape', 'ndim', 'dtype', 'size', 'strides', 'nbytes']

    def __getattr__(self, attribute):
        if attribute in self._array_attr:
            if self.data is None:
                raise AttributeError('Image data not set.')
            return getattr(self.data, attribute)


    def __copy__(self):
        return Image(self.data.copy(), self.vox,
                self.transform.copy(), self.grad.copy(),
                copy.copy(self.comments))

    def copy(self):
        ''' Copy image in memory. '''
        return self.__copy__()


    @classmethod
    def empty_as(cls, hdr):
        ''' Create empty image based off the header of another image. '''
        return cls(None, hdr.vox, hdr.transform, hdr.grad, hdr.comments)


    @property
    def vox(self):
        ''' Image voxel size. '''
        if self.data is None:
            return self._vox
        else:
            n = min(len(self._vox), self.ndim)
            return self._vox[:n] + (self.ndim - n) * (1.,)

    @vox.setter
    def vox(self, v):
        ''' Set voxel size. '''
        self._vox = tuple(map(float, v))


    @property
    def nvox(self):
        ''' Get number of voxels in the image. '''
        if self.data is None:
            return 0
        else:
            return np.prod(self.shape[:3])

    def load(self, filename, header_only=False, memmap=False):
        ''' Load MRtrix .mif or .mif.gz file. '''
        is_gz = False
        if memmap:
            if is_gz:
                raise IOError('memory mapping of .mif.gz files is not supported')
            if header_only:
                raise IOError('memory mapping and header_only are mutually exclusive')
        if filename.endswith('.mif'):
            ofun = open
        elif filename.endswith('.mif.gz'):
            ofun = gzip.open
            is_gz = True
        else:
            raise IOError('file extension not supported: ' + str(filename))
        # read image header
        with ofun(filename, 'rt', encoding='latin-1') as f:
            fl = ''
            tr_count = 0
            while fl != 'END':
                fl = f.readline().strip()
                if fl.startswith('dim'):
                    imsize = tuple(map(int, fl.split(':')[1].strip().split(',')))
                elif fl.startswith('vox'):
                    self.vox = fl.split(':')[1].strip().split(',')
                elif fl.startswith('layout'):
                    layout = fl.split(':')[1].strip().split(',')
                elif fl.startswith('datatype'):
                    dtstr = fl.split(':')[1].strip()
                    dt = np.dtype(_dtdict.get(dtstr, 'u1'))
                elif fl.startswith('file'):
                    offset = int(fl.split('.')[1].strip());
                elif fl.startswith('transform'):
                    self.transform[tr_count,:] = np.array(fl.split(':')[1].strip().split(','), dtype=float)
                    tr_count = tr_count + 1
                elif fl.startswith('labels'):
                    self.labels = fl.split(':')[1].strip()
                elif fl.startswith('units'):
                    self.units = fl.split(':')[1].strip()
                elif fl.startswith('comments'):
                    self.comments.append(fl[9:].strip())
                elif fl.startswith('dw_scheme'):
                    gbrow = np.array(fl.split(':')[1].strip().split(','), dtype=float)
                    if self.grad is None:
                        self.grad = gbrow
                    else:
                        self.grad = np.vstack([self.grad, gbrow])
        if not header_only:
            # read image data
            with ofun(filename, 'rb') as f:
                if is_gz:
                    f.seek(offset, 0)
                    buf = f.read()
                    image = np.frombuffer(buffer=buf, dtype=dt)
                elif memmap:
                    image = np.memmap(f, mode='r', dtype=dt, offset=offset)  # offset required in memmap
                else:
                    f.seek(offset, 0)
                    image = np.fromfile(file=f, dtype=dt)
                if (dtstr == 'Bit'):
                    image = np.unpackbits(image)
                s, o = self._layout_to_strides(layout, imsize, dt)
                self.data = np.ndarray(shape=imsize, dtype=dt, buffer=image, strides=s, offset=o)
        return self


    def save(self, filename):
        ''' Save image to MRtix .mif file. '''
        if self.data is None:
            raise RuntimeError('Image data not set.')
        if not filename.endswith('.mif'):
            raise IOError('only .mif file type supported for writing')
        # write image header
        with open(filename, 'w', encoding='latin-1') as f:
            f.write('mrtrix image\n')
            f.write('dim: ' + self._to_csv(self.shape) + '\n');
            f.write('vox: ' + self._to_csv(self.vox) + '\n')
            f.write('layout: ' + self._to_csv(self.layout, precision='%s') + '\n')
            f.write('datatype: ' + _dtdict_inv[self.dtype.descr[0][1]] + '\n')
            f.write(self._to_csv2D(self.transform[:3], 'transform: '))
            if self.labels is not None:
                f.write('labels: ' + self._to_csv(self.labels, precision='%s') + '\n')
            if self.units is not None:
                f.write('units: ' + self._to_csv(self.units, precision='%s') + '\n')
            if self.comments:
                f.write('comments: %s\n' * len(self.comments) % tuple(self.comments))
            if self.grad is not None:
                f.write(self._to_csv2D(self.grad, line_prefix='dw_scheme: '))
            f.flush()
            offset = f.tell() + 13
            offset += int(np.floor(np.log10(offset))) + 1
            f.write('file: . {:d}\n'.format(offset))
            f.write('END\n')
            f.flush()
        # write image data
        with open(filename, 'ab') as f:
            self.data.ravel(order='K').tofile(f)
        return self


    def _layout_to_strides(self, layout, size, dtype):
        strides = [0 for l in layout]
        stride, offset = int(dtype.itemsize), 0
        for dim in sorted(range(len(layout)), key=lambda k: int(layout[k][1:])):
            if layout[dim][0] is '-':
                strides[dim] = -stride
                offset += (size[dim]-1) * stride
            else:
                strides[dim] = stride
            stride *= size[dim]
        return strides, offset


    @property
    def layout(self):
        ''' Data layout in output file.
        Currently, only positive strides are supported due to numpy limitations.
        '''
        #return tuple(('-' if self.strides[s]<0 else '+') + str(s) for s in np.argsort(np.argsort(np.abs(self.strides))))
        return tuple('+'+str(s) for s in np.argsort(np.argsort(np.abs(self.strides))))


    def _to_csv(self, a, precision=None):
        if isinstance(precision, str):
            fmt = ','.join([precision]*len(a))
            return fmt % tuple(a)
        if precision is None:
            if isinstance(a, np.ndarray):
                precision = np.finfo(a.dtype).precision
            else:
                precision = 16
        elif not isinstance(precision, int) or precision < 0:
            raise ValueError('precision needs to be non-negative integer, got ' + str(precision))
        return ','.join(['%.'+str(precision)+'g']*len(a)) % tuple(a)


    def _to_csv2D(self, a, line_prefix='', postfix='\n', precision=None):
        if not isinstance(a, np.ndarray):
            a = np.asanyarray(a)
        if a.ndim != 2:
            raise ValueError('require 2D array, got array of shape ' + str(a.shape))
        if precision is None:
            if isinstance(a, np.ndarray):
                precision = np.finfo(a.dtype).precision
            else:
                precision = 16
        elif not isinstance(precision, int) or precision < 0:
            raise ValueError('precision needs to be non-negative integer, got ' + str(precision))
        fmt = [','.join(['%.'+str(precision)+'g'] * a.shape[1])] * a.shape[0]
        fmt = line_prefix+(postfix+line_prefix).join(fmt)
        return fmt % tuple(a.ravel()) + postfix


    def __str__(self):
        out = 'mrtrix image:'
        if self.data is not None:
            out += '\n  dimensions: ' + self._to_csv(self.shape) + '\n'
            out += '  voxel size: ' + self._to_csv(self.vox) + '\n'
            out += '  datatype: ' + _dtdict_inv[self.dtype.descr[0][1]] + '\n'
            tx, ty, tz = str(self.transform[:3,:]).split('\n')
            out += '  transform: ' + tx + '\n'
            out += '             ' + ty + '\n'
            out += '             ' + tz
            if self.grad is not None:
                out += '\n  gradient table: {:d} x {:d}'.format(*self.grad.shape)
        else:
            out += ' empty'
        return out


    def __iter__(self):
        self._pos = 0
        return self


    def __next__(self):
        if self._pos >= self.nvox:
            raise StopIteration
        out = self.data[self._pos % self.shape[0],
                        self._pos//self.shape[0] % self.shape[1],
                        self._pos//np.prod(self.shape[:2]) % self.shape[2]]
        self._pos += 1
        return out


def load_mrtrix(filename, **kwargs):
    ''' Load image in mrtrix format. '''
    img = Image()
    img.load(filename, **kwargs)
    return img


def save_mrtrix(filename, image):
    ''' Save image in mrtrix format. '''
    image.save(filename)




