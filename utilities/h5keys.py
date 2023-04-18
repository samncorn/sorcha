import h5py
import argparse

parser = argparse.ArgumentParser(
    prog='h5keys',
    description='Read first level key names from an hdf5 file written by pandas.',
)

parser.add_argument('filename')
args = parser.parse_args()

h5file = h5py.File(args.filename, 'r')

for group in h5file:
    print(group)