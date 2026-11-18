import h5py

f = h5py.File('dca.hdf5', 'r')
chi = f['step_1/CT-AUX-SOLVER-functions/spin-ZZ-chi-accumulated'][:]
chi_static = chi[0,0,...]
cc = chi_static.reshape(4,5,4,5)

chi = sum(cc, axis=(1, 3))
print(chi)


print(cc)
