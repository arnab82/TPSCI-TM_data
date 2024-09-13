from pyscf import gto, scf
import numpy as np
from pyblock2._pyscf.ao2mo import integrals as itg
from pyblock2.driver.core import DMRGDriver, SymmetryTypes
import numpy as np
ncore=120
ncas=38
n_elec=32#fockspace= [(5, 2), (2, 5), (3, 3), (3, 3), (3, 3)]
spin=0

orb_sym=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0,0,0,0,0,0,0,0]
h1e=np.load("ints_h1.npy")
g2e=np.load("ints_h2.npy")
ecore=np.load("ints_h0.npy")
print("NCAS = %d NCASELEC = %d" % (ncas, n_elec))
driver = DMRGDriver(scratch="./tmp", symm_type=SymmetryTypes.SU2,
                    stack_mem=4 << 30, n_threads=4)
# idx = driver.orbital_reordering(np.abs(h1e), np.abs(g2e)) #this does not change the energy much
# print('reordering = ', idx)
# h1e = h1e[idx][:, idx]
# g2e = g2e[idx][:, idx][:, :, idx][:, :, :, idx]
driver.initialize_system(n_sites=ncas, n_elec=n_elec, spin=spin, orb_sym=orb_sym)
mpo = driver.get_qc_mpo(h1e=h1e, g2e=g2e, ecore=ecore, iprint=0)
ket = driver.get_random_mps(tag="KET", bond_dim=250, nroots=1)
bond_dims = [250] * 8 + [500] * 8 + [750] * 8 + [1000] * 8
noises = [1e-4] * 8 + [1e-5] * 18 + [0]
thrds = [1e-8] * 32

energies = driver.dmrg(mpo, ket, n_sweeps=84, bond_dims=bond_dims, noises=noises,
    thrds=thrds, iprint=1, twosite_to_onesite=32)

bond_dims = [1150] * 8 
noises =  [0]
thrds = [1e-8] * 8

energies = driver.dmrg(mpo, ket, n_sweeps=144, bond_dims=bond_dims, noises=noises,
    thrds=thrds, iprint=1, twosite_to_onesite=30)
print('State-averaged MPS Variational energies = ',energies)

bond_dims = [1000] *48 + [800] *48 + [600] *48
noises = [0] * 72
thrds = [1e-8] * 2

ket_orig = driver.copy_mps(ket, tag='KET-ORIG')
ket = driver.adjust_mps(ket, dot=2)[0]
energies = driver.dmrg(mpo, ket, n_sweeps=160, bond_dims=bond_dims, noises=noises,
                     tol=0, thrds=thrds, iprint=1)
import scipy.stats

ds, dws, eners = driver.get_dmrg_results()
print(eners)
print('BOND DIMS         = ', ds[3::4])
print('Discarded Weights = ', dws[3::4])
print('Energies          = ', eners[3::4, 0])
reg = scipy.stats.linregress(dws[3::4], eners[3::4, 0])
emin, emax = min(eners[3::4, 0]), max(eners[3::4, 0])
print('DMRG energy (extrapolated) = %20.15f +/- %15.10f' %
            (reg.intercept, abs(reg.intercept - emin) / 5))
import matplotlib.pyplot as plt
from matplotlib import ticker
de = emax - emin
x_reg = np.array([0, dws[-1] + dws[3]])
ax = plt.gca()
ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1E}"))
plt.plot(x_reg, reg.intercept + reg.slope * x_reg, '--', linewidth=1, color='#426A5A')
plt.plot(dws[3::4], eners[3::4, 0], ' ', marker='s', mfc='white', mec='#426A5A', color='#426A5A', markersize=5)
plt.text(dws[3] * 0.25, emax + de * 0.1, "$E(M=\\infty) = %.6f \\pm %.6f \\mathrm{\\ Hartree}$" %
    (reg.intercept, abs(reg.intercept - emin) / 5), color='#426A5A', fontsize=12)
plt.text(dws[3] * 0.25, emax - de * 0.0, "$R^2 = %.6f$" % (reg.rvalue ** 2),
    color='#426A5A', fontsize=12)
plt.xlim((0, dws[-1] + dws[3]))
plt.ylim((reg.intercept - de * 0.1, emax + de * 0.2))
plt.xlabel("Largest Discarded Weight")
plt.ylabel("Sweep Energy (Hartree)")
plt.subplots_adjust(left=0.17, bottom=0.1, right=0.95, top=0.95)
