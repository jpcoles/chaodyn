import sys
from pylab import figure, subplot, show, plot
from numpy import loadtxt, mean, ptp, amin

figure()
ss = []
ss.append(subplot(6,2,1));  ss[-1].set_ylabel(r"$T$");
ss.append(subplot(6,2,2));  ss[-1].set_ylabel(r"$U$");
ss.append(subplot(6,2,3));  ss[-1].set_ylabel(r"$T+U$");
ss.append(subplot(6,2,4));  ss[-1].set_ylabel(r"$2T/U$");
ss.append(subplot(6,2,5));  ss[-1].set_ylabel(r"$T/(r\cdot\nabla\phi)$");
ss.append(subplot(6,2,6));  ss[-1].set_ylabel(r"$p$");
ss.append(subplot(6,2,7));  ss[-1].set_ylabel(r"$R$ [kpc]");
ss.append(subplot(6,2,8));  ss[-1].set_ylabel(r"$L$");
ss.append(subplot(6,2,9));  ss[-1].set_ylabel(r"$L_x$");
ss.append(subplot(6,2,10)); ss[-1].set_ylabel(r"$L_y$");
ss.append(subplot(6,2,11)); ss[-1].set_ylabel(r"$L_z$");
ss.append(subplot(6,2,12)); ss[-1].set_ylabel(r"sdf");

if len(sys.argv) == 1:
    files = ['sphere-%sENERGY' % tag]
else:
    files = sys.argv[1:]

ds = [[] for s in ss]
for f in files:
    X = loadtxt(f, unpack=True)

    E  = X[3,0]
    CM = X[13,0]
    if CM == 0: CM = 1
    print E

    ds[0]  += [X[0], X[1]/E]
    ds[1]  += [X[0], X[2]/E]
    ds[2]  += [X[0], X[3]/E]
    ds[3]  += [X[0], 2*X[1]/X[2]]
    ds[4]  += [X[0,1:], X[12,1:]]
    ds[5]  += [X[0], X[11]]
    ds[6]  += [X[0], X[13]/CM]
    ds[7]  += [X[0], X[7]]
    ds[8]  += [X[0], X[4]]
    ds[9]  += [X[0], X[5]]
    ds[10] += [X[0], X[6]]
    #ds[11] += [X[0,1:], X[3,1:]/X[3,:-1]]
    #ds[11] += [X[0,1:], X[1,:-1]/X[1,1:]]
    ds[11] += [X[0,1:], X[2,:-1]/X[2,1:]]

for s,d in zip(ss,ds):
    s.plot(*d)

show()

