import sys
from pylab import figure, subplot, show, plot
from numpy import loadtxt

figure()
s1=subplot(611); s1.set_ylabel("K");
s2=subplot(612); s2.set_ylabel("P");
s3=subplot(613); s3.set_ylabel("K+P");
s4=subplot(614); s4.set_ylabel("2K/P");
s5=subplot(615); s5.set_ylabel("M");
s6=subplot(616); s6.set_ylabel("L");

if len(sys.argv) == 1:
    files = ['sphere-%sENERGY' % tag]
else:
    files = sys.argv[1:]

d1 = []
d2 = []
d3 = []
d4 = []
d5 = []
d6 = []
for f in files:
    X = loadtxt(f, unpack=True)

    E = X[3,0]
    print E

    d1.append(X[0]); d1.append(X[1]/E)
    d2.append(X[0]); d2.append(X[2]/E)
    d3.append(X[0]); d3.append((X[1]+X[2]) / E)
    d4.append(X[0]); d4.append(2*X[1]/X[2])
    d5.append(X[0]); d5.append(X[11])
    d6.append(X[0]); d6.append(X[6])

s1.plot(*d1)
#s1.axhline(1e63)
s2.plot(*d2)
#s2.axhline(-3e63)
s3.plot(*d3)
#s3.axhline(X[1,0] + X[2,0])
#s3.axhline(-5.665e63)
s4.plot(*d4)
#s4.axhline(-1e35)
s5.plot(*d5)
s6.plot(*d6)

show()

