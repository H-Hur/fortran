sub WL2B20s6_6 als.5 psfs.5 psf.5 sub.5 veri-

psf sub.5 WL2b20s6_6.ap psfs.5 psf.6 psfs.6 psg.6 veri- int-
all WL2B20s6_6 psg.6 psf.6 als.6 rjt.6 sub.66 veri-
sub WL2B20s6_6 als.6 psfs.6 psf.6 sub.6 veri-

psf sub.6 WL2b20s6_6.ap psfs.6 psf.7 psfs.7 psg.7 veri- int-
all WL2B20s6_6 psg.7 psf.7 als.7 rjt.7 sub.77 veri-
sub WL2B20s6_6 als.7 psfs.7 psf.7 sub.7 veri-

psf sub.7 WL2b20s6_6.ap psfs.7 WL2b20s6_6.psf WL2b20s6_6.psfs WL2b20s6_6.psg veri- int-
del als.* 
del rjt.* 
del sub.*.fits 
del psfs.* 
del psf.*.fits 
del psg.* 
del ap.* 
del coo.*

all WL2B20s6_6 WL2b20s6_6.ap WL2b20s6_6.psf als.1 rjt.1 sub.1 veri-
psort WL2b20s6_6.psfs yc
psort als.1 yc
sub WL2B20s6_6 als.1 WL2b20s6_6.psfs WL2b20s6_6.psf  WL2b20s6_6.sub.apc veri-
disp WL2b20s6_6.sub.apc 1
!./apco
mv apc.coo WL2b20s6_6.apc.coo
tvm 1 WL2b20s6_6.apc.coo rad=18.6 col=208
tvm 1 WL2b20s6_6.apc.coo rad=30
tvm 1 WL2b20s6_6.apc.coo rad=35
vi WL2b20s6_6.apc.coo
