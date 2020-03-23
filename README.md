To running RTM demo, you need two file:

1 Acoustic_wave.m    

This is a microseismic forward modeling code with acoustic wave. Finiite difference mathod is used to solve the wave equation.
V1 is the input velocity model. The default siza is 300* 300, which can be change.
There are the meaning of parameters:
x0: horizontal location of source
z0: vertical location of source
dx,dz: intervel of x,z
dt: time intervel
T: length of record

You need to make sure that V*dt>dx/4 to run this code.

Seis is the output of the seismic record, you can plot it by colunm or use imagesc to see all signal.


2 RTM.m    

You need a input velocity model, seismic signal and channel number to run this code.
seismic signal is the output of Acoustic_wave.m and you should use an estimated velocity model as input. For synthetic data, you can directly use same velocity model as the forward modeling part.
Channel is an input paremeter of the active channel, if you want to use all data in Seis, just use 1:n as input. N should be the totol column of Seis.

