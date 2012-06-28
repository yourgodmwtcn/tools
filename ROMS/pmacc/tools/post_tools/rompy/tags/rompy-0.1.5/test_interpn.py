#!/usr/bin/env python
import numpy as np

from rompy import utils

a = np.linspace(-10.,10.,200)
x = np.tile(a.reshape([1,1,len(a)]),[len(a),len(a),1])
y = np.tile(a.reshape([1,len(a),1]),[len(a),1,len(a)])
z = np.tile(a.reshape([len(a),1,1]),[1,len(a),len(a)])

#(x, y, z) = utils.meshgrid(a,a,a)

# print('x')
# print(x[0,0,:])
# print('y')
# print(y[0,:,0])
# print('z')
# print(z[:,0,0])

#ai = np.linspace(-10.,10,100)
#(xi,yi,zi) = utils.meshgrid(ai,np.array([0.0]),ai)

ai = np.array([-5., 0., 5.])
xi = np.tile(ai.reshape([1,1,len(ai)]),[len(ai),len(ai),1])
yi = np.tile(ai.reshape([1,len(ai),1]),[len(ai),1,len(ai)])
zi = np.tile(ai.reshape([len(ai),1,1]),[1,len(ai),len(ai)])

#data = np.sin(x) + np.cos(y) + z
#di = np.sin(xi) + np.cos(yi) + zi

data = x + 2*y + 3*z
di = xi + 2*yi + 3*zi

#utils.interp_3d(x,y,z,data,xi,yi,zi)

#print(data)

#print(di)

d_interp1 = utils.interp_3d(x,y,z,data,ai,ai,ai)
d_interp2 = utils.interp_3d(x,y,z,data,xi,yi,zi)

print(di - d_interp1)
print(di - d_interp2)
