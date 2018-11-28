import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from scipy.ndimage import rotate

'''
dat=fits.open('psf1.fits')
#21theta 21phi 5energies 1defocus 256x256 pixels
image=dat[0].data[10][20][1][0]
plt.imshow(image,origin='lower')
plt.show()
sys.exit()

fig = plt.figure(figsize=(21,21))
fig.subplots_adjust(hspace=0, wspace=0)
k=1
for i in range(21):
    for j in range(21):
        print(i,j)
        ax=fig.add_subplot(21,21,k)
        k+=1
        image=dat[0].data[i][j][0][0]
        plt.imshow(image,origin='lower')
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        #ax.text(125, 200, str((i,j)),fontsize=10, color='w', ha='center')
plt.tight_layout()
plt.savefig('psf.png',format='png')
dat.close()
'''

theta=9
phi=90
roll=270
phi_eq=phi

elev=theta*np.sin(phi_eq*np.pi/180.)
azim=theta*np.cos(phi_eq*np.pi/180.)
ii=10+int(round(elev))
jj=10+int(round(azim))

dat=fits.open('psf1_rebinned.fits')
#21theta 21phi 5energies 1defocus 256x256 pixels
image=dat[0].data[ii][jj][1][0]
#image2=dat[0].data[ii][jj][3][0]

plt.imshow(image,origin='lower')
plt.show()
#plt.imshow(image2,origin='lower')
#plt.show()
'''
#data_orig = misc.face()
data_orig = image
x0,y0 = 128,128 # left eye; (xrot,yrot) should point there

def rot(image, xy, angle):
    im_rot = rotate(image,angle) 
    org_center = (np.array(image.shape[:2][::-1])-1)/2.
    rot_center = (np.array(im_rot.shape[:2][::-1])-1)/2.
    org = xy-org_center
    a = np.deg2rad(angle)
    new = np.array([org[0]*np.cos(a) + org[1]*np.sin(a),
            -org[0]*np.sin(a) + org[1]*np.cos(a) ])
    return im_rot, new+rot_center

plt.imshow(data_orig,origin='lower')
plt.show()

data_rot, (x1,y1) = rot(data_orig, np.array([x0,y0]), roll)
print(data_rot.shape)
'''
im_rot = rotate(image,roll) 
plt.imshow(im_rot,origin='lower')
plt.show()

'''
plt.imshow(image,origin='lower')
plt.show()

larger_image=np.zeros((400,400))
n=0
for jj in range(72,328):
	m=0
	for ii in range(72,328):
		larger_image[ii][jj]=image[m][n]
		m=m+1
	n=n+1

plt.imshow(larger_image,origin='lower')
plt.show()

alfa=roll
r=np.array([[np.cos(alfa*np.pi/180.),-np.sin(alfa*np.pi/180.)],[np.sin(alfa*np.pi/180.),np.cos(alfa*np.pi/180.)]])
newimage=np.zeros_like(larger_image)

plt.imshow(newimage,origin='lower')
plt.show()

n=0
for jj in range(72,328):
	m=0
	for ii in range(72,328):
		x=r[0][0]*(ii-127.5)+r[0][1]*(127.5-jj)+127.5
		y=r[1][0]*(ii-127.5)+r[1][1]*(jj-127.5)+127.5
		#newimage[ii][jj]=image[m][n]
		if (int(round(x)) < 256) and (int(round(y)) < 256):
			newimage[ii][jj]=image[int(round(x))][int(round(y))]
		m=m+1
	n=n+1


for i in range(len(larger_image)):
	for j in range(len(larger_image[i])):
		x=r[0][0]*(i-200)+r[0][1]*(200-j)+200
		y=r[1][0]*(i-200)+r[1][1]*(j-200)+200
		#if (int(round(x)) < 256) and (int(round(y)) < 256):
		#	print(int(round(x)),int(round(y)))
		newimage[i][j]=larger_image[int(round(x))][int(round(y))]

plt.imshow(newimage,origin='lower')
plt.show()


'''


'''
for i, angle in enumerate([66,-32,90]):
    data_rot, (x1,y1) = rot(data_orig, np.array([x0,y0]), angle)
    axes.flatten()[i+1].imshow(data_rot,origin='lower')
    axes.flatten()[i+1].scatter(x1,y1,c="r" )
    axes.flatten()[i+1].set_title("Rotation: {}deg".format(angle))

plt.show()
'''

