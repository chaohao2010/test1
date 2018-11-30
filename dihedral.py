#!/usr/bin/env python
# coding: utf-8

# # Calculate the dihedral of four atoms

# In[5]:


import math

def Dihedral_Calculation(atomA, atomB, atomC, atomD):
    #each atom is list which contains the x,y,z coordinates
    i = 0
    pi = math.pi
    x = [atomA[0], atomB[0], atomC[0], atomD[0]]
    y = [atomA[1], atomB[1], atomC[1], atomD[1]]
    z = [atomA[2], atomB[2], atomC[2], atomD[2]]
    x_vect = []; y_vect = []; z_vect = []
    x_Nvect = []; y_Nvect = []; z_Nvect = []

    #calculate 3 vectors
    for i in range(3):
        x_vect.append(x[i+1] - x[i])
        y_vect.append(y[i+1] - y[i])
        z_vect.append(z[i+1] - z[i])
    #calculate 2 normal vectors
    for i in range(2):
        x_Nvect.append(y_vect[i] * z_vect[i+1] - z_vect[i] * y_vect[i+1])
        y_Nvect.append(z_vect[i] * x_vect[i+1] - x_vect[i] * z_vect[i+1])
        z_Nvect.append(x_vect[i] * y_vect[i+1] - y_vect[i] * x_vect[i+1])

    #calculate cos of the dihedrals
    cos_dihedral = ((x_Nvect[0] * x_Nvect[1] + y_Nvect[0] * y_Nvect[1] + z_Nvect[0] * z_Nvect[1]) /
                    (x_Nvect[0] * x_Nvect[0] + y_Nvect[0] * y_Nvect[0] + z_Nvect[0] * z_Nvect[0]) ** 0.5 /
                    (x_Nvect[1] * x_Nvect[1] + y_Nvect[1] * y_Nvect[1] + z_Nvect[1] * z_Nvect[1]) ** 0.5)

    if cos_dihedral == 1:
        dihedral = 0
    elif cos_dihedral == -1:
        dihedral = pi
    else:
        dihedral = math.acos(cos_dihedral)

    #covert the radian measure to degree measure
    dihedral = dihedral / pi * 180.0

    #check the dihedral is positive or negetive
    x_cross = y_vect[0] * z_vect[2] - z_vect[0] * y_vect[2]
    y_cross = z_vect[0] * x_vect[2] - x_vect[0] * z_vect[2]
    z_cross = x_vect[0] * y_vect[2] - y_vect[0] * x_vect[2]

    if (x_cross * x_vect[1] + y_cross * y_vect[1] + z_cross * z_vect[1]) > 0:
        dihedral = -1.0 * dihedral
    return dihedral

def main():
    atomA = [29.053, 40.319, 5.276]
    atomB = [30.195, 39.477, 4.840]
    atomC = [29.684, 38.597, 3.706]
    atomD = [30.535, 37.881, 2.972]
    dihedral = Dihedral_Calculation(atomA, atomB, atomC, atomD)
    print(dihedral)

if __name__ == '__main__':
    main()


# In[ ]:




