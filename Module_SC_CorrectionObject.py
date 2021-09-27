# This script is for creating the module of the object of intepolations of strong coupling beta,
# Tc, effective mass, Fermi velocity. The data sheet comes form table. II of
# PRB 101.024517

# This is version 0.0 (v.0.0) of the intepolation object (class)

# author: quang. zhang (github@timohyva)


# zeta3 = 1.2020569;


def region_judgement(p):
    if (p >= 0) and (p < 2):
         pk=0;pk1=2;k=0
    if (p >= 2) and (p < 4):
         pk=2;pk1=4;k=1
    if (p >= 4) and (p < 6):
         pk=4;pk1=6;k=2
    if (p >= 6) and (p < 8):
         pk=6;pk1=8;k=3
    if (p >= 8) and (p < 10):
         pk=8;pk1=10;k=4
    if (p >= 10) and (p < 12):
         pk=10;pk1=12;k=5
    if (p >= 12) and (p < 14):
         pk=12;pk1=14;k=6
    if (p >= 14) and (p < 16):
         pk=14;pk1=16;k=7
    if (p >= 16) and (p < 18):
         pk=16;pk1=18;k=8     
    if (p >= 18) and (p < 20):
         pk=18;pk1=20;k=9
    if (p >= 20) and (p < 22):
         pk=20;pk1=22;k=10
    if (p >= 22) and (p < 24):
         pk=22;pk1=24;k=11
    if (p >= 24) and (p < 26):
         pk=24;pk1=26;k=12
    if (p >= 26) and (p < 28):
         pk=26;pk1=28;k=13
    if (p >= 28) and (p < 30):
         pk=28;pk1=30;k=14     
    if (p >= 30) and (p < 32):
         pk=30;pk1=32;k=15
    if (p >= 32) and (p < 34):
         pk=32;pk1=34;k=16

    return [pk, pk1, k] #return list


def intepolation_presure_SC(pk,pk1,p,fk,fk1):
    K = (fk1-fk)/(pk1-pk);print("K =",K, "pk1 =",pk1, "pk =",pk, "fk1 =",fk1, "fk =",fk)
    fip = K * (p-pk) + fk;print("fip =",fip)

    # return the cip or Tcp
    return fip
         
class BETA:

     # class of beta object
    

     def __init__(self,name):
         self.name = name
         print(" beta oject is crated ")

     def c1_function(self,pressure,pk,pk1,k):
         c1 = [-0.0098, -0.0127, -0.0155, -0.0181, -0.0207, -0.0231, -0.0254, -0.0275, -0.0295, -0.0314, -0.0330, -0.0345, -0.0358, -0.0370, -0.0381, -0.0391, -0.0402, -0.0413]
         self.c1p =  intepolation_presure_SC(pk,pk1,pressure,c1[k],c1[k+1])

     def c2_function(self,pressure,pk,pk1,k):
         c2 = [-0.0419, -0.0490, -0.0562, -0.0636, -0.0711, -0.0786, -0.0861, -0.0936, -0.1011, -0.1086, -0.1160, -0.1233, -0.1306, -0.1378, -0.1448, -0.1517, -0.1583, -0.1645]
         self.c2p =  intepolation_presure_SC(pk,pk1,pressure,c2[k],c2[k+1])

     def c3_function(self,pressure,pk,pk1,k):
         c3 = [-0.0132, -0.0161, -0.0184, -0.0202, -0.0216, -0.0226, -0.0233, -0.0239, -0.0243, -0.0247, -0.0249, -0.0252, -0.0255, -0.0258, -0.0262, -0.0265, -0.0267, -0.0268]
         self.c3p =  intepolation_presure_SC(pk,pk1,pressure,c3[k],c3[k+1])

     def c4_function(self,pressure,pk,pk1,k):
         c4 = [-0.0047, -0.0276, -0.0514, -0.0760, -0.1010, -0.1260, -0.1508, -0.1751, -0.1985, -0.2208, -0.2419, -0.2614, -0.2795, -0.2961, -0.3114, -0.3255, -0.3388, -0.3518]
         self.c4p =  intepolation_presure_SC(pk,pk1,pressure,c4[k],c4[k+1])

     def c5_function(self,pressure,pk,pk1,k):
         c5 = [-0.0899, -0.1277, -0.1602, -0.1880, -0.2119, -0.2324, -0.2503, -0.2660, -0.2801, -0.2930, -0.3051, -0.3167, -0.3280, -0.3392, -0.3502, -0.3611, -0.3717, -0.3815]
         self.c5p =  intepolation_presure_SC(pk,pk1,pressure,c5[k],c5[k+1])

     def tc_function(self,pressure,pk,pk1,k):
         Tc = [0.929, 1.181, 1.388, 1.560, 1.705, 1.828, 1.934, 2.026, 2.106, 2.177, 2.239, 2.293, 2.339, 2.378, 2.411, 2.438, 2.463, 2.486]
         self.tcp =  intepolation_presure_SC(pk,pk1,pressure,Tc[k],Tc[k+1])

     def mStar_function(self,pressure,pk,pk1,k):
         Ms = [2.80, 3.05, 3.27, 3.48, 3.68, 3.86, 4.03, 4.20, 4.37, 4.53, 4.70, 4.86, 5.02, 5.18, 5.34, 5.50, 5.66, 5.82] # in unit of helium-3 atom
         self.ms =  intepolation_presure_SC(pk,pk1,pressure,Ms[k],Ms[k+1])

     def vFermi_function(self,pressure,pk,pk1,k):
         VF = [59.03, 55.41, 52.36, 49.77, 47.56, 45.66, 44.00, 42.51, 41.17, 39.92, 38.74, 37.61, 36.53, 35.50, 34.53, 33.63, 32.85, 32.23]
         self.vf = intepolation_presure_SC(pk,pk1,pressure,VF[k],VF[k+1])
         
         



