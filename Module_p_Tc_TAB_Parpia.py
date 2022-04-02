





import matplotlib.pyplot as plt
import numpy as np

import csv

##########################################################
# import the Tc csv data of Parpia experment
# with open('Tc_Parpia.csv', newline='') as f:
#      reader = csv.reader(f)
#      data = list(reader)
list1 = np.genfromtxt("Tc_Parpia.csv",delimiter=",")

# print("\n",data,"\n")


# list1 = list(zip(*data))
# print(list1)
# print(list(list1[0]))
# print(list1[1])

Tc_array = np.array([])
for ii in list1[:,0]:
    Tc_array = np.append(Tc_array,ii)

print("\nParpia's Tc-line temperture is \n", Tc_array)  
    
p_Tc_array = np.array([])
for ii in list1[:,1]:
    p_Tc_array = np.append(p_Tc_array,ii)    

print("\nParpia's Tc-line pressure is \n", p_Tc_array)

###########################################################
###########################################################

# import the Tc csv data of Parpia experment
# with open('TAB_Parpia.csv', newline='') as ff:
#      reader1 = csv.reader(ff)
#      data1 = list(reader1)

list11 = np.genfromtxt('TAB_Parpia.csv',delimiter=',')


# print("\n\n",data1,"\n")


# list11 = list(zip(*data1))
# print("\n",list11)
# print(list(list1[0]))
# print("\n",list11[:,1])

TAB_array = np.array([])
for ii in list11[:,0]:
    TAB_array = np.append(TAB_array, ii)
TAB_array = np.sort(TAB_array)[::-1] # sort as decreasing    

print("\nParpia's TAB-line temperture is \n", TAB_array)  
    
p_TAB_array = np.array([])
for jj in list11[:,1]:
    p_TAB_array = np.append(p_TAB_array, jj)
p_TAB_array = np.sort(p_TAB_array)    

print("\nParpia's TAB-line pressure is \n", p_TAB_array)

#############################################################
#############################################################

p = np.arange(22.,29.1,0.1) # bar

print(" \n Pre looks like ", p)

TAB = np.interp(p, p_TAB_array, TAB_array)

Tc = np.interp(p, p_Tc_array, Tc_array)

tAB = TAB/Tc

print("\n TAB looks like ", TAB, "\n Tc looks like ", Tc, "\n tAB looks like ", tAB)



# fig, ax = plt.subplots(1,1)

# ax.scatter(Tc_array, p_Tc_array)
# ax.plot(Tc, p)
# ax.scatter(TAB_array, p_TAB_array)
# ax.plot(TAB, p)
# ax.set_xlabel(r"$T/mK$")
# ax.set_ylabel(r"$p/bar$")
# ax.set_title(r"digtized experimantal TAB, Tc from Parpia manuscript.pdf")

# ax.legend(labels=("Tc-Parpia","TAB-Parpia"), loc="lower right")

# fig.savefig("digtized_experimantal_TAB_Tc_Parpia_manuscript.pdf")

# plt.show()


