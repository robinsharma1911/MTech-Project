from  tkinter import *
import math
import numpy as np
import pandas as pd
window = Tk()
window.title('software')
window.geometry('2000x1000')

A11 = 0
A22 = 0
A12 = 0
A66 = 0
A16 = 0
A26 = 0
B11 = 0
B12 = 0
B22 = 0
B66 = 0
B26 = 0
B16 = 0
D11 = 0
D22 = 0
D12 = 0
D66 = 0
D16 = 0
D26 = 0
sum = 0
var = 0

orientation = []
thickness = []
Z_elements_from_top = []
fail_layer = []
list_stress_values_X = []
list_stress_values_Y = []
list_stress_values_XY = []
list_stress_values_L = []
list_stress_values_T = []
list_stress_values_LT = []
list_location = []

list_stress_values_top = []
list_stress_values_bottom = []
fail_layer_number = []

tensile_test = []



def func():
    for i in range(1):
        height = Entry(leftframe)
        height.grid()
        angle = Entry(leftframe)
        angle.grid()
        height = height.get()
        angle = angle.get()
            
        height = float(input("Enter thickness of layer"+str([i+1])+" from top:   "))
        thickness.append(height)
        angle = float(input("Enter orientation angle of layer"+str([i+1])+" in degrees:   "))
        orientation.append(angle)



def func1():
    global A11
    global A22
    global A12
    global A66
    global A16
    global A26
    global B11
    global B12
    global B22
    global B66
    global B26
    global B16
    global D11
    global D22
    global D12
    global D66
    global D16
    global D26
    global N_x
    global N_y
    global N_xy
    global M_x
    global M_y
    global M_xy

    A11,A12,A22,A16,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66 = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0

    vol_frac_fiber = float(vol_frac_fibera.get())
    young_mod_fiber = float(young_mod_fibera.get())
    young_mod_matrix = float(young_mod_matrixa.get())
    shear_mod_fiber = float(shear_mod_fibera.get())
    shear_mod_matrix = float(shear_mod_matrixa.get())
    poisson_ratio_fiber = float(poisson_ratio_fibera.get())
    poisson_ratio_matrix = float(poisson_ratio_matrixa.get())

    #  Longitudinal Modulus
    uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
    print("\nLongitudinal modulus of composite: ", uni_longi_mod)

    # Transverse Modulus
    ratio = (young_mod_fiber / young_mod_matrix)
    eta = (ratio - 1) / (ratio + 2)
    trans_mod = young_mod_matrix * (1 + (2 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
    print("Transverse Modulus value is: ", trans_mod)

    # Shear Modulus
    shear_mod_ratio = (shear_mod_fiber / shear_mod_matrix)
    eta = (shear_mod_ratio - 1) / (shear_mod_ratio + 1)
    shear_mod = shear_mod_matrix * (1 + (1 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
    print("Shear Modulus value (in Mpa) is: ", shear_mod)

    # Poisson's Ratio
    poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
    print("Major Poisson's Ratio is:   ", poisson_LT)

    uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
    poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
    poisson_TL = (young_mod_matrix / uni_longi_mod) * poisson_LT
    print("Minor Poisson's Ratio is:   ", poisson_TL)

    for i in range(1):
        # long_elastic_mod = float(input("\nEnter longitudinal elastic modulus of layer"+str([i+1])+":  "))
        # trans_elastic_mod = float(input("Enter transverse elastic modulus of layer"+str([i+1])+":  "))
        # shear_mod = float(input("Enter shear modulus of layer"+str([i+1])+":  "))
        # poisson_LT = float(input("Enter major poisson ratio of layer"+str([i+1])+":  "))
        # poisson_TL = float(input("Enter minor poisson ratio of layer"+str([i+1])+":  "))

        theta = orientation[i]
        sine = (math.sin(math.radians(theta)))
        cosine = (math.cos(math.radians(theta)))

        Q11 = uni_longi_mod / (1 - (poisson_LT * poisson_TL))
        Q22 = trans_mod / (1 - (poisson_LT * poisson_TL))
        Q12 = (poisson_LT * trans_mod) / (1 - (poisson_LT * poisson_TL))
        Q66 = shear_mod

        if (i+1) in fail_layer_number:
            Q11_bar, Q22_bar, Q12_bar, Q26_bar,Q16_bar,Q66_bar = 0,0,0,0,0,0
        else:
            Q11_bar = Q11*(cosine**4) + Q22*(sine**4) + 2*(Q12 + (2*Q66))*(sine**2)*(cosine**2)
            Q22_bar = Q11*(sine**4)+Q22*(cosine**4) + 2*(Q12 + (2*Q66))*(sine**2)*(cosine**2)
            Q12_bar = (Q11 + Q22 - (4*Q66)) * (sine**2) * (cosine**2) + Q12 * ((cosine**4) + (sine**4))
            Q66_bar = (Q11 + Q22 - (2*Q12) - (2*Q66)) * (sine**2) * (cosine**2) + Q66*((cosine**4) + (sine**4))
            Q16_bar = (Q11 + Q12 - (2*Q66)) * (cosine**3) * (sine) - (Q22 - Q12 - (2*Q66)) * (sine**3) * (cosine)
            Q26_bar = (Q11 + Q12 - (2*Q66)) * (cosine) * (sine**3) - (Q22 - Q12 - (2*Q66)) * (sine) * (cosine**3)

        A11 = A11 + thickness[i]*float(Q11_bar)
        A22 = A22 + thickness[i]*float(Q22_bar)
        A12 = A12 + thickness[i]*float(Q12_bar)
        A66 = A66 + thickness[i]*float(Q66_bar)
        A16 = A16 + thickness[i]*float(Q16_bar)
        A26 = A26 + thickness[i]*float(Q26_bar)

        B11 = B11 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q11_bar)
        B22 = B22 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q22_bar)
        B12 = B12 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q12_bar)
        B66 = B66 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q66_bar)
        B16 = B16 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q16_bar)
        B26 = B26 + (((Z_elements_from_top[i + 1] ** 2) - (Z_elements_from_top[i] ** 2)) / 2) * float(Q26_bar)

        D11 = D11 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q11_bar)
        D22 = D22 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q22_bar)
        D12 = D12 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q12_bar)
        D66 = D66 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q66_bar)
        D16 = D16 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q16_bar)
        D26 = D26 + (((Z_elements_from_top[i + 1] ** 3) - (Z_elements_from_top[i] ** 3)) / 3) * float(Q26_bar)


    A_matrix = np.array([
        [A11 , A12 , A16],
        [A12 , A22 , A26],
        [A16 , A26 , A66] ])

    B_matrix = np.array([
        [B11 , B12 , B16],
        [B12 , B22 , B26],
        [B16 , B26 , B66] ])

    D_matrix = np.array([
        [D11 , D12 , D16],
        [D12 , D22 , D26],
        [D16 , D26 , D66] ])

    print("\n [A]  Matrix:")
    print("-------------------------------------------------------------------")
    print(A_matrix)
    print("-----------------------------------------------------------------\n")

    print("\n [B]  Matrix:")
    print("-------------------------------------------------------------------")
    print(B_matrix)
    print("-----------------------------------------------------------------\n")

    print("\n [D]  Matrix:")
    print("-------------------------------------------------------------------")
    print(D_matrix)
    print("-----------------------------------------------------------------\n")


    # NM_matrix = np.matrix([
    #     [N_x],[N_y],[N_xy],[M_x],[M_y],[M_xy]
    # ])

    Stiffness_matrix = np.matrix([
        [A11 , A12 , A16 , B11 , B12 , B16],
        [A12 , A22 , A26 , B12 , B22 , B26],
        [A16 , A26 , A66 , B16 , B26 , B66],
        [B11 , B12 , B16 , D11 , D12 , D16],
        [B12 , B22 , B26 , D12 , D22 , D26],
        [B16 , B26 , B66 , D16 , D26 , D66]
    ])
    Inverse_stiffness_matrix = Stiffness_matrix.getI()    # Inverse of a matrix

    # Strain_curvature_matrix = Inverse_stiffness_matrix * NM_matrix
    # print("\n------------------------------- Strains and curvatures matrix for mid plane:-------------------------")
    # print(Strain_curvature_matrix)
    # print("----------------------------------------------------------------------------------------------------\n")

    print("\n------------------------------- Inverse Stiffness Matrix for elastic constants:-------------------------")
    print(Inverse_stiffness_matrix)
    print("----------------------------------------------------------------------------------------------------\n")


# """==================================================================================================================================
# -------------------------------------Strain Energy Release Rate----------------------------------------------"""
#
# print("Note: Calculation is valid for bend specimens with S/W = 4.\n")
# load = float(input("Enter load in KN:   "))
# spec_thickness = float(input("Enter specimen thickness in 'cm':   "))
# spec_width = float(input("Enter specimen depth(width) in 'cm':    "))
# crack_length = float(input("Enter crack length in 'cm':    \n"))
# x = crack_length / spec_width
#
# F_x = (6*(x**0.5)*(1.99-x*(1-x)*(2.15-(3.93*x) + (2.7*x*x)))) / ((1 + (2*x))*(1-x)**1.5)
#
# strain_energy_release_rate = (F_x * load)/ (spec_thickness*(spec_width**0.5))
#
# print("Strain Energy Release Rate(Kq):   ",strain_energy_release_rate)
#
#
# """------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------Plain Strain Fracture Toughness(Gic)-----------------------------------------"""
# #
# # Corr_energy = float(input("Enter Corrected Energy Value:   "))
# #
# # dA_dx = ((16*x*x)/((1-x)**2))*((-33.717)+(159.232*x) - (338.856*x*x) - (339.26*x*x*x) - (128.36*x*x*x*x)) + ((32*x)/((1-x)**3)) * (8.9 - (33.717*x) + (79.616*x*x) - (112.952*x*x*x) + (84.815*x*x*x*x) - (25.672*x*x*x*x*x))
# #
# # A = ((16*x*x)/((1-x)**2)) * (8.9 - (33.7.7*x) + (79.616*x*x) - (112.952*x*x*x) + (84.815*x*x*x*x) - (25.672*x*x*x*x*x))
# #
# # phai = ( A + 18.64 ) / dA_dx
# #
# # plain_strain_fracture_toughness = Corr_energy / (spec_thickness * spec_width * phai)
# #
# # print("\nPlain Strain Fracture Toughness Value:    ", plain_strain_fracture_toughness)
#
#
# """---------------------------------------------------------------------------------------------------------------------"""
# print("\nEnter Mechanical properties for each layer.\n")
# long_elastic_mod = float(input("Enter longitudinal elastic modulus:  "))
# trans_elastic_mod = float(input("Enter transverse elastic modulus:  "))
# shear_mod = float(input("Enter shear modulus:  "))
# poisson_LT = float(input("Enter major poisson ratio:  "))
# poisson_TL = float(input("Enter minor poisson ratio:  "))

#"))
#
# layers = int(input("Enter number of layers in laminate:   "))


leftframe = Frame(window)
leftframe.grid()
vol_frac_fiber_lbl = Label(leftframe,text="Enter volume fraction of fiber:  ")
vol_frac_fiber_lbl.grid(row = 0, column = 0)
young_mod_fiber_lbl = Label(leftframe,text="Enter Young's Modulus (in Mpa) of fiber:  ")
young_mod_fiber_lbl.grid(row = 1, column = 0)
young_mod_matrix_lbl = Label(leftframe,text="Enter Young's Modulus (in Mpa) of Matrix:  ")
young_mod_matrix_lbl.grid(row = 2, column = 0)
shear_mod_fiber_lbl = Label(leftframe,text="Enter Shear Modulus (in Mpa) of fiber:  ")
shear_mod_fiber_lbl.grid(row = 3, column = 0)
shear_mod_matrix_lbl = Label(leftframe,text="Enter Shear Modulus (in Mpa) of Matrix:  ")
shear_mod_matrix_lbl.grid(row = 4, column = 0)
poisson_ratio_fiber_lbl = Label(leftframe,text="Enter poisson's ratio of fiber:  ")
poisson_ratio_fiber_lbl.grid(row = 5, column = 0)
poisson_ratio_matrix_lbl = Label(leftframe,text="Enter poisson's ratio of matrix:  ")
poisson_ratio_matrix_lbl.grid(row = 6, column = 0)


vol_frac_fibera = Entry(leftframe)
vol_frac_fibera.grid(row = 0, column = 1)
young_mod_fibera = Entry(leftframe)
young_mod_fibera.grid(row = 1, column = 1)
young_mod_matrixa = Entry(leftframe)
young_mod_matrixa.grid(row = 2, column = 1)
shear_mod_fibera = Entry(leftframe)
shear_mod_fibera.grid(row = 3, column = 1)
shear_mod_matrixa =Entry(leftframe)
shear_mod_matrixa.grid(row = 4, column = 1)
poisson_ratio_fibera = Entry(leftframe)
poisson_ratio_fibera.grid(row = 5, column = 1)
poisson_ratio_matrixa = Entry(leftframe)
poisson_ratio_matrixa.grid(row = 6, column = 1)





output_button = Button(leftframe,text="output",command=func)
output_button.grid(row=12,column=1)



# vol_frac_fiber = float(input("Enter volume fraction of fiber:  "))
# young_mod_fiber = float(input("Enter Young's Modulus (in Mpa) of fiber:  "))
# young_mod_matrix = float(input("Enter Young's Modulus (in Mpa) of Matrix:  "))
# shear_mod_fiber = float(input("Enter Shear Modulus (in Mpa) of fiber:  "))
# shear_mod_matrix = float(input("Enter Shear Modulus (in Mpa) of Matrix:  "))
# poisson_ratio_fiber = float(input("Enter poisson's ratio of fiber:  "))
# poisson_ratio_matrix = float(input("Enter poisson's ratio of matrix: " ))


func()
# Calculation of Z values for each layer according to thickness
for j in range(0, len(thickness)):
    sum = sum + thickness[j]
median = sum / 2

for j in range(0, len(thickness)):
    if j == 0:
        Z_elements_from_top.append(-median)

    var = var + thickness[j]
    z_value = var - median
    Z_elements_from_top.append(z_value)

print("\nOrientation angle from top to bottom: ",orientation)
print("Thickness of each layer from top to bottom: ",thickness)
print("Z elements required for [B] and [D] matrix:  ",Z_elements_from_top)

func1()



window.mainloop()




