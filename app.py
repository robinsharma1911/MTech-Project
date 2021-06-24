from flask import Flask , render_template,request

import math
import numpy as np
import pandas as pd
  
app = Flask(__name__) #creating the Flask class object   
 
@app.route('/', methods = ['POST', 'GET']) 
def home():  
    if request.method == 'POST':
        #program   
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

#______________________GUI______________________________________-

        vol_frac_fiber = float(request.form['vol_frac_fiber'])
        young_mod_fiber = float(request.form['young_mod_fiber'])
        young_mod_matrix = float(request.form['young_mod_matrix'])
        shear_mod_fiber = float(request.form['shear_mod_fiber'])
        shear_mod_matrix = float(request.form['shear_mod_matrix'])
        poisson_ratio_fiber = float(request.form['poisson_ratio_fiber'])
        poisson_ratio_matrix = float(request.form['poisson_ratio_matrix'])


        if request.form['height1'] and request.form['angle1']:
            height1 = float(request.form['height1'])
            angle1 = float(request.form['angle1'])
            thickness.append(height1)
            orientation.append(angle1)

        if request.form['height2'] and request.form['angle2']:
            height2 = float(request.form['height2'])
            angle2 = float(request.form['angle2'])
            thickness.append(height2)
            orientation.append(angle2)

        if request.form['height3']  and request.form['angle3']:
            height3 = float(request.form['height3'])
            angle3 = float(request.form['angle3'])
            thickness.append(height3)
            orientation.append(angle3)
        if request.form['height4']  and request.form['angle4']:
            height4 = float(request.form['height4'])
            angle4 = float(request.form['angle4'])
            thickness.append(height4)
            orientation.append(angle4)

        if request.form['height5'] and request.form['angle5']:
            height5 = float(request.form['height5'])
            angle5 = float(request.form['angle5'])
            thickness.append(height5)
            orientation.append(angle5)
        if request.form['height6'] and request.form['angle6']:
            height6 = float(request.form['height6'])
            angle6 = float(request.form['angle6'])
            thickness.append(height6)
            orientation.append(angle6)
      
#__________________________GUI_END_____________________________________

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

     

        A11,A12,A22,A16,A26,A66,B11,B12,B16,B22,B26,B66,D11,D12,D16,D22,D26,D66 = 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0


        #  Longitudinal Modulus
        uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
        # print("\nLongitudinal modulus of composite: ", uni_longi_mod)

        # Transverse Modulus
        ratio = (young_mod_fiber / young_mod_matrix)
        eta = (ratio - 1) / (ratio + 2)
        trans_mod = young_mod_matrix * (1 + (2 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
        # print("Transverse Modulus value is: ", trans_mod)

        # Shear Modulus
        shear_mod_ratio = (shear_mod_fiber / shear_mod_matrix)
        eta = (shear_mod_ratio - 1) / (shear_mod_ratio + 1)
        shear_mod = shear_mod_matrix * (1 + (1 * eta * vol_frac_fiber)) / (1 - (eta * vol_frac_fiber))
        # print("Shear Modulus value (in Mpa) is: ", shear_mod)

        # Poisson's Ratio
        poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
        # print("Major Poisson's Ratio is:   ", poisson_LT)

        uni_longi_mod = young_mod_matrix * (vol_frac_fiber * ((young_mod_fiber / young_mod_matrix) - 1) + 1)
        poisson_LT = (poisson_ratio_fiber * vol_frac_fiber) + (poisson_ratio_matrix * (1 - vol_frac_fiber))
        poisson_TL = (young_mod_matrix / uni_longi_mod) * poisson_LT
        # print("Minor Poisson's Ratio is:   ", poisson_TL)

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


        A_matrix = np.matrix([
            [A11 , A12 , A16],
            [A12 , A22 , A26],
            [A16 , A26 , A66] ])

        B_matrix = np.matrix([
            [B11 , B12 , B16],
            [B12 , B22 , B26],
            [B16 , B26 , B66] ])

        D_matrix = np.matrix([
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
        Inverse_A_matrix = A_matrix.getI()
        tensile_modulus = 1 / (sum * Inverse_A_matrix[0,0])
        print("Sum:  ",sum)
        # Strain_curvature_matrix = Inverse_stiffness_matrix * NM_matrix
        # print("\n------------------------------- Strains and curvatures matrix for mid plane:-------------------------")
        # print(Strain_curvature_matrix)
        # print("----------------------------------------------------------------------------------------------------\n")

        print("\n------------------------------- Inverse Stiffness Matrix for elastic constants:-------------------------")
        print("\nLongitudinal modulus of composite: ", uni_longi_mod)
        print("Transverse Modulus value is: ", trans_mod)
        print("Shear Modulus value (in Mpa) is: ", shear_mod)
        print("Major Poisson's Ratio is:   ", poisson_LT)
        print("Minor Poisson's Ratio is:   ", poisson_TL)
        print("tensile_modulus Modulus value is: ", tensile_modulus)



        
        print(Inverse_stiffness_matrix)

        print("----------------------------------------------------------------------------------------------------\n")


        return render_template('index.html',matrix1 = Inverse_stiffness_matrix,uni_longi_mod =uni_longi_mod,trans_mod = 
        trans_mod, shear_mod = shear_mod, poisson_LT = poisson_LT, poisson_TL = poisson_TL, tensile_modulus = tensile_modulus ) 

    return render_template('index.html') 
  
if __name__ =='__main__':  
    app.run(debug = True, port="8000" )  