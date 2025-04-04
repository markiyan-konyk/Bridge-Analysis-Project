from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
# Length of Bridge = 1250mm
# Length of Train = 960mm
# Distance to Travel = 290mm

V1 = 0
M1 = 0

V2 = 0
M2 = 0

def findI(h,b):
    # h is the height of the supports
    # b is the width of the flange
    # 1.26 is the thicknes of everything
    A1 = 1.26*h
    A2 = 1.26*b
    th = 1.26
    ybar = (((h/2)*(2*A1))+((h+0.63)*(A2)))/((2*A1)+(A2))
    
    
    I1 = (th*(h**3))/12
    I2 = ((th**3)*b)/12
    yd1 = ybar - (h/2)
    yd2 = (h+1.26) - ybar
    I = (2*I1) + I2 + (A1*(yd1**2)) + (A2*(yd2**2))
    I = I/1000000
    
    return I


def graphInertia():
    all_I = []
    h_axis = []
    w_axis = []
    
    for i in range(60,201,20):
        h_axis.append(i)
        w_axis.append(300-i)
            
    for i in range(len(h_axis)):
        all_I.append(findI(h_axis[i],w_axis[i]))
        
    ax = plt.axes(projection ='3d')

    ax.plot3D(h_axis,w_axis,all_I)
    ax.set_xlabel("Height (mm)")
    ax.set_ylabel("Width (mm)")
    ax.set_zlabel("Moment of Inertia ( *10^6 mm4)")
    plt.show()
    
    
def LoadCase1(F):  ## F is the weight of each wheel
    global V1
    global M1
    
    P = [0,0,0,0,0,0]
    P[0],P[1],P[2],P[3],P[4],P[5] = -(F/6),-(F/6),-(F/6),-(F/6),-(F/6),-(F/6)
    print(sum(P))
    print(P)

    # The locomotive wagon (heaviest) is on the right side for these calculations

    all_fbd = {}
    
    # Make a dictionary full of zeros
    for i in range(291):
        all_fbd[i] = {}
        for x in range(1251):
            all_fbd[i][x] = 0

    for i in range(291):
        all_fbd[i][52+i],all_fbd[i][228+i],all_fbd[i][392+i],all_fbd[i][568+i],all_fbd[i][732+i],all_fbd[i][908+i] = P[0],P[1],P[2],P[3],P[4],P[5]
        # Find Fb (right side) and put it into the dictionary
        all_fbd[i][1250] = -(P[0]*(52+i) + P[1]*(228+i) + P[2]*(392+i) + P[3]*(568+i) + P[4]*(732+i) + P[5]*(908+i))/1250
        # Find Fa (left side) and put it into the dictionary
        all_fbd[i][0] = (-sum(P)) - all_fbd[i][1250]
    #print(all_fbd)

    all_sfd = {} # Dictionary like all_fbd
    max_shear = [] # The maximum shear possible with the train moving, we will graph this
    shear = 0 # A value to put into each placeholder
    
    # Make a dictionary full of zeros
    for i in range(291):
        all_sfd[i] = {}
        for x in range(1251):
            all_sfd[i][x] = 0
            
    for i in range(291):
        for x in range(1251):
            shear += all_fbd[i][x]
            all_sfd[i][x] = shear
        shear = 0
    #print(all_sfd)
    a = [] # a is a placeholder for all of the possible shears at each point of the bridge
    
    for x in range(1251):
        for i in range(291):
            a.append(all_sfd[i][x])
        if x <= 625:
            max_shear.append(max((a)))
        else:
            max_shear.append(min((a)))
        #print(a)
        a = []
        
    #print(max_shear)
        
    #example = list(all_sfd[290].values())
    #plt.plot(example)
    #plt.show()
    
    # Plot the max max_shear
    plt.plot(max_shear)
    plt.axhline(y = 0, color = 'r', linestyle = '-')
    plt.xlabel("Position from the Left (mm)")
    plt.ylabel("Max Shear Force (N)")
    plt.title("Shear - N Train")
    plt.grid()
    plt.show()
        
        
    all_bmd = {}
    max_moment = []
    
    for i in range(291):
        all_bmd[i] = {}
        for x in range(1251):
            all_bmd[i][x] = 0
    distance = 0
    for i in range(291):
        for x in range(1251):
            if x == 0:
                all_bmd[i][0] = 0
            else:
                all_bmd[i][x] = all_bmd[i][x-1] + all_sfd[i][x]
                
    b = [] # a is a placeholder for all of the possible moments at each point of the bridge
    
    for x in range(1251):
        for i in range(291):
            b.append(all_bmd[i][x])
        max_moment.append(max(b))
        b = []
                
    #print(all_bmd)  
    
    plt.plot(max_moment)
    plt.axhline(y = 0, color = 'r', linestyle = '-')
    plt.xlabel("Position from the Left (mm)")
    plt.ylabel("Max Moments (Nmm)")
    plt.title("MAX Moments from Train moving")
    plt.grid()
    plt.show() 
    
    V1 = max(np.absolute(max_shear))
    M1 = max(max_moment) 

def LoadCase2(F):  ## F is the weight of each wheel
    global V2
    global M2
    
    P = [0,0,0,0,0,0]
    P[0],P[1],P[2],P[3],P[4],P[5] = -(F*(1/7.17)),-(F*(1/7.17)), -(F*(1.1/7.17)),-(F*(1.1/7.17)), -(F*(1.485/7.17)),-(F*(1.485/7.17))
    print(sum(P))
    print(P)

    # The locomotive wagon (heaviest) is on the right side for these calculations

    all_fbd = {}
    
    # Make a dictionary full of zeros
    for i in range(291):
        all_fbd[i] = {}
        for x in range(1251):
            all_fbd[i][x] = 0

    for i in range(291):
        all_fbd[i][52+i],all_fbd[i][228+i],all_fbd[i][392+i],all_fbd[i][568+i],all_fbd[i][732+i],all_fbd[i][908+i] = P[0],P[1],P[2],P[3],P[4],P[5]
        # Find Fb (right side) and put it into the dictionary
        all_fbd[i][1250] = -(P[0]*(52+i) + P[1]*(228+i) + P[2]*(392+i) + P[3]*(568+i) + P[4]*(732+i) + P[5]*(908+i))/1250
        # Find Fa (left side) and put it into the dictionary
        all_fbd[i][0] = (-sum(P)) - all_fbd[i][1250]
    #print(all_fbd)

    all_sfd = {} # Dictionary like all_fbd
    max_shear = [] # The maximum shear possible with the train moving, we will graph this
    shear = 0 # A value to put into each placeholder
    
    # Make a dictionary full of zeros
    for i in range(291):
        all_sfd[i] = {}
        for x in range(1251):
            all_sfd[i][x] = 0
            
    for i in range(291):
        for x in range(1251):
            shear += all_fbd[i][x]
            all_sfd[i][x] = shear
        shear = 0
    #print(all_sfd)
    a = [] # a is a placeholder for all of the possible shears at each point of the bridge
    
    for x in range(1251):
        for i in range(291):
            a.append(all_sfd[i][x])
        if x <= 625:
            max_shear.append(max((a)))
        else:
            max_shear.append(min((a)))
        #print(a)
        a = []
        
    #print(max_shear)
        
    #example = list(all_sfd[290].values())
    #plt.plot(example)
    #plt.show()
    
    # Plot the max max_shear
    plt.plot(max_shear)
    plt.axhline(y = 0, color = 'r', linestyle = '-')
    plt.xlabel("Position from the Left (mm)")
    plt.ylabel("Max Shear Force (N)")
    plt.title("Shear - N Train")
    plt.grid()
    plt.show()
        
        
    all_bmd = {}
    max_moment = []
    
    for i in range(291):
        all_bmd[i] = {}
        for x in range(1251):
            all_bmd[i][x] = 0
    distance = 0
    for i in range(291):
        for x in range(1251):
            if x == 0:
                all_bmd[i][0] = 0
            else:
                all_bmd[i][x] = all_bmd[i][x-1] + all_sfd[i][x]
                
    b = [] # a is a placeholder for all of the possible moments at each point of the bridge
    
    for x in range(1251):
        for i in range(291):
            b.append(all_bmd[i][x])
        max_moment.append(max(b))
        b = []
                
    #print(all_bmd)  
    
    plt.plot(max_moment)
    plt.axhline(y = 0, color = 'r', linestyle = '-')
    plt.xlabel("Position from the Left (mm)")
    plt.ylabel("Max Moments (Nmm)")
    plt.title("MAX Moments from Train moving")
    plt.grid()
    plt.show() 
    
    V2 = max(max_shear)
    M2 = max(max_moment) 


def FOS0():
    # Tensile strength = 30MPa
    # Compressive Strength = 6MPa
    # Shear Strength (MBoard) = 4MPa
    # Young's Modulus = 4000MPa
    # Poisson Ratio = 0.2
    E = 4000
    
    # Contact Cement Shear Strength = 2Mpa
    
    #Set all the variables
    global V1
    global M1
    global V2
    global M2
 
    FOSten = 0
    FOScompr = 0
    # Case 1 is the buckling on the edge of the flange
    FOSflex1 = 0
    # Case 2 is the buckling on the middle of the flange
    FOSflex2 = 0
    # Case 3 is the buckling on the web
    FOSflex3 = 0
    FOSsh = 0
    FOSshglue = 0
    FOSshear_buck = 0
    
    # Get ybar
    th = 1.27  # Thickness
    A1 = th*(75-(2*th)) # Area of web without the tab
    A2 = th*100 # Area of top flange
    A3 = th*80 # Area of bottom flange
    AGlueTap = 6.27 * th # Area of glue tab
    height = 75 + th
    ybar = ((2*((75/2))*(A1))+((75+(th/2))*(A2))+((th/2)*A3)+(2*AGlueTap*(75-(th/2))))/((2*A1)+(A2)+(A3)+(2*AGlueTap))
    # Becaues the hand calculations roudn ybar, this will be rounded here too
    
    ybar = round(ybar,1)
    #Get I
    I1 = (th*((75-th-th)**3))/12
    I2 = ((th**3)*100)/12
    I3 = ((th**3)*80)/12
    Itab = ((th**3)*(5+th))/12
    yd1 = ybar - ((75/2))
    yd2 = (75+(th/2)) - ybar
    yd3 = (ybar - (th/2))
    ydtab = (75 - (th/2)) - ybar
    I = (2*I1) + I2 + I3 + (2*Itab) + (2*(A1*(yd1**2))) + (A2*(yd2**2)) + (A3*(yd3**2)) + (2*(AGlueTap*ydtab))
    
    # Get Q
    ## Maximum shear occurs at centroidal axis so
    AQshape = (75-th-ybar)*th
    yd = ((A2*(75+(th/2)-ybar)) + (2*AQshape*((75-th-ybar)/2)) + (2*AGlueTap*(75-(th/2)-ybar)))/ (A2 + (2*AQshape) + (2*AGlueTap))
    
    Qcent = yd * ((2*AQshape) + (2*AGlueTap) + A2)
    Qglue = A2 * (75+(th/2)-ybar)
    
    FOSten = 30 / ((M1*ybar)/I)
    FOScompr = 6 / ((M1*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 78.73
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/78.73)**2)) / ((M1*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 10.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/10.635)**2)) / ((M1*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (th+80+(th/2)-ybar) = 38.26
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(75-(th/2)-ybar))**2)) / ((M1*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V1*Qcent)/(I*th*2)
    tauglue = (V1*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 400 mm of distance between diaphragms is counted and b = 80mm
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/80)**2)+((th/400)**2))) / taucent
    print(M1, "Nmm")
    print(V1, "N")
    print(ybar, "mm")
    print((I/1000000), " *10^6 mm4")
    print(Qglue)
    print(Qcent)
    
    print("-----LOAD CASE 1------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M1*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M1*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    FOSten = 30 / ((M2*ybar)/I)
    FOScompr = 6 / ((M2*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 78.73
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/78.73)**2)) / ((M2*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 10.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/10.635)**2)) / ((M2*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (th+80+(th/2)-ybar) = 38.26
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(75-(th/2)-ybar))**2)) / ((M2*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V2*Qcent)/(I*th*2)
    tauglue = (V2*Qglue)/(I*((th*2)+10))
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 400 mm of distance between diaphragms is counted and b = 80mm
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/80)**2)+((th/400)**2))) / taucent
    print("-----LOAD CASE 2------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M2*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M2*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")

def FOS1():
    # Tensile strength = 30MPa
    # Compressive Strength = 6MPa
    # Shear Strength (MBoard) = 4MPa
    # Young's Modulus = 4000MPa
    # Poisson Ratio = 0.2
    E = 4000
    
    # Contact Cement Shear Strength = 2Mpa
    
    #Set all the variables
    global V1
    global M1
    global V2
    global M2
 
    FOSten = 0
    FOScompr = 0
    # Case 1 is the buckling on the edge of the flange
    FOSflex1 = 0
    # Case 2 is the buckling on the middle of the flange
    FOSflex2 = 0
    # Case 3 is the buckling on the web
    FOSflex3 = 0
    FOSsh = 0
    FOSshglue = 0
    FOSshear_buck = 0
    
    # Get ybar
    th = 1.27  # Thickness
    A1 = th*(60-(2*th)) # Area of web without the tab
    A2 = th*150 # Area of top flange
    A3 = th*110 # Area of bottom flange
    AGlueTap = 6.27 * th # Area of glue tab
    height = 60 + th
    ybar = ((2*((60/2))*(A1))+((60+(th/2))*(A2))+((th/2)*A3)+(2*AGlueTap*(60-(th/2))))/((2*A1)+(A2)+(A3)+(2*AGlueTap))
    # Becaues the hand calculations roudn ybar, this will be rounded here too
    
    ybar = round(ybar,1)
    #Get I
    I1 = (th*((60-th-th)**3))/12
    I2 = ((th**3)*100)/12
    I3 = ((th**3)*80)/12
    Itab = ((th**3)*(5+th))/12
    yd1 = ybar - ((60/2))
    yd2 = (60+(th/2)) - ybar
    yd3 = (ybar - (th/2))
    ydtab = (60 - (th/2)) - ybar
    I = (2*I1) + I2 + I3 + (2*Itab) + (2*(A1*(yd1**2))) + (A2*(yd2**2)) + (A3*(yd3**2)) + (2*(AGlueTap*ydtab))
    
    # Get Q
    ## Maximum shear occurs at centroidal axis so
    AQshape = (60-th-ybar)*th
    yd = ((A2*(60+(th/2)-ybar)) + (2*AQshape*((60-th-ybar)/2)) + (2*AGlueTap*(60-(th/2)-ybar)))/ (A2 + (2*AQshape) + (2*AGlueTap))
    
    Qcent = yd * ((2*AQshape) + (2*AGlueTap) + A2)
    Qglue = A2 * (60+(th/2)-ybar)
    
    FOSten = 30 / ((M1*ybar)/I)
    FOScompr = 6 / ((M1*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 110
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/110)**2)) / ((M1*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 20-0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/19.365)**2)) / ((M1*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (60 - (th/2) - ybar)
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(60-(th/2)-ybar))**2)) / ((M1*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V1*Qcent)/(I*th*2)
    tauglue = (V1*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 400 mm of distance between diaphragms is counted and b = 80mm
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(60-th))**2)+((th/200)**2))) / taucent
    print(M1, "Nmm")
    print(V1, "N")
    print(ybar, "mm")
    print((I/1000000), " *10^6 mm4")
    print(Qglue)
    print(Qcent)
    
    print("-----LOAD CASE 1------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M1*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M1*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    FOSten = 30 / ((M2*ybar)/I)
    FOScompr = 6 / ((M2*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 110
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/110)**2)) / ((M2*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 20-0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/19.365)**2)) / ((M2*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (60 - (th/2) - ybar)
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(60-(th/2)-ybar))**2)) / ((M2*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V2*Qcent)/(I*th*2)
    tauglue = (V2*Qglue)/(I*((th*2)+10))
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 400 mm of distance between diaphragms is counted and b = 60-th
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(60-th))**2)+((th/200)**2))) / taucent
    print("-----LOAD CASE 2------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M2*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M2*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    

def FOS2():
    # Tensile strength = 30MPa
    # Compressive Strength = 6MPa
    # Shear Strength (MBoard) = 4MPa
    # Young's Modulus = 4000MPa
    # Poisson Ratio = 0.2
    E = 4000
    
    # Contact Cement Shear Strength = 2Mpa
    
    #Set all the variables
    global V1
    global M1
    global V2
    global M2
 
    FOSten = 0
    FOScompr = 0
    # Case 1 is the buckling on the edge of the flange
    FOSflex1 = 0
    # Case 2 is the buckling on the middle of the flange
    FOSflex2 = 0
    # Case 3 is the buckling on the web
    FOSflex3 = 0
    FOSsh = 0
    FOSshglue = 0
    FOSshear_buck = 0
    
    # Get ybar
    th = 1.27  # Thickness
    A1 = th*(180-th-th) # Area of web without the tab
    A2 = th*100 # Area of top flange
    AGlueTap = 6.27 * th # Area of glue tab
    height = 180 + th
    ybar = ((2*((180/2))*(A1))+((180+(th/2))*(A2))+(2*AGlueTap*(180-(th/2))))/((2*A1)+(A2)+(2*AGlueTap))
    # Becaues the hand calculations roudn ybar, this will be rounded here too
    
    ybar = round(ybar,1)
    #Get I
    I1 = (th*((180-th-th)**3))/12
    I2 = ((th**3)*100)/12
    Itab = ((th**3)*(5+th))/12
    yd1 = ybar - ((180/2))
    yd2 = (180+(th/2)) - ybar
    ydtab = (180 - (th/2)) - ybar
    I = (2*I1) + I2 + (2*Itab) + (2*(A1*(yd1**2))) + (A2*(yd2**2)) + (2*(AGlueTap*ydtab))
    
    # Get Q
    ## Maximum shear occurs at centroidal axis so
    AQshape = (180-th-ybar)*th
    yd = ((A2*(180+(th/2)-ybar)) + (2*AQshape*((180-th-ybar)/2)) + (2*AGlueTap*(180-(th/2)-ybar)))/ (A2 + (2*AQshape) + (2*AGlueTap))
    
    Qcent = yd * ((2*AQshape) + (2*AGlueTap) + A2)
    Qglue = A2 * (180 + (th/2) - ybar)
    
    FOSten = 30 / ((M1*ybar)/I)
    FOScompr = 6 / ((M1*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 60-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(60-th))**2)) / ((M1*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(th/2))**2)) / ((M1*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (th+80+(th/2)-ybar) = 38.26
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(180-(th/2)-ybar))**2)) / ((M1*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V1*Qcent)/(I*th*2)
    tauglue = (V1*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (180-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(180-th))**2)+((th/200)**2))) / taucent
    print(M1, "Nmm")
    print(V1, "N")
    print(ybar, "mm")
    print((I/1000000), " *10^6 mm4")
    print(Qglue)
    print(Qcent)
    
    print("-----LOAD CASE 1------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M1*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M1*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    FOSten = 30 / ((M2*ybar)/I)
    FOScompr = 6 / ((M2*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 60-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(60-th))**2)) / ((M2*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(th/2))**2)) / ((M2*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab is (th+80+(th/2)-ybar) = 38.26
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(180-(th/2)-ybar))**2)) / ((M2*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V2*Qcent)/(I*th*2)
    tauglue = (V2*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (180-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(180-th))**2)+((th/200)**2))) / taucent
    print("-----LOAD CASE 2------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M2*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M2*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    
def FOS3():
     # Tensile strength = 30MPa
    # Compressive Strength = 6MPa
    # Shear Strength (MBoard) = 4MPa
    # Young's Modulus = 4000MPa
    # Poisson Ratio = 0.2
    E = 4000
    
    # Contact Cement Shear Strength = 2Mpa
    
    #Set all the variables
    global V1
    global M1
    global V2
    global M2
 
    FOSten = 0
    FOScompr = 0
    # Case 1 is the buckling on the edge of the flange
    FOSflex1 = 0
    # Case 2 is the buckling on the middle of the flange
    FOSflex2 = 0
    # Case 3 is the buckling on the web
    FOSflex3 = 0
    FOSsh = 0
    FOSshglue = 0
    FOSshear_buck = 0
    
    # Get ybar
    th = 1.27  # Thickness
    A1 = th*(160-th-th) # Area of web without the tab
    A2 = th*100*2 # Area of top flange
    AGlueTap = 6.27 * th # Area of glue tab
    height = 180 + th
    ybar = ((2*((160/2))*(A1))+((160+(th/2))*(A2))+(2*AGlueTap*(160-(th/2))))/((2*A1)+(A2)+(2*AGlueTap))
    # Becaues the hand calculations roudn ybar, this will be rounded here too
    
    ybar = round(ybar,1)
    #Get I
    I1 = (th*((160-th-th)**3))/12
    I2 = (((th*2)**3)*100)/12
    Itab = ((th**3)*(5+th))/12
    yd1 = ybar - ((160/2))
    yd2 = (160+(th/2)) - ybar
    ydtab = (160 - (th/2)) - ybar
    I = (2*I1) + I2 + (2*Itab) + (2*(A1*(yd1**2))) + (A2*(yd2**2)) + (2*(AGlueTap*ydtab))
    
    # Get Q
    ## Maximum shear occurs at centroidal axis so
    AQshape = (160-th-ybar)*th
    yd = ((A2*(160+(th/2)-ybar)) + (2*AQshape*((160-th-ybar)/2)) + (2*AGlueTap*(160-(th/2)-ybar)))/ (A2 + (2*AQshape) + (2*AGlueTap))
    
    Qcent = yd * ((2*AQshape) + (2*AGlueTap) + A2)
    Qglue = A2 * (160 + (th/2) - ybar)
    
    FOSten = 30 / ((M1*ybar)/I)
    FOScompr = 6 / ((M1*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 100-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(100-th))**2)) / ((M1*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(th/2))**2)) / ((M1*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab 
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(160-(th/2)-ybar))**2)) / ((M1*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V1*Qcent)/(I*th*2)
    tauglue = (V1*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (160-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(160-th))**2)+((th/200)**2))) / taucent
    print(M1, "Nmm")
    print(V1, "N")
    print(ybar, "mm")
    print((I/1000000), " *10^6 mm4")
    print(Qglue)
    print(Qcent)
    
    print("-----LOAD CASE 1------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M1*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M1*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    FOSten = 30 / ((M2*ybar)/I)
    FOScompr = 6 / ((M2*(height-ybar))/I)
    ## The width of the middle flange from center of web to centre is 100-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(100-th))**2)) / ((M2*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 0.635
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(th/2))**2)) / ((M2*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab 
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(160-(th/2)-ybar))**2)) / ((M2*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V2*Qcent)/(I*th*2)
    tauglue = (V2*Qglue)/(I*((th*2)+10))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (160-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(160-th))**2)+((th/200)**2))) / taucent
    print("-----LOAD CASE 2------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M2*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M2*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    

def FOS4():
     # Tensile strength = 30MPa
    # Compressive Strength = 6MPa
    # Shear Strength (MBoard) = 4MPa
    # Young's Modulus = 4000MPa
    # Poisson Ratio = 0.2
    E = 4000
    
    # Contact Cement Shear Strength = 2Mpa
    
    #Set all the variables
    global V1
    global M1
    global V2
    global M2
 
    FOSten = 0
    FOScompr = 0
    # Case 1 is the buckling on the edge of the flange
    FOSflex1 = 0
    # Case 2 is the buckling on the middle of the flange
    FOSflex2 = 0
    # Case 3 is the buckling on the web
    FOSflex3 = 0
    FOSsh = 0
    FOSshglue = 0
    FOSshear_buck = 0
    
    # Get ybar
    th = 1.27  # Thickness
    A1 = th*(160-th-th) # Area of web without the tab
    A2 = th*100*2 # Area of top flange
    AGlueTap = (11+th) * th # Area of glue tab
    height = 180 + th
    ybar = ((2*((160/2))*(A1))+((160+(th/2))*(A2))+(2*AGlueTap*(160-(th/2))))/((2*A1)+(A2)+(2*AGlueTap))
    # Becaues the hand calculations roudn ybar, this will be rounded here too
    
    ybar = round(ybar,1)
    #Get I
    I1 = (th*((160-th-th)**3))/12
    I2 = (((th*2)**3)*100)/12
    Itab = ((th**3)*(5+th))/12
    yd1 = ybar - ((160/2))
    yd2 = (160+(th/2)) - ybar
    ydtab = (160 - (th/2)) - ybar
    I = (2*I1) + I2 + (2*Itab) + (2*(A1*(yd1**2))) + (A2*(yd2**2)) + (2*(AGlueTap*ydtab))
    
    # Get Q
    ## Maximum shear occurs at centroidal axis so
    AQshape = (160-th-ybar)*th
    yd = ((A2*(160+(th/2)-ybar)) + (2*AQshape*((160-th-ybar)/2)) + (2*AGlueTap*(160-(th/2)-ybar)))/ (A2 + (2*AQshape) + (2*AGlueTap))
    
    Qcent = yd * ((2*AQshape) + (2*AGlueTap) + A2)
    Qglue = A2 * (160 + (th/2) - ybar)
    
    FOSten = 30 / ((M1*ybar)/I)
    FOScompr = 6 / ((M1*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 75-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(75-th))**2)) / ((M1*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 100-75-th
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/((100-75-th)/2))**2)) / ((M1*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab 
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(160-(th/2)-ybar))**2)) / ((M1*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V1*Qcent)/(I*th*2)
    tauglue = (V1*Qglue)/(I*((th*2)+22))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (160-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(160-th))**2)+((th/296)**2))) / taucent
    print(M1, "Nmm")
    print(V1, "N")
    print(ybar, "mm")
    print((I/1000000), " *10^6 mm4")
    print(Qglue)
    print(Qcent)
    
    print("-----LOAD CASE 1------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M1*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M1*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    FOSten = 30 / ((M2*ybar)/I)
    FOScompr = 6 / ((M2*(height-ybar))/I)
    # The width of the middle flange from center of web to centre is 75-th
    FOSflex1 = (((4*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(75-th))**2)) / ((M2*(height-ybar))/I)
    # The distance from the edge of the flange to the middle of the web is 100-75-th
    FOSflex2 = (((0.425*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/((100-75-th)/2))**2)) / ((M2*(height-ybar))/I)
    # From the centroidal axis to the middle of the glue tab 
    FOSflex3 = (((6*((np.pi)**2)*E)/(12*(1-(0.2**2))))*((th/(160-(th/2)-ybar))**2)) / ((M2*(height-ybar))/I)
    # Lets get the shear stresses first
    taucent = (V2*Qcent)/(I*th*2)
    tauglue = (V2*Qglue)/(I*((th*2)+22))
    
    FOSsh = 4 / taucent
    FOSshglue = 2 / tauglue
    
    # 200 mm of distance between diaphragms is counted and b = (160-th)
    FOSshear_buck = (((5*((np.pi)**2)*E)/(12*(1-(0.2**2)))) * (((th/(160-th))**2)+((th/296)**2))) / taucent
    print("-----LOAD CASE 2------")
    print("FOS tension: ", FOSten)
    print("FOS compression: ", FOScompr)
    print("FOS flex buckling 1: ", FOSflex1)
    print("FOS flex buckling 2: ", FOSflex2)
    print("FOS flex buckling 3: ", FOSflex3)
    print("FOS shear: ", FOSsh)
    print("FOS shear glue: ", FOSshglue)
    print("FOS shear buckling: ", FOSshear_buck)
    print("---------------------")
    print("ybar: ", ybar, " mm")
    print("I: ", I, " *10^6 mm4")
    print("First Moment of Inertia of Centroid: ", Qcent, " mm3")
    print("First Moment of Inertia at Glue: ",Qglue, " mm3")
    print("Compression Stress: ", ((M2*(height-ybar))/I), " MPa")
    print("Tension Stress: ", ((M2*ybar)/I), " MPa")
    print("Shear Stress at centroid: ", taucent, " MPa" )
    print("Shear Stress at glue: ", tauglue, " Mpa")
    
    
    

LoadCase1(400)
LoadCase2(478)
FOS0()
g   raphInertia()