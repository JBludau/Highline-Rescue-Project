#! python

from sympy import symbols, Eq, solve, re, im
from scipy import optimize as optimize
import math
import matplotlib.pyplot as plt
import numpy as np
import csv

def calculateForcesAthleteCentered(Weight,Length,Modulus):

    x = symbols('x')

    eq1 = Eq(lhs=x, rhs=(Weight*Weight * (1+x/Modulus)**2 / ((1+x/Modulus)**2 -1))**0.5)

    sol = solve(eq1)
    Force = re(sol[0])
    print(f"Force is {Force} N")
    alpha = math.asin(Weight/2/Force)
    print(f"Resulting angle {alpha/math.pi*180} degrees")
    TotLength = Length*(1+Force/Modulus)
    print(f"Totallength of line {TotLength} m")
    Sag = math.sin(alpha)*0.5*TotLength
    print(f"Sag in middle {Sag} m")
    print("finished")

    return Force,alpha,TotLength,Sag

#Creates vector of equations that are used to calculate the forces for a free highline with an athlete in a non-centered position  as described in highline_hoist_contact.pdf
def athleteAtNonCenterPosition(x,l_0, l_L, l_R, F, E, PreTension):
    fun = np.ones(7,float);
    fun[0] = (math.tan(x[0])-x[2]/l_L)
    fun[1] = (math.tan(x[1])-x[2]/l_R)
    fun[2] = (math.cos(x[0])-l_L/x[3])
    fun[3] = (math.cos(x[1])-l_R/x[4])
    fun[4] = (x[3]+x[4]-x[5])
    fun[5] = ((math.sin((x[0]))+math.sin((x[1]))) * (x[6]+PreTension) - F)
    fun[6] = (x[5]-l_0*(1+x[6]/E))

    return fun


#Creates vector of equations that are used to calculate the forces for a highline in contact with a hoist cable  with an athlete in a non-centered position  as described in highline_hoist_contact.pdf
def athleteAtNonCenterPositionContactHeli(x,l_0, l_L, l_R, F_A, E,PreTension, F_R, h_y, h_z):
    fun = np.ones(11,float);
    fun[0] = math.tan(x[0])-x[2]/l_L
    fun[1] = (math.tan(x[1])-x[2]/l_R)
    fun[2] =(math.cos(x[0])-l_L/x[3])
    fun[3] =(math.cos(x[1])-l_R/x[4])
    fun[4] =(x[3]+x[4]-x[5])
    fun[5] =((math.sin((x[0]))+math.sin((x[1]))) * (x[6]+PreTension) - (F_A/math.cos(x[9])))
    fun[6] =(x[5]-l_0*(1+((x[6])**2+x[7]**2)**0.5/E))

    fun[7]= (math.tan(x[9]) - x[7]/F_A)
    fun[8]=(math.sin(x[10])*x[8]-x[7])
    fun[9]=(math.cos(x[10])*x[8] - F_R)

    fun[10]=(math.tan(x[10])-(h_y-math.sin(x[9])*x[2])/(h_z-x[2]*(1-math.cos(x[9]))))

    return fun


if __name__ == "__main__":

    #mechanical parameters
    Lengths_pink = [30,40,60,80,100,120,160,200]
    Modulus_pink = 30000
    Lengths_y2k = [60,80,100,120,160,200,300,400]
    Modulus_y2k = 420000

    #athlete related parameters (out of center)
    positions = [0.2,0.3,0.4,0.5]
    weights = [60,80,100,120]

    #heli related parameters
    heli_y=[1,2,3,4,5,6,7,8,9,10]
    heli_z=[3,5,7,9,11,13,15,17,19,21,24,27,30]
    F_R = 100*9.81

    #highline related parameters
    Modulus = Modulus_y2k
    Lengths = Lengths_y2k
    PreTensions=[500,1000,1500,2000,2500,3000,3500]

    #output config
    base_name = 'results_y2k'
    results = []

    for weight in weights:
        for length in Lengths:
            for PreTension in PreTensions:
                for position in positions:
                            #if not converged the result should not be tabulated ... try/except is not the best idea for this but it works for now (e.g. division by 0 might appear in the solvers internals)
                            try:
                                l_0 = length
                                l_L = position*l_0
                                l_R = 1.0*(1-position)*l_0
                                F_A = weight*9.81
                                E = Modulus
                    
                                #start values for iterative solver
                                guess = np.ones(7, float)
                                guess[0] = 4/180*math.pi
                                guess[1] = 4/180*math.pi
                                guess[2] = 2
                                guess[3] = l_L
                                guess[4] = l_R
                                guess[5] = l_0
                                guess[6] = 2000

                                # solution step ... the Levenberg Marquardt proved to be the better solver
                                #  sol = optimize.root(athleteAtNonCenterPositionContactHeli, guess,args=(l_0,l_L,l_R,F_A,E), method='krylov',tol=0.000001, options = {'maxiter':10000000})
                                sol = optimize.root(athleteAtNonCenterPosition, guess,args=(l_0,l_L,l_R,F_A,E,PreTension), method='lm')

                                results.append([sol.success,weight,length,position,l_0,l_R,l_L,F_A,E,PreTension,sol.x[0]/math.pi*180,sol.x[1]/math.pi*180,sol.x[2],sol.x[3],sol.x[4],sol.x[5],sol.x[6]+PreTension])
                                print(f"success with params {weight},{length},{position},{l_0},{l_R},{l_L},{F_A},{E}")

                            except:
                                print(f"failed with params {weight},{length},{position},{l_0},{l_R},{l_L},{F_A},{E}")

    #print finish notification to user and write results to file
    fields=['#success','weight','length','position','l_0','l_R','l_L','F_A','E','PreTension','alpha_L','alpha_R','d_0','S_L','S_R','S','TensionForce']
    with open(base_name+'.csv', 'w') as f:
        write = csv.writer(f)
        write.writerow(fields)
        for item in results:
            write.writerow(item)


    # do same calculations for contact with helicopter
    count_success = 0
    count_fail = 0
    results = []
    for weight in weights:
        for length in Lengths:
            for PreTension in PreTensions:
                for position in positions:
                    for h_z in heli_z:
                        for h_y in heli_y:

                            try:
                                l_0 = length
                                l_L = position*l_0
                                l_R = 1.0*(1-position)*l_0
                                F_A = weight*9.81
                                E = Modulus

                                guess = np.ones(11, float)
                                guess[0] = 4/180*math.pi
                                guess[1] = 4/180*math.pi
                                guess[2] = 2
                                guess[3] = l_L
                                guess[4] = l_R
                                guess[5] = l_0
                                guess[6] = 2000
                                guess[7] = 100
                                guess[8] = 800
                                guess[9] = 3/180*math.pi
                                guess[10] = 5/180*math.pi

                                # solution step ... the Levenberg Marquardt proved to be the better solver
                                #  sol = optimize.root(athleteAtNonCenterPositionContactHeli, guess,args=(l_0,l_L,l_R,F_A,E,F_R,h_y,h_z), method='krylov',tol=0.000001, options = {'maxiter':10000000})
                                sol = optimize.root(athleteAtNonCenterPositionContactHeli, guess,args=(l_0,l_L,l_R,F_A,E,PreTension,F_R,h_y,h_z), method='lm')

                                results.append([sol.success,weight,length,position,l_0,l_R,l_L,F_A,E,PreTension,F_R,h_y,h_z,sol.x[0]/math.pi*180,sol.x[1]/math.pi*180,sol.x[2],sol.x[3],sol.x[4],sol.x[5],sol.x[6]+PreTension,sol.x[7],sol.x[8],sol.x[9]/math.pi*180,sol.x[10]/math.pi*180])
                                if(sol.success):
                                    count_success = count_success +1
                                else:
                                    count_fail = count_fail + 1
                                print(f"success with params {weight},{length},{position},{l_0},{l_R},{l_L},{F_A},{E},{PreTension},{F_R},{h_y},{h_z}")
                            except:
                                print(f"failed with params {weight},{length},{position},{l_0},{l_R},{l_L},{F_A},{E},{F_R},{h_y},{h_z}")
                                


    #inform user and print to file
    print(f"ratio of success/fail = {count_success} / {count_fail}")

    fields=['#success','weight','length','position','l_0','l_R','l_L','F_A','E','PreTension','F_R','h_y','h_z','alpha_L','alpha_R','d_0','S_L','S_R','S','TensionForce','F_C','F_H','Psi','beta']
    with open(base_name + '_heli.csv', 'w') as f:
        write = csv.writer(f)
        write.writerow(fields)
        for item in results:
            write.writerow(item)




#  plt.show()
