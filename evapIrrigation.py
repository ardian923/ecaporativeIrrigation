#!/usr/bin/python

# Evaporative Irrigation modelling write by Dr.Ardiansyah for research purpose
# Free to adopt for further research. Contact : ardi.plj@gmail.com 

# we're going to use Microsoft Excel spreadsheet and xlrd xlwt module to link to Python
from xlrd import open_workbook
from xlwt import Workbook
from xlutils.copy import copy

# example writing and reading
################################
workbook_name = 'evapIrrigation.xls'
sheet_name    = '1-PlotEvaporativeIrrigation'
rb = open_workbook(workbook_name, formatting_info=True)
wb = copy(rb)

# for writing to excel
ws = wb.get_sheet(0)		# sheet indexing starts from 0

# for reading from excel
rs = rb.sheet_by_name(sheet_name)
################################

import numpy as np
from scipy.integrate import odeint
from scipy import interpolate as intr
import matplotlib as mpl
import matplotlib.pyplot as plt

i_row = 5

# initial conditions
hplot	= rs.cell_value(i_row,2) 		# water level in center of plot, measured from bottom root zone
hcon	= rs.cell_value(i_row+1,2) 		# water level in controller pipe, measured from bottom root zone
flood	= rs.cell_value(i_row+2,2) 		# 
SPin	= rs.cell_value(i_row+3,2) 		# 
SPout	= rs.cell_value(i_row+4,2) 		# 
Aplot	= rs.cell_value(i_row+5,2) 		# 
Acon	= rs.cell_value(i_row+6,2) 		# 
Re	= rs.cell_value(i_row+7,2) 		# 
ETo	= rs.cell_value(i_row+8,2) 		# 
Evap	= rs.cell_value(i_row+9,2) 		# 
Irr	= rs.cell_value(i_row+10,2) 		# 
Dr	= rs.cell_value(i_row+11,2) 		# 
Ks	= rs.cell_value(i_row+12,2) 		#
dx      = rs.cell_value(i_row+13,2) 		# 
Lp	= rs.cell_value(i_row+15,2) 		#
Wp      = rs.cell_value(i_row+16,2) 		# 
rc      = rs.cell_value(i_row+17,2) 		# 

rootZone  = 300 #mm, rootzone depth 30 cm
thetaInit = 0.4 # initial water content
watDepthInit = thetaInit * rootZone # mm, initial water depth, initial hplot
thetaC    = 0.4 # critical water content

def cropCoefficient(HST):
    #growth stage
    L_ini = 30
    L_dev = 30
    L_mid = 40
    L_late= 20
    #kc on growth stage
    kc_ini = 1.05
    kc_mid = 1.2
    kc_end = 0.8
    #calculate kc on growth stage
    if (HST <= L_ini):
        kc = kc_ini
    elif (HST <= (L_ini+L_dev)) and (HST > L_ini):
        kc = kc_ini + ((HST-L_ini)/(L_dev*1.0)) * (kc_mid - kc_ini) 
    elif (HST <= (L_ini+L_dev+L_mid)) and (HST > (L_ini+L_dev)):
        kc = kc_mid
    else: #(HST <= (L_ini+L_dev+L_mid+L_late)) and (HST > (L_ini+L_dev+L_mid)):
        kc = kc_mid - ((HST-(L_ini+L_dev+L_mid))/(L_late*1.0)) * (kc_mid - kc_end)
    return kc

def cropEvapoTranspiration(HST, ETo):
    kc = cropCoefficient(HST)
    ETc = kc * ETo
    return ETc

# solve the system dy/dt = f(y, t)
# dhplot/dt = (Irr + Re + SPin) - (SPout + Dr + ETc)
# dhcon/dt  = (Re - Evap + (Ks/dx) * (hplot - hcon))

# solve the system dy/dt = f(y, t)
def dhplot(t, hplot, hcon):
    #f = (1/Aplot)*(Irr + Re + SPin) - (SPout + Dr + ETc)
    f = (Irr + Re + SPin) - (SPout + Dr + ETc)
    return f

def dhcon(t, hplot, hcon):
    #f = (1/Acon) * (Re - Evap + ((2*Ks)/(rc*dx)) * hcon * (hplot - hcon))
    f = (Re - Evap - SPout + Lf)
    return f

def RK4(f1, f2):	# Rungge-Kuta Method
    k1 = dt * dhplot(t, hplot, hcon)
    k2 = dt * dhcon(t, hplot, hcon)
    l1 = dt * dhplot(t + dt/2., hplot + k1/2., hcon + k2/2.)
    l2 = dt * dhcon(t + dt/2., hplot + k1/2., hcon + k2/2.)
    m1 = dt * dhplot(t + dt/2., hplot + l1/2., hcon + l2/2.)
    m2 = dt * dhcon(t + dt/2., hplot + l1/2., hcon + l2/2.)
    n1 = dt * dhplot(t + dt   , hplot + m1   , hcon + m2   )
    n2 = dt * dhcon(t + dt   , hplot + m1   , hcon + m2   )
    df1 = (k1 + 2*l1 + 2*m1 + n1)/6.
    df2 = (k2 + 2*l2 + 2*m2 + n2)/6.
    return df1, df2

if __name__ == '__main__':  ## Run simulation as standalone program
    dt = 1
    HST_end = 121
    HST =  np.arange(1., HST_end, dt) #create time array from day 1 to day 44

    #allocate space for array variables
    arhplot =  np.zeros(HST.size)
    arhcon  =  np.zeros(HST.size)
    arKc    =  np.zeros(HST.size)
    arETc   =  np.zeros(HST.size)
    arcumETc=  np.zeros(HST.size)
    arIrr   =  np.zeros(HST.size)
    arRe    =  np.zeros(HST.size)
    arDr    =  np.zeros(HST.size)
    arSPin  =  np.zeros(HST.size)
    arSPout =  np.zeros(HST.size)
    cumETc  = 0
    Lf      = 0


    for j in range(HST.size):
        t = HST[j]
        Kc  = cropCoefficient(t)
        ETc = cropEvapoTranspiration(t, ETo)
        cumETc = cumETc + ETc

#        if (hplot <= (thetaC*rootZone)): #if water depth in plot less than critical water content
        if (hplot <= (thetaC*rootZone)): #if water depth in plot less than critical water content
            Irr = ETc                    #give irrigation as much as ETc
        else:
            Irr = 0
        df1, df2 = RK4(lambda t, hplot, hcon: f1(t, hplot, hcon), lambda t, hplot, hcon: f2(t, hplot, hcon))
        t, hplot, hcon = t + dt, hplot + df1, hcon + df2
        #if hcon less than zero, set to zero
        if (hcon <= 0):
            hcon = 0
        Lf = ((2*Ks)/(rc*dx)) * hcon * (hplot - hcon)
        print (t,',', hplot, hcon)
        arhplot[j] = hplot
        arhcon[j]  = hcon
        arKc[j]    = Kc
        arETc [j]  = ETc
        arcumETc [j]  = cumETc
        arIrr[j]   = Irr
        arRe[j]    = Re 
        arDr[j]	   = Dr
        arSPin[j]     = SPin
        arSPout[j]    = SPout 
        #print to excel spreadsheet
        ws.write(j+6, 9, HST[j]); ws.write(j+6, 10, arhplot[j]); ws.write(j+6, 11, arhcon[j]);
        ws.write(j+6, 12, arKc[j]); ws.write(j+6, 13,arETc[j]); ws.write(j+6, 14, arcumETc[j]); ws.write(j+6, 15, arIrr[j]);
        ws.write(j+6, 16, arDr[j]); ws.write(j+6, 17,arRe[j]); ws.write(j+6, 18, arSPin[j]);
        ws.write(j+6, 19, arSPout[j])
        ws.write(j+6, 20, Lf)
    # save to excel files
    wb.save(workbook_name)

    # plot graph
    plt.rcParams['figure.figsize'] = 8, 7
    mpl.rcParams.update({'font.size': 20})

    # plot HST - hplot, hcon
    fig1, ax1 = plt.subplots()
    ax1.plot(HST, arhplot, 'r-', label='hplot')
    ax1.plot(HST, arhcon, 'k-', label='hcon')
    ax1.legend(loc='upper right', shadow=False, frameon=False)
    ax1.set_xlabel('HST ($hari$)')
    ax1.set_ylabel('Level Air ($mm$)')

    # Plot Water balance in paddy plot
    fig2, ax2 = plt.subplots()
    ax2.plot(HST, arETc, 'g--', label='ETc')
    ax2.plot(HST, arIrr, 'b--', label='Irr')
    ax2.plot(HST, arDr, 'k-', label='Dr')
    ax2.plot(HST, arSPin, 'r--', label='SPin')
    ax2.plot(HST, arSPout, 'r-.', label='SPout')
    ax2.plot(HST, arRe, 'b--', label='Re')
    ax2.legend(loc='upper left', shadow=False, frameon=False, ncol=2)
    ax2.set_ylim([0, 8])
    ax2.set_xlabel('HST ($hari$)')
    ax2.set_ylabel('Ketebalan Air (mm)', color='k')
    plt.show()

