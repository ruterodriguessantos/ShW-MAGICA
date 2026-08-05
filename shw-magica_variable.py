#!/usr/bin/env python
# coding: utf-8
'''
**************************************************************************
ShW-MAGICA

--> Spatially varying electric field implementation without explicit shield wire modelling <--

Python implementation of the Lehtinen–Pirjola modified method for computing
Geomagnetically Induced Currents (GICs) in power transmission networks,
including the explicit modelling of shield wires.

This code is based on the GEOMAGICA framework, originally developed from:

    • P. Weidelt – thin-sheet approximation algorithm
    • C. Beggan, K. Turnbull and A. McKay (British Geological Survey) –
      adaptation and original Python implementation
    • R. Bailey – adaptation for the Austrian power network
https://github.com/geomagpy/GEOMAGICA

The present version extends and adapts the original implementation for
the Portuguese transmission network (ShW-MAGICA). The main modifications
introduced are:

    • Electric field interpolation using Delaunay triangulation.
    • Fixed 350 m integration step for line integral calculations.
    • Implementation of the LPm method to avoid numerical infinite
      grounding resistances.
    • Extension of the network and circuit model to explicitly represent
      shield wires.
    • Inclusion of shield wire equivalent circuit parameters for GIC
      computation.

Contributors to this version:
    • Rute Rodrigues dos Santos
    • Alexandra Pais
    • Joana Alves Ribeiro
    • Fernando Pinheiro

Developed at the University of Coimbra
CITEUC 


**************************************************************************
'''
import os
import sys
import getopt
import numpy as np
import pandas as pd
from scipy import interpolate
from math import radians, tan, atan, atan2, cos, sin, acos, asin
from math import sqrt, pi, log
import pickle
import numpy, scipy.io
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import scipy
import timeit
import tkinter as tk
from tkinter.filedialog import askopenfilename
try:
    import IPython
except:
    pass

def input_output(X):

    testfield       = False  
    
    networkpath = 'Data/GRID.txt'
    connectionspath =  'Data/CONNECTIONS.txt'

    temp=Efield_pkl.split('/')
    filename1=temp[-1]
    SavePath=Efield_pkl[0:-len(filename1)]
         
    out_name        =   'variable_Efield_'+filename1[:-4]                  
    
    return testfield,networkpath,connectionspath,Efield_pkl,Points_pkl,out_name, SavePath;


#####################################################################
#                       FUNCTIONS                                   #
#####################################################################

def grc_azimuth(lonlat1, lonlat2):
    """
    Function to compute the geographic distance on a prolate ellipsiod such
    as a planet or moon. This computes the distance accurately - not that it
    makes much difference in most cases, as the location of the points of interest is
    generally relatively poorly known.
    
    This function is based on the formula from the Wikipedia page, verified
    against the Geoscience Australia website calculator
    
    Author: Ciaran Beggan
    Rewritten from matlab into python by R Bailey, ZAMG.
    
    Returns azimuth between two points.
    """

    a, b = 6378.137, 6356.752
    f = (a-b)/a

    u1 = atan((1.-f)*tan(pi/180.*(lonlat1[1])))
    u2 = atan((1.-f)*tan(pi/180.*(lonlat2[1])))

    L = pi/180.*(lonlat2[0] - lonlat1[0])
    Lambda, converge, iters = L, False, 0

    while not converge and iters < 20:
        sinsig = sqrt((cos(u2)*sin(Lambda))**2. + (cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(Lambda))**2.)
        cossig = sin(u1)*sin(u2) + cos(u1)*cos(u2)*cos(Lambda)
        sig = atan2(sinsig, cossig)
    
        sinalpha = (cos(u1)*cos(u2)*sin(Lambda))/sinsig
        cossqalpha = 1. - sinalpha**2.
        cos2sigm = cossig - (2.*sin(u1)*sin(u2))/cossqalpha

        C = (f/16.) * cossqalpha*(4. + f*(4.-3.*cossqalpha))
    
        calclambda = L + (1.-C)*f*sinalpha*(sig + C*sinalpha*(cos2sigm + C*cossig*(-1. + 2.*cos2sigm) ))
    
        if (abs(Lambda - calclambda) < 10.**(-12.)):
            converge = True
            Lambda = calclambda
        else:
            iters = iters + 1
            Lambda = calclambda

    usq = cossqalpha * ((a**2. - b**2.)/b**2.)
    A = 1. + usq/16384. * (4096. + usq*(-768. + usq*(320. - 175.*usq)))
    B = usq/1024.* (256. + usq*(-128. + usq*(74. - 47.*usq)))
    delsig = B * sinsig * (cos2sigm + 0.25 * B *(cossig *(-1. + 2.*cos2sigm) -(1./6.)*B * cos2sigm*(-3. + 4.*sinalpha**2.)*(-3.+4.*cos2sigm**2.)   ))
    s = b*A*(sig - delsig)
    a1 = atan2(cos(u2)*sin(Lambda), cos(u1)*sin(u2) - sin(u1)*cos(u2)*cos(Lambda) )
    #a2 = atan2(cos(u1)*sin(Lambda), -sin(u1)*cos(u2) + cos(u1)*sin(u2)*cos(Lambda) )    
    #print usq, A, B, delsig, s, a1
 
    if np.isnan(a1):
        a1 = 0.

    #a1 = -a1
    if a1 < 0.:
        a1 = 2.*pi + a1

    return a1


def grc_distance(lat1, lon1, lat2, lon2, result='km'):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees) using the Haversine method.
    Combination of:
    http://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points
    &
    http://gis.stackexchange.com/questions/29239/calculate-bearing-between-two-decimal-gps-coordinates
    """

    # Convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # Haversine formula for distance:
    # Source: Wiki https://en.wikipedia.org/wiki/Haversine_formula
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    #c =  2.* asin( sqrt( sin((dlat)/2.)**2. + cos(lat1)*cos(lat2)* sin((dlon)/2.)**2. ) )
    r = 6371. # Radius of earth in kilometers. Use 3956 for miles
    
    if dlat == 0. and dlon == 0.:
        return 0.

    # Great circle distance:
    c = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(dlon))

    if result == 'km':
        return c * r
    elif result == 'rad':
        return c


#####################################################################
#                       MAIN PROGRAM                                #
#####################################################################

# ===============================================================
# 0) READ IN OPTIONS
# ===============================================================
    
usage = ("-------------------------------------------------------------------",
             "DESCRIPTION:",
             "  Python script for modelling GICs using the Lehtinen and Pirjola",
             "  (1985) method of GIC computation applied to Horton et al. (2012)",
             "  example network model. The input geoelectric field can be a basic",
             "  1 V/km field or a file for a spatially varying field (for",
             "  Austria) can be defined instead. After computation, the values",
             "  of GIC per network node and per transformer are printed.",
             "  Details on the Horton grid are provided in the folder 'network'.",
             "-------------------------------------------------------------------",
             "OPTIONS:",
             "  -e/--efile:       Defines input geoelectric field file path that",
             "                    replaces 1 V/km electric field option.",
             "                    python GIC_Model_Horton.py -e <efilepath>",
             "  -h/--help:        Prints this helpful text.",
             "                    python GIC_Model_Horton.py -h",
             "-------------------------------------------------------------------",
             "EXAMPLE USAGE:",
             "  - Basic model with 1 V/km geoelectric field values:",
             "    $ python GIC_Model_Horton.py",
             "  - With defined geoelectric field input file from thin-sheet code:",
             "    $ python GIC_Model_Horton.py -e Efiles/E_39_2017-09-07T23:25:00.txt",
             "-------------------------------------------------------------------",
              )
    

#!/usr/bin/env python
'''
********************************************************************
Basic GIC calculation method following Lehtinen and Pirjola 1985.
Adapted from ComputeGIC.m from C. Beggan's (BGS) adaptation of
K. Turnbull's Fortran code.

For usage to calculate expected GIC values in the grid from the 
Horton et al. (2012) paper with a 1 V/km geoelectric field execute:
    $ python GIC_Model_Horton.py

NOTES:
-   Due to differences in distance calculations and Python Numpy
    methods, final values are off the exact Horton values by
    a few percent in most cases.
-   Lehtinen-Pirjola method needs more null points added in to 
    equate to nodal admission matrix method, hence the representation
    of lines between HV, LV, ground and switch nodes.
-   Note that in this case the transformer resistance is explicitly 
    given as a connection between nodes.

Created by R Bailey (ZAMG, Austria) on 2015-08-03.
#********************************************************************
'''

testfield,networkpath,connectionspath,Efield_pkl,Points_pkl,out_name,SavePath      =input_output('x')

# ===============11================================================
# 1) DEFINE NETWORK CONFIGURATION
# =================================================================

# Read station and transformer data:
network = open(networkpath, 'r')
netdata = network.readlines()

# Define number of nodes:
nnodes = len(netdata)

# Read data into arrays:
npf = np.float32
geolat, geolon = np.zeros(nnodes, dtype=npf), np.zeros(nnodes, dtype=npf)
country = []
sitecode, sitename = [], []
sitenum = np.zeros(nnodes, dtype=np.int32)
res_earth = np.zeros(nnodes, dtype=npf)
res_trans = np.zeros(nnodes, dtype=npf)
# Station-to-index and index-to-station dicts:
s2i, i2s = {}, {}       

for i in range(nnodes):
    # READ NODE DATA:
    # ---------------
    data = netdata[i].split("\t")
    #print(data[0],data[2])
    # Dictionary for station to index:
    # In this way, the index for stations starts in '0' and ends in 'nnodes-1'
    s2i[data[2]] = int(data[0]) - 1
    # Dictionary for index to station:
    i2s[int(data[0]) - 1] = data[2]
    # Site number:
    sitenum[i] = int(data[0]) - 1           # -1 to simplify python indices
    # Site readable names:
    sitename.append(data[1])
    # Site code names:
    sitecode.append(data[2])
    # Site country:
    country.append(data[3])
    # Geographic latitude of node:
    geolat[i] = float(data[4])
    # Geographic longitude of node:
    geolon[i] = float(data[5])
    # Earthing resistance of each node:
    res_earth[i] = float(data[6])
    # Transformer resistances (see notes):
    res_trans[i] = float(data[7])
    
# READ CONNECTION DATA:
# ---------------------
connectionsfile = open(connectionspath, 'r')
conndata = connectionsfile.readlines()

# Number of connections:
nconnects = len(conndata)

#alt jaribeiro --> RUTE ========================================
substation, TransType, Wind, TrueLength, TypeShieldWire, Phase_Exception = [], [], [], [], [], []
nodefrom = np.zeros(nconnects, dtype=np.int32)
nodeto = np.zeros(nconnects, dtype=np.int32)
res_line = np.zeros(nconnects, dtype=npf)
voltage_lines = np.zeros(nconnects, dtype=np.int32)

for i in range(0, nconnects):
    conns = conndata[i].split("\t")
    ##      # Connection starts at site:
    ##      nodefrom[i] = int(s2i[conns[1]])
    ##       # Connection ends at site:
    ##       nodeto[i] = int(s2i[conns[2]])
    ##       # Line resistance between connecting sites:
    ##       # (Divide by 3 for full transformer representation as paper values are given per phase)
    ##       res_line[i] = float(conns[5])/3.
    ##       # Voltage level of line:
    ##       voltage_lines[i] = float(conns[6]) 
    ####--------------------------------------------------------------------------------------
    #alt jaribeiro
    # Substation name:
    substation.append (conns[0])
    # Transformer type:
    TransType.append(conns[1])
    # Transformer winding:
    Wind.append(conns[2])
    # Connection starts at site:
    nodefrom[i] = int(s2i[conns[4]])
    # Connection ends at site:
    nodeto[i] = int(s2i[conns[5]])
    # True line length
    TrueLength.append(conns[6])
    # Type of shield wire cable(
    TypeShieldWire.append(conns[7])
    # Two phase exception == '1'
    Phase_Exception.append (conns[11])
    # Line resistance between connecting sites:
    # jaribeiro ALT            
    # (Divide by 3 for full transformer representation as paper values are given per phase, and divide for 2 for connections lines in the exception file)
    if Phase_Exception[i] == 'X':
        res_line[i] = float(conns[12])/2.
    else:
        res_line[i] = float(conns[12])/3. 
     
    # Voltage level of line:
    voltage_lines[i] = float(conns[13])
####---------------------------------------------------------------------------------------

# Set inf values to very high resistance for computation:
valInf = 1.e8
res_trans[res_trans==np.inf] = valInf
res_earth[res_earth==np.inf] = 0
res_line[res_line==np.inf] = 0

# ===============================================================
# 2) ESTABLISH MATRICES FOR LP1985 METHOD
# ===============================================================

# Define resistance matrices:
resis = np.zeros((nnodes, nnodes), dtype=npf)
connections = np.zeros((nnodes,nnodes), dtype=npf)

# I changed a bit the original if statement, in order to guarantee that matrix resis
# is filled symmetrically
#               A.P.

for i in range(0,nconnects):
    x, y = int(nodefrom[i]), int(nodeto[i])
    if (resis[x,y] > 0. and res_line[i] > 0.):
        # If res already added to this line, add another line in parallel:
        resis[x,y] = 1./(1./resis[x,y] + 1./res_line[i])
        resis[y,x] = resis[x,y]
    else:
        resis[x,y] = res_line[i]
        resis[y,x] = resis[x,y]
    if resis[x,y] > 0.:
        connections[x,y] = 1./resis[x,y]
        connections[y,x] = 1./resis[x,y]

# Calculate matrix of distance between each point:
dists = np.zeros((nnodes, nnodes), dtype=npf)
azi = np.zeros((nnodes, nnodes), dtype=npf)

# I changed a bit the original for statement, in order to guarantee that dist(j,i)
# and azi(j,i) are correctly inserted. Before, azi(j,i) was made equal to azi(i,j)!...
#               A.P.

for i in range(0,nnodes):
    for j in range(0,nnodes):
        lati, loni = geolat[i], geolon[i]
        latj, lonj = geolat[j], geolon[j]
        if resis[i, j] != 0.:
            dists[i, j] = grc_distance(lati, loni, latj, lonj, result='km')
            dists[j, i] = dists[i, j]
            #dists[j, i] = grc_distance(latj, lonj, lati, loni, result='km')
            try:
                azi[i, j] = grc_azimuth([loni, lati], [lonj, latj])
                if azi[i, j] >= 0. and azi[i, j] <= np.pi:
                    azi[j, i] = azi[i, j] + np.pi
                if azi[i, j] > np.pi and azi[i, j] <= 2*np.pi:
                    azi[j, i] = azi[i, j] - np.pi
                #azi[j, i] = grc_azimuth([loni, lati], [lonj, latj])
            except ZeroDivisionError:
                azi[i, j] = np.nan
                azi[j, i] = np.nan
                
                
# Create earth impedance matrix:
res_earth[res_trans != 0]=0  # set zero to substation that don't have transformers connected to ground
adm_earth=res_earth.copy()
adm_earth[adm_earth != 0] = 1./(adm_earth[adm_earth != 0]) 
earthadmit = np.diag(adm_earth)


# Calculate network admittance matrix:
# LP1984 eq. (10): **Y**
netadmit = -1.*connections + np.diag(sum(connections))

systemmat = netadmit + earthadmit
print(np.linalg.det(systemmat))
# ===============================================================
# 3) DEFINE LOCATION VARIABLES
# ===============================================================

# Boundaries of Portuguese territory where the grid is located:
nbound_NA = 39.0
ebound_NA = -7.0
sbound_NA = 37.0
wbound_NA = -9.0


    
# Geoelectric field provided, computed elsewhere using some conductivity model
Efield_in=open(Efield_pkl,'rb')      
Points_in=open(Points_pkl,'rb')

#
# Reading the Einduced field pickle file
# Some parameters are determined:
# nfr --> number of time series points
# npt --> number of Earth surface points where induced field is given
#
Edata=pickle.load(Efield_in)

nfr = Edata[0][0].size
aux = Edata[0].size
npt = aux // nfr

Ex  = np.zeros(  (npt,nfr)   ,dtype=float )
Ey  = np.zeros(  (npt,nfr)   ,dtype=float )
Ex=Edata[0]
Ey=Edata[1]

mLt  = []
mLn  = []
mLt[0:npt],mLn[0:npt]=pickle.load(Points_in)

S_lt=mLt[102]
W_lt=mLt[26]  
W_ln=mLn[62]
#
# generate  X(Npoints, Ndims) ndarray of floats -->
#           --> Data point coordinates --> NDinterpolator input
#      
pts=[]
for j in range(0,len(mLt)):
      pts.append(   [ mLn[j] , mLt[j] ]   )   
          

#
# alt Rute      
#
# Before cycle over time, create all path for lines
#
# Now, we don't have the same number of steps per line (200)
#                              steps = dist_line * 3
# pathlatsteps and pathlonsteps are matrix of all steps size(nconnects,max_steps)
# nconnects --> total number of lines
# max_steps --> dist(longest line)*3
#

max_steps = 320 
pathlatsteps, pathlonsteps = np.zeros((nconnects,max_steps), dtype=npf), np.zeros((nconnects,max_steps), dtype=npf)  

for i in range(0,nconnects):                            

    slat = geolat[nodefrom[i]]
    slon = geolon[nodefrom[i]]

    flat = geolat[nodeto[i]]
    flon = geolon[nodeto[i]]

    steps = dists[nodefrom[i],nodeto[i]] * 3   
    print(steps)

    if steps != 0. :
        
        # in case of very short lines
        if steps <10: # menos de 3km colocar step fixo !!!!!!!!!!!
            steps=10

        vector = np.linspace(slat, flat, int(steps))
        N = max_steps-vector.size
        new_vector = np.pad(vector, (0, N), 'constant')
        pathlatsteps[i,:] = new_vector

        vector = np.linspace(slon, flon, int(steps))
        N = max_steps-vector.size
        new_vector = np.pad(vector, (0, N), 'constant')
        pathlonsteps[i,:] = new_vector



E_n_sub, E_e_sub = np.zeros((nnodes,nfr), dtype=npf), np.zeros((nnodes,nfr), dtype=npf)      

#
# The cycle over time points starts here.
# While running, the GIC times values are appended to list 'global_results'
#
    
global_results=[]
current_sources_results=[]
results_fem=[]
    
for k in range(0,nfr):

    valuesEx = Ex[0:npt,k]
    valuesEy = Ey[0:npt,k]

    en_int_N = scipy.interpolate.NearestNDInterpolator(pts,valuesEx)
    ee_int_N = scipy.interpolate.NearestNDInterpolator(pts,valuesEy)

    en_int_D = scipy.interpolate.LinearNDInterpolator(pts,valuesEx)
    ee_int_D = scipy.interpolate.LinearNDInterpolator(pts,valuesEy)

    # ===============================================================
    # 5) INTEGRATE FIELD ALONG LINES
    # ===============================================================
    # Number of steps in path integration (more is slower but more exact):

    Vn_tot, Ve_tot = np.zeros(nconnects, dtype=npf), np.zeros(nconnects, dtype=npf)
    E_e, E_n = np.zeros((nconnects,max_steps), dtype=npf), np.zeros((nconnects,max_steps), dtype=npf)            

    for i in range(0,nconnects):

            number_steps_lon=np.count_nonzero(pathlonsteps[i,:])
            for p in range(0, number_steps_lon-1): #len(pathlonsteps[i,:])): 

                Plat=pathlatsteps[i,p]
                Plon=pathlonsteps[i,p]
                if (Plat < S_lt) or (Plat < W_lt and Plon < W_ln):
                    E_n[i,p] = en_int_N(pathlonsteps[i,p], pathlatsteps[i,p])
                    E_e[i,p] = ee_int_N(pathlonsteps[i,p], pathlatsteps[i,p])

                else:
                    E_n[i,p] = en_int_D(pathlonsteps[i,p], pathlatsteps[i,p])
                    E_e[i,p] = ee_int_D(pathlonsteps[i,p], pathlatsteps[i,p])             

            # Integrate to get V = int(E*dL), use cylindrical coordinates:
            intline = dists[nodefrom[i], nodeto[i]]
            intazi = azi[nodefrom[i], nodeto[i]]
            steps = number_steps_lon #intline * 3   # int

            for j in range(int(steps)-1):
                # For a North field:
                vnseg = (0.5 * ( E_n[i,j] + E_n[i,j+1] ) * cos(intazi)*(intline/steps) )   # comparar estes vetores
                # For an East field:
                veseg = (0.5 * ( E_e[i,j] + E_e[i,j+1] ) * sin(intazi)*(intline/steps) )
                
                
                Vn_tot[i] += vnseg
                Ve_tot[i] += veseg


    # ===============================================================
    # 6) USE V AND SYSTEM MATRIX TO DETERMINE CURRENT (J)
    # ===============================================================
    # Create Nvoltage matrix for northward field:
    Nvoltage = np.zeros((nnodes, nnodes), dtype=npf)
    for l in range(nconnects):
        Nvoltage[nodefrom[l],nodeto[l]] = Vn_tot[l]
        Nvoltage[nodeto[l],nodefrom[l]] = -Vn_tot[l]


    # Create Evoltage matrix for eastward field:
    Evoltage = np.zeros((nnodes, nnodes), dtype=npf)
    for l in range(nconnects):
        Evoltage[nodefrom[l],nodeto[l]] = Ve_tot[l]
        Evoltage[nodeto[l],nodefrom[l]] = -Ve_tot[l]
    
    results_fem.append([Vn_tot, Ve_tot])
    
    # Nvoltage and Evoltage were not made symmetric. 
    # Making them symmetric makes it simpler to compute the free current vector (see below)
    #                               A.P.

    # Use Ohm's law to calculate the current along each pathlength:
    # LP1985 eq. (14): J = V/R
    Nlinecurr = np.zeros((nnodes, nnodes), dtype=npf)
    for m in range(nnodes):
        for n in range(nnodes):
            if resis[m,n] > 0.:
                newval = Nvoltage[m,n]/resis[m,n]
                if not np.isnan(newval):
                    Nlinecurr[m,n] = newval
                else:
                    Nlinecurr[m,n] = 0.

    Elinecurr = np.zeros((nnodes, nnodes), dtype=npf)
    for m in range(nnodes):
        for n in range(nnodes):
            if resis[m,n] > 0.:
                newval = Evoltage[m,n]/resis[m,n]
                if not np.isnan(newval):
                    Elinecurr[m,n] = newval
                else:
                    Elinecurr[m,n] = 0.

    # Nlinecurr and Elinecurr are not made symmetric
    # Making them symmetric makes it simpler to compute the free current vector (see below)
    #                               A.P.


    # Total line current from both components:
    Tlinecurr = Nlinecurr + Elinecurr

    # Total current at each node due to the E field:
    # By making the Nsourcevec and Esourcevec matrices symmetric, the computation of 
    # the free current vector is made simpler (no need to use two matrices, only one)
    #                               A.P.
    Nsourcevec = np.sum(Nlinecurr, axis=0) 
    Esourcevec = np.sum(Elinecurr, axis=0)

    # GIC formed at each node due to the E field:
    netconN = np.dot(np.linalg.inv(systemmat), Nsourcevec)
    netconE = np.dot(np.linalg.inv(systemmat), Esourcevec)
    netconT = netconN + netconE

#****************************************************************

    results = [netconN, netconE]                # V PER NODE ! not GIC
    global_results.append([netconN, netconE]) 
    current_sources_results.append([Nlinecurr, Elinecurr])    

#****************************************************************
#COPY FROM JAR
#****************************************************************
transf_results=[]
subs_results=[]
# This calculation may be required later.

# Calculate current flowing through transformers (NOT line current):
# ALT jaribeiro
for k in range(0,len(global_results)):
    Isub_N = np.multiply(global_results[k][0],adm_earth)
    Isub_E = np.multiply(global_results[k][1],adm_earth)
    subs_results.append([Isub_N, Isub_E])
    Iline_N, Iline_E = np.zeros(nconnects), np.zeros(nconnects)

# Note that the algebraic sign of currents
# is + for current from nodefrom[l] to nodeto[l].
# It means that, if we want to obtain signs according to the standard convention,
# we have to guarantee in the connections input file that, in transformers, connection 
# always goes TO the earthed node (gray nodes)
#                         A.P.

    for l in range(0,nconnects):
        m=nodefrom[l]
        n=nodeto[l]
        newval = 1/res_line[l]
        if not np.isinf(newval):
            # i_mn = j_mk + (v_n − v_m)y_nm 
            Iline_N[l] = current_sources_results[k][0][m,n] + np.dot((global_results[k][0][m]-global_results[k][0][n]),(1/res_line[l])) 
            Iline_E[l] = current_sources_results[k][1][m,n] + np.dot((global_results[k][1][m]-global_results[k][1][n]),(1/res_line[l])) 
        else:
            Iline_N[l] = 0
            Iline_E[l] = 0


#****************************************************************
    transresults = [Iline_N, Iline_E]   # GIC PER TRANSFORMER
# ALT jaribeiro
    transf_results.append([Iline_N, Iline_E])
#****************************************************************


# ===============================================================
# 7) PRINT RESULTS ON SCREEN
# ===============================================================
#
# The list global_results is converted into a numpy array in the end
# dimensions are (nfr,2,nnodes)
#

#Alt jaribeiro
GICsubs=np.array(subs_results)
# ALT jaribeiro
GICtran=np.array(transf_results)

# ===============================================================
# 8) SAVE GIC RESULTS ON PICKLE FILE
# ===============================================================
#
        
with open(SavePath+'GICs_'+ out_name+'.pkl', 'wb') as f:
        pickle.dump([sitename,GICsubs],f)

# #Alt jaribeiro
with open(SavePath+'GICtrf_'+out_name + '.pkl', 'wb') as f:
            pickle.dump([substation,TransType,Wind,GICtran],f) 

