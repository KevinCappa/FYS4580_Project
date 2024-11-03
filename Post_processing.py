#####################################################################################################################################################
# FYS4580 PROJECT H24 - POST PROCESSING
#####################################################################################################################################################

import openmc # type: ignore
import openmc.deplete # type: ignore
import openmc.model # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from OpenMC_Project import temperatures

#####################################################################################################################################################
# NEUTRON DISTRIBUTION MESH FILTER
#####################################################################################################################################################

# Fuel Pin Cell
sp = openmc.StatePoint('statepoint.100.h5')
tally = sp.tallies[1]
flux = tally.get_slice(scores=['flux'])
prompt = tally.get_slice(scores=['prompt-nu-fission'])
flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
prompt.std_dev.shape = (100, 100)
prompt.mean.shape = (100, 100)

fig = plt.figure(figsize=(16, 12), dpi=500)

ax1 = plt.subplot(121)
ax1.set_title('Neutron flux distribution', size=25, y=1.05)
im1 = ax1.imshow(flux.mean)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad = 0.05)
fig.colorbar(im1, cax=cax, orientation = 'vertical')

ax2 = plt.subplot(122)
ax2.set_title('Fission site heatmap', size=25, y=1.05)
im2 = ax2.imshow(prompt.mean)

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad = 0.05)
fig.colorbar(im2, cax=cax, orientation = 'vertical')

plt.savefig("Pictures/Mesh_Fuel_Pin_Cell", dpi=500)

# Reactor Pressure Vessel
sp = openmc.StatePoint('statepoint.100.h5')
tally = sp.tallies[2]
flux = tally.get_slice(scores=['flux'])
prompt = tally.get_slice(scores=['prompt-nu-fission'])
flux.std_dev.shape = (100, 100)
flux.mean.shape = (100, 100)
prompt.std_dev.shape = (100, 100)
prompt.mean.shape = (100, 100)

fig = plt.figure(figsize=(16, 12), dpi=500)

ax1 = plt.subplot(121)
ax1.set_title('Neutron flux distribution', size=25, y=1.05)
im1 = ax1.imshow(flux.mean)

divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='5%', pad = 0.05)
fig.colorbar(im1, cax=cax, orientation = 'vertical')

ax2 = plt.subplot(122)
ax2.set_title('Fission site heatmap', size=25, y=1.05)
im2 = ax2.imshow(prompt.mean)

divider = make_axes_locatable(ax2)
cax = divider.append_axes('right', size='5%', pad = 0.05)
fig.colorbar(im2, cax=cax, orientation = 'vertical')

plt.savefig("Pictures/Mesh_RPV", dpi=500)

#####################################################################################################################################################
# NEUTRON ENERGY DISTRIBUTION
#####################################################################################################################################################

erange = np.logspace(base=10, start=-5, stop=7, num=501)

tally2 = sp.tallies[3]
flx = tally2.mean.ravel()
plt.figure()
plt.loglog(erange[:-1], flx)
plt.grid()
plt.xlabel("Energy eV")
plt.ylabel("Flux [n/cm-src]")
plt.title("Neutron energy spectrum")

plt.savefig("Pictures/Energy")

#####################################################################################################################################################
# KEFF - FIVE-FACTOR FORMULA
#####################################################################################################################################################

# Get the fission and absorption rate tallies
fiss_rate = sp.get_tally(name='fiss. rate')
abs_rate = sp.get_tally(name='abs. rate')

# Get the leakage tally
leak = sp.get_tally(name='leakage')
leak = leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)

# (This method gives somehow the wrong keff)

# # Compute k-infinity using tally arithmetic                                                  
# keff = fiss_rate / (abs_rate + leak)
# keffcalc = keff.get_pandas_dataframe()
# print(keffcalc.iloc[0]["mean"], keffcalc.iloc[0]["std. dev."])

# FAST FISSION FACTOR
therm_fiss_rate = sp.get_tally(name='therm. fiss. rate')
fast_fiss = fiss_rate / therm_fiss_rate
epsilon = fast_fiss.get_pandas_dataframe()
epsilon_value = round(epsilon.iloc[0]["mean"],5)
epsilon_value_err = round(epsilon.iloc[0]["std. dev."],5)

# RESONANCE ESCAPE PROBABILITY
therm_abs_rate = sp.get_tally(name='therm. abs. rate')
thermal_leak = sp.get_tally(name='thermal leakage')
thermal_leak = thermal_leak.summation(filter_type=openmc.MeshSurfaceFilter, remove_filter=True)
res_esc = (therm_abs_rate + thermal_leak) / (abs_rate + thermal_leak)
p = res_esc.get_pandas_dataframe()
p_value = round(p.iloc[0]["mean"],5)
p_value_err = round(p.iloc[0]["std. dev."],5)

# THERMAL UTILIZATION FACTOR
fuel_therm_abs_rate = sp.get_tally(name='fuel therm. abs. rate')
therm_util = fuel_therm_abs_rate / therm_abs_rate
f = therm_util.get_pandas_dataframe()
f_value = round(f.iloc[0]["mean"],5)
f_value_err = round(f.iloc[0]["std. dev."],5)

# NEUTRON REPRODUCTION FACTOR
eta = therm_fiss_rate / fuel_therm_abs_rate
eta = eta.get_pandas_dataframe()
eta_value = round(eta.iloc[0]["mean"],5)
eta_value_err = round(eta.iloc[0]["std. dev."],5)

# NON-LEAKAGE PROBABILITY
# Fast Neutrons
p_fnl = (abs_rate + thermal_leak) / (abs_rate + leak)
PNL_F = p_fnl.get_pandas_dataframe()
PNL_F_value = round(PNL_F.iloc[0]["mean"],5)
PNL_F_value_err = round(PNL_F.iloc[0]["std. dev."],5)
# Thermal Neutrons
p_tnl = therm_abs_rate / (therm_abs_rate + thermal_leak)
PNL_T = p_tnl.get_pandas_dataframe()
PNL_T_value = round(PNL_T.iloc[0]["mean"],5)
PNL_T_value_err = round(PNL_T.iloc[0]["std. dev."],5)

# KEFF
Keff_value = p_value * epsilon_value * f_value * eta_value* PNL_F_value * PNL_T_value

print(f"Fast fission factor = {epsilon_value}",
      f"Resonance escape probability = {p_value}",
      f"Thermal utilization factor = {f_value}",
      f"Neutron reproduction factor = {eta_value}",
      f"Non-leakege probability (fast) = {PNL_F_value}",
      f"Non-leakege probability (thermal) = {PNL_T_value}",
      f"Calculated Neutron multiplication number, k = {Keff_value:.5f}", sep="\n")

#####################################################################################################################################################
# BETA-DELAYED NEUTRONS
#####################################################################################################################################################

del_beta_n = sp.get_tally(name='del.nu.fis.')
tot_fission_n = sp.get_tally(name='tot.nu.fis.')
beta = del_beta_n / tot_fission_n
beta = beta.get_pandas_dataframe()
beta_value = round(beta.iloc[0]["mean"],5)
beta_value_err = round(beta.iloc[0]["std. dev."],5)

print(f"The number of beta-delayed neutrons are {beta_value:.7f} +\- {beta_value_err:.7f}")

#####################################################################################################################################################
# CVS DATA FILE
#####################################################################################################################################################

keff = sp.keff                               
keff_nom = keff.nominal_value
keff_std = keff.std_dev

# data = [keff_nom, keff_std, Keff_value, epsilon_value, p_value, f_value, eta_value, PNL_F_value,  PNL_T_value]
# csvfile = open('data.csv', 'a') 
# csvwriter = csv.writer(csvfile)
# csvwriter.writerow(data)
# csvfile.close()

print(f"OpenMC Simulated Neutron Multiplication Number, k effective = {keff_nom:.5f} +\- {keff_std:.5f}")
print(f"Difference (Simulation k - Calculated k) = {keff_nom-Keff_value:.5f}")

#####################################################################################################################################################
# DOPPLER BROADENING
#####################################################################################################################################################

x = [] 
y = [] 
y_error = []
  
with open('temperatures.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',')
    next(csvfile) 
    for row in lines: 
        x.append(row[0])

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y.append(round(float(row[0]),5))

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_error.append(round(float(row[1]),5))

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y, color = 'g', linestyle = 'dashed', marker = 'o',label = "Fuel Temperature Data") 
plt.errorbar(x, y, yerr = y_error, fmt ='o', capsize=3)  
plt.xticks(rotation = 25) 
plt.xlabel('Temperature(°K)', size=25, labelpad=25) 
plt.ylabel('Neutron Multiplication Number k', size=25, labelpad=25) 
plt.title('Doppler Broadening', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Pictures/Doppler") 

'''
# Doppler Broadening for Five-Factor Formula

y_Keff = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_Keff.append(round(float(row[2]),5))

y_epsilon = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_epsilon.append(round(float(row[3]),5))

y_p = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_p.append(round(float(row[4]),5))

y_f = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_f.append(round(float(row[5]),5))

y_eta = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_eta.append(round(float(row[6]),5))

y_PNLF = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_PNLF.append(round(float(row[7]),5))

y_PNLT = []

with open('data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    for row in lines: 
        y_PNLT.append(round(float(row[8]),5))

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y_Keff, color = 'green', linestyle = 'solid', marker = 'o',label = "Calculated k") 
plt.plot(x, y_epsilon, color = 'red', linestyle = 'solid', marker = 'o',label = "Fast fission factor") 
plt.plot(x, y_p, color = 'blue', linestyle = 'dotted', marker = 'o',label = "Resonance escape probability") 
plt.plot(x, y_f, color = 'yellow', linestyle = 'solid', marker = 'o',label = "Thermal utilization factor") 
plt.plot(x, y_eta, color = 'black', linestyle = 'solid', marker = 'o',label = "Neutron reproduction number") 
plt.plot(x, y_PNLF, color = 'pink', linestyle = 'solid', marker = 'o',label = "Non-Leakage probability (fast)") 
plt.plot(x, y_PNLT, color = 'orange', linestyle = 'solid', marker = 'o',label = "Non-Leakage probability (thermal)") 
plt.xticks(rotation = 25) 
plt.xlabel('Temperature(°K)', size=25, labelpad=25) 
plt.ylabel('Five-Factor Formula', size=25, labelpad=25) 
plt.title('Doppler Broadening', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Pictures/Doppler_Five")  
'''

#####################################################################################################################################################
# DEPLETION
#####################################################################################################################################################

results = openmc.deplete.Results("./depletion_results.h5")
time, k = results.get_keff()
time /= (24 * 60 * 60)  # convert back to days from seconds
plt.figure()
plt.errorbar(time, k[:, 0], yerr=k[:, 1], capsize=3, fmt="blue", marker = 'o', ecolor = "black")
plt.xlabel("Time [d]")
plt.ylabel("$k_{eff}\pm \sigma$")
plt.grid() 
# plt.legend(fontsize = 10) 
plt.savefig("Pictures/Depletion_keff",  dpi=500)

_time, u235 = results.get_atoms("1", "U235")
_time, pu239 = results.get_atoms("1", "Pu239")
plt.figure()
plt.plot(time, u235, marker = 'o', label="U235")
plt.plot(time, pu239, marker = 'o', label="Pu239")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - U235 and Pu239")
plt.grid() 
plt.legend(fontsize = 10) 
plt.savefig("Pictures/Depletion_U235_P239")

_time, u238 = results.get_atoms("1", "U238")
_time, pu241 = results.get_atoms("1", "Pu241")
plt.figure()
plt.plot(time, u238, marker = 'o', label="U238")
plt.plot(time, pu241, marker = 'o', label="Pu241")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - U238 and Pu241")
plt.grid() 
plt.legend(fontsize = 10) 
plt.savefig("Pictures/Depletion_U238_Pu241")

_time, xe135 = results.get_atoms("1", "Xe135")
_time, sm149 = results.get_atoms("1", "Sm149")
plt.figure()
plt.plot(time, xe135, marker = 'o', label="Xe135")
plt.plot(time, sm149, marker = 'o', label="Sm149")
plt.xlabel("Time [d]")
plt.ylabel("Number of atoms - Xe135 and Sm149")
plt.grid() 
plt.legend(fontsize = 10) 
plt.savefig("Pictures/Depletion_Xe135_Sm149")

_time, u235_fission = results.get_reaction_rate("1", "U235", "fission")
_time, u238_fission = results.get_reaction_rate("1", "U238", "fission")
_time, Pu239_fission = results.get_reaction_rate("1", "Pu239", "fission")
_time, Pu241_fission = results.get_reaction_rate("1", "Pu241", "fission")
plt.figure()
plt.plot(time, u235_fission, marker = 'o', label="U235")
plt.plot(time, u238_fission, marker = 'o', label="U238")
plt.plot(time, Pu239_fission, marker = 'o', label="Pu239")
plt.plot(time, Pu241_fission, marker = 'o', label="Pu241")
plt.xlabel("Time [d]")
plt.ylabel("Fission reactions / s")
plt.grid() 
plt.legend(fontsize = 10) 
plt.savefig("Pictures/Depletion_Fission")

#####################################################################################################################################################
# END
#####################################################################################################################################################

# TO DO LIST: 
# Make better plots
# Run Depletion again
# Write report