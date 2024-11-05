#####################################################################################################################################################
# FYS4580 PROJECT H24 - POST PROCESSING
#####################################################################################################################################################

import openmc # type: ignore
import openmc.model # type: ignore
import matplotlib.pyplot as plt # type: ignore
import numpy as np
import csv
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

plt.savefig("Mesh_Fuel_Pin_Cell", dpi=500)

#####################################################################################################################################################
# NEUTRON ENERGY DISTRIBUTION
#####################################################################################################################################################

erange = np.logspace(base=10, start=-5, stop=7, num=501)

tally2 = sp.tallies[2]
flx = tally2.mean.ravel()
plt.figure()
plt.loglog(erange[:-1], flx)
plt.grid()
plt.xlabel("Energy eV")
plt.ylabel("Flux [n/cm-src]")
plt.title("Neutron energy spectrum")

plt.savefig("Energy")

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
# csvfile = open('pincell_data.csv', 'a') 
# csvwriter = csv.writer(csvfile)
# csvwriter.writerow(data)
# csvfile.close()

print(f"OpenMC Simulated Neutron Multiplication Number, k effective = {keff_nom:.5f} +\- {keff_std:.5f}")
print(f"Difference (Simulation k - Calculated k) = {keff_nom-Keff_value:.5f}")

#####################################################################################################################################################
# PLOTS FOR VARIABLE PIN SIZES
#####################################################################################################################################################

x = [0.75, 0.80, 0.85, 0.90, 0.95, 1.00] 

# Light Water Plot 
y_Keff = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            y_Keff.append(round(float(row[0]),5))
            rownum +=1
      else: break
      

y_epsilon = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_epsilon.append(round(float(row[3]),5))
            rownum +=1
      else: break

y_p = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_p.append(round(float(row[4]),5))
            rownum +=1
      else: break

y_f = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_f.append(round(float(row[5]),5))
            rownum +=1
      else: break

y_eta = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile) 
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_eta.append(round(float(row[6]),5))
            rownum +=1
      else: break

y_PNL = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_PNL.append(round(float(row[7]),5))
            rownum +=1
      else: break

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y_Keff, color = 'green', linestyle = 'solid', marker = 'o',label = "Calculated k") 
plt.plot(x, y_epsilon, color = 'red', linestyle = 'solid', marker = 'o',label = "Fast fission factor") 
plt.plot(x, y_p, color = 'blue', linestyle = 'dotted', marker = 'o',label = "Resonance escape probability") 
plt.plot(x, y_f, color = 'yellow', linestyle = 'solid', marker = 'o',label = "Thermal utilization factor") 
plt.plot(x, y_eta, color = 'black', linestyle = 'solid', marker = 'o',label = "Neutron reproduction number") 
plt.plot(x, y_PNL, color = 'pink', linestyle = 'solid', marker = 'o',label = "Non-Leakage probability (fast)") 
plt.xticks(rotation = 25) 
plt.xlabel('Fuel Pin Size [cm]', size=25, labelpad=25) 
plt.ylabel('Five-Factor Formula', size=25, labelpad=25) 
plt.title('Size Impact on k-effective', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Light_Water")

# Light Water Plot for k-effective and Thermal utilization factor
y_Keff = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            y_Keff.append(round(float(row[0]),5))
            rownum +=1
      else: break

y_f = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7: 
            y_f.append(round(float(row[5]),5))
            rownum +=1
      else: break

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y_Keff, color = 'red', linestyle = 'solid', marker = 'o',label = "k-effective") 
plt.plot(x, y_f, color = 'blue', linestyle = 'solid', marker = 'o',label = "Thermal utilization factor") 
plt.xticks(rotation = 25) 
plt.xlabel('Fuel Pin Size [cm]', size=25, labelpad=25) 
plt.ylabel('Five-Factor Formula', size=25, labelpad=25) 
plt.title('k-effective and Thermal utilization factor', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Light_Water_kandf")

# Heavy Water Plot 
y_Keff = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13: 
            y_Keff.append(round(float(row[0]),5))
            rownum +=1
      else: break
      
y_epsilon = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13:  
            y_epsilon.append(round(float(row[3]),5))
            rownum +=1
      else: break

y_p = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13:  
            y_p.append(round(float(row[4]),5))
            rownum +=1
      else: break

y_f = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13:  
            y_f.append(round(float(row[5]),5))
            rownum +=1
      else: break

y_eta = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile) 
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13:  
            y_eta.append(round(float(row[6]),5))
            rownum +=1
      else: break

y_PNL = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <7:
            rownum +=1
      elif rownum <13: 
            y_PNL.append(round(float(row[7]),5))
            rownum +=1
      else: break

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y_Keff, color = 'green', linestyle = 'solid', marker = 'o',label = "Calculated k") 
plt.plot(x, y_epsilon, color = 'red', linestyle = 'solid', marker = 'o',label = "Fast fission factor") 
plt.plot(x, y_p, color = 'blue', linestyle = 'dotted', marker = 'o',label = "Resonance escape probability") 
plt.plot(x, y_f, color = 'yellow', linestyle = 'solid', marker = 'o',label = "Thermal utilization factor") 
plt.plot(x, y_eta, color = 'black', linestyle = 'solid', marker = 'o',label = "Neutron reproduction number") 
plt.plot(x, y_PNL, color = 'pink', linestyle = 'solid', marker = 'o',label = "Non-Leakage probability (fast)") 
plt.xticks(rotation = 25) 
plt.xlabel('Fuel Pin Size [cm]', size=25, labelpad=25) 
plt.ylabel('Five-Factor Formula', size=25, labelpad=25) 
plt.title('Size Impact on k-effective', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Heavy_Water")

# Graphite Plot 
y_Keff = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19: 
            y_Keff.append(round(float(row[0]),5))
            rownum +=1
      else: break
      
y_epsilon = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19:  
            y_epsilon.append(round(float(row[3]),5))
            rownum +=1
      else: break

y_p = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19:  
            y_p.append(round(float(row[4]),5))
            rownum +=1
      else: break

y_f = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19:  
            y_f.append(round(float(row[5]),5))
            rownum +=1
      else: break

y_eta = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile) 
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19:  
            y_eta.append(round(float(row[6]),5))
            rownum +=1
      else: break

y_PNL = []

with open('pincell_data.csv','r') as csvfile: 
    lines = csv.reader(csvfile, delimiter=',') 
    next(csvfile)
    rownum = 1
    for row in lines:
      if rownum <13:
            rownum +=1
      elif rownum <19: 
            y_PNL.append(round(float(row[7]),5))
            rownum +=1
      else: break

plt.figure(figsize=(16, 12), dpi=500) 
plt.plot(x, y_Keff, color = 'green', linestyle = 'solid', marker = 'o',label = "Calculated k") 
plt.plot(x, y_epsilon, color = 'red', linestyle = 'solid', marker = 'o',label = "Fast fission factor") 
plt.plot(x, y_p, color = 'blue', linestyle = 'dotted', marker = 'o',label = "Resonance escape probability") 
plt.plot(x, y_f, color = 'yellow', linestyle = 'solid', marker = 'o',label = "Thermal utilization factor") 
plt.plot(x, y_eta, color = 'black', linestyle = 'solid', marker = 'o',label = "Neutron reproduction number") 
plt.plot(x, y_PNL, color = 'pink', linestyle = 'solid', marker = 'o',label = "Non-Leakage probability (fast)") 
plt.xticks(rotation = 25) 
plt.xlabel('Fuel Pin Size [cm]', size=25, labelpad=25) 
plt.ylabel('Five-Factor Formula', size=25, labelpad=25) 
plt.title('Size Impact on k-effective', fontsize = 40, y=1.05) 
plt.grid() 
plt.legend(fontsize = 20) 
plt.savefig("Graphite")

#####################################################################################################################################################
# END
#####################################################################################################################################################