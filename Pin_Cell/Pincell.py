##############################################################################################################################################################################################################################################################################################################################################################################################
# FYS4580 PROJECT H24 - PART 1                                                                                                                                                       # COMMENTS
##############################################################################################################################################################################################################################################################################################################################################################################################

import openmc # type: ignore
import openmc.deplete # type: ignore
import openmc.model # type: ignore
import matplotlib.pyplot as plt
import numpy as np

##############################################################################################################################################################################################################################################################################################################################################################################################
# MATERIALS
##############################################################################################################################################################################################################################################################################################################################################################################################

uranium_enrichment = 0.048
urox = openmc.Material(name = "Fuel UO2")                                                                                                                                   # Defining Uranium dioxide material
urox.add_nuclide("U235", uranium_enrichment)                                                                                                                                              # ("name", percentage)
urox.add_nuclide("U238", 1 - uranium_enrichment)
urox.add_element("O", 2)                                                                                                                                                    # Ignore natural distribution of oxigen isotopes (simplification))
urox.set_density("g/cm3", 10.98)
urox.temperature = 1474

helium = openmc.Material(name = "Helium")                                                                                                                                   # Defining Helium material
helium.add_element("He", 1)
helium.set_density("g/cm3", 0.00000178)
helium.temperature = 691

zircaloy4 = openmc.Material(name = "Zircaloy")                                                                                                                              # Defining Zircaloy material (simplified)
zircaloy4.add_element("Zr", 0.985)
zircaloy4.add_element("Sn", 0.015)
zircaloy4.set_density("g/cm3", 6.56)
zircaloy4.temperature = 691

graphite = openmc.Material(name = "Graphite")
graphite.add_element("C", 1)
graphite.set_density("g/cm3", 1.85)
graphite.add_s_alpha_beta("c_Graphite")

light_water = openmc.Material(name = "Water")                                                                                                                               # Defining Water material
light_water.add_element("H", 2)
light_water.add_element("O", 1)
light_water.set_density("g/cm3", 1)
light_water.temperature = 902

heavy_water = openmc.Material(name = "Heavy Water")                                                                                                                               # Defining Water material
heavy_water.add_nuclide("H2", 2)
heavy_water.add_element("O", 1)
heavy_water.set_density("g/cm3", 1.107)
heavy_water.temperature = 902

B10_concentration = 0.2
boric_acid = openmc.Material(name = "B(OH)3")
boric_acid.add_nuclide("B10", B10_concentration)
boric_acid.add_nuclide("B11", 1 - B10_concentration)
boric_acid.add_element("O", 3)
boric_acid.add_element("H", 3)
boric_acid.set_density("g/cm3", 1.435)
boric_acid.temperature = 902

Boric_acid_concentration = 0.01
borated_water = openmc.Material.mix_materials([light_water, boric_acid], [1 - Boric_acid_concentration, Boric_acid_concentration], 'wo')
borated_water.add_s_alpha_beta("c_H_in_H2O")


##############################################################################################################################################################################################################################################################################################################################################################################################
# GEOMETRY
##############################################################################################################################################################################################################################################################################################################################################################################################

# Fuel Rod Diameter
frd = 0.95 + 0.05                                                                                                                                                                 # All of the reactor geometry dimentions will be decided as a function of the fuel rod diamater

# Define Boundary Type
boundary = "reflective" #"vacuum" 

# Top and Bottom Planes
z_min = openmc.ZPlane(-1.83, boundary_type = boundary)                                                                                                                         # Defining a minimum amd a maximum in the z-plane which will determine the height of the cylinders and applying boundary conditions: transmissive, vacuum or reflektive (for neutrons)
z_max = openmc.ZPlane(1.83, boundary_type = boundary)

# Fuel Rods
fuel_pin = -openmc.ZCylinder(r = frd/2 - 0.05, boundary_type = "transmission")                                                                                                                                       # Defining cylinders with infinte z axis and r (always given in cm)
helium_gap = -openmc.ZCylinder(r = frd/2 - 0.03, boundary_type = "transmission")                                                                                                                                     # - sign defines interior of cylinder, + sign defines exterior. Helium_gap is a thin cylindrical shell 2mm thick, so r = 0.40 + 0.02
cladding_region = -openmc.ZCylinder(r = frd/2, boundary_type = "transmission")                                                                                                                                 # r = 0.40 + 0.02 + 0.03

fuel_pin = fuel_pin & +z_min & -z_max                                                                                                                                       # Defining the cylinders' height / + sign means upwards while - sign means downwards
helium_gap = helium_gap & +z_min & -z_max
cladding_region = cladding_region & +z_min & -z_max

helium_gap = helium_gap & ~fuel_pin                                                                                                                                         # Making the cylindres hollow, ~ sign removes the following volume
cladding_region = cladding_region & ~helium_gap & ~fuel_pin

# Outer box enviroment
box = -openmc.model.RectangularPrism(width = 2*0.95, height = 2*0.95, boundary_type = boundary)                                                                                   # Defining an external box
box = box & +z_min & -z_max

fuel_box = box & ~cladding_region & ~helium_gap & ~fuel_pin                                                                                                                 # None of the geometry elements is overlapping

##############################################################################################################################################################################################################################################################################################################################################################################################
# CELLS
##############################################################################################################################################################################################################################################################################################################################################################################################

# Fuel Rods Cells
fuel_cell = openmc.Cell(fill = urox, region = fuel_pin)                                                                                                                     # Filling the geometry elements with fuel material
helium_cell = openmc.Cell(fill = helium, region = helium_gap)
cladding_cell = openmc.Cell(fill = zircaloy4 , region = cladding_region)
water_cell = openmc.Cell(fill = borated_water, region = fuel_box)

cell_list_fuel = [fuel_cell, helium_cell, cladding_cell, water_cell] 

fuel_universe = openmc.Universe(cells = cell_list_fuel)                                                                                                                      # Need to put all the geometry elements in a "universe" to be sure that the program will run and the elements are well defined, it's not allowed to put the same cell element in two different "universes" 

##############################################################################################################################################################################################################################################################################################################################################################################################
# SIMULATION SETTINGS 
##############################################################################################################################################################################################################################################################################################################################################################################################

materials = openmc.Materials([urox, helium, zircaloy4, light_water, heavy_water, graphite, boric_acid, borated_water])
materials.export_to_xml()                                                                                                                                                   # OpenMc does not run a python script so everything is saved into xml files. (The smallest set of xml files neede is: materials.xml, geometry.xml og settings.xml)

geometry = openmc.Geometry()
geometry.root_universe = fuel_universe
geometry.export_to_xml()

settings = openmc.Settings()
settings.batches = 100                                                                                                                                                      # How many iterations/batches of the same neutron generation to simulate
settings.inactive = 10                                                                                                                                                      # Number of inactive batches = How many iterations are used to determine sample spots (MC simulation)
settings.particles = 10000                                                                                                                                                  # Number of neutrons in each generation
settings.temperature = {"method":"interpolation", "tolerance":100}                                                                                                                                   # Number of generations per batch (not used here)

##############################################################################################################################################################################################################################################################################################################################################################################################
# OpenMc PARTICLE SAMPLE SITES 
##############################################################################################################################################################################################################################################################################################################################################################################################

bounds = [-1,-1,-3,1,1,3]           
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])                                                                                         
settings.source = openmc.IndependentSource(space=uniform_dist, constraints={'fissionable': True})                                                                           # Generating neutrons only in fissible material (e.g. not in water but in urox)            
settings.export_to_xml()

##############################################################################################################################################################################################################################################################################################################################################################################################
# PLOTS AND PICTURES
##############################################################################################################################################################################################################################################################################################################################################################################################
'''
fuel_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", graphite : "midnightblue"})
plt.savefig("Pincell_plot_xy")

fuel_universe.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", graphite : "midnightblue"})
plt.savefig("Pincell_plot_yz")
'''
##############################################################################################################################################################################################################################################################################################################################################################################################
# TALLY ARITHMETICS 
##############################################################################################################################################################################################################################################################################################################################################################################################

# Defining Tallies File
tallies_file = openmc.Tallies()                                                                                                                                             # For our tallies to be included we have to put them into an inp

# Defining a Spatial Mesh for Fuel Pin Cell
pin_mesh = openmc.RegularMesh()                                                                                                                                                 # Making a mesh instance
pin_mesh.dimension = [100,100]                                                                                                                                                  # Setting the bins of the mesh to 100 x 100 bins
pin_mesh.lower_left = [-1,-1]                                                                                                
pin_mesh.upper_right = [1, 1]                                                                                                  
pin_mesh_filter = openmc.MeshFilter(pin_mesh)

# Spatial Distribution of the Prompt Fission Sites and the Spatial Neutron Flux Distribution for Fuel Pin Cell
pin_tally = openmc.Tally(name = "Pin Tally")                                                                                                                                        # Making a tally instance
pin_tally.filters = [pin_mesh_filter]                                                                                                                                               # Designatingen the mesh for the tally
pin_tally.scores = ["flux", 'prompt-nu-fission']                                                                                                                                # Designating the objects to be tallied
tallies_file.append(pin_tally)  

# Energy Distribution of the Neutron Flux in the Fuel on a log-log Scale.
energy_filter = openmc.EnergyFilter(np.logspace(np.log10(10e-5),np.log10(20e6),501))
energy = openmc.Tally(name = "Energy")                                                                                                                                      # Making a new tally only looking at flux
energy.filters = [energy_filter]                                                                                                                                            # Designating the energy filter
energy.scores = ["flux"]
tallies_file.append(energy)

# Five-Factor Formula Coefficients
# --------------------------------

# Fast Fission Factor tallies
therm_fiss_rate = openmc.Tally(name='therm. fiss. rate')
therm_fiss_rate.scores = ['nu-fission']
therm_fiss_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_fiss_rate)

# Resonance Escape Probability tallies
therm_abs_rate = openmc.Tally(name='therm. abs. rate')
therm_abs_rate.scores = ['absorption']
therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625])]
tallies_file.append(therm_abs_rate)

# Thermal Flux Utilization tallies
fuel_therm_abs_rate = openmc.Tally(name='fuel therm. abs. rate')
fuel_therm_abs_rate.scores = ['absorption']
fuel_therm_abs_rate.filters = [openmc.EnergyFilter([0., 0.625]), openmc.CellFilter([fuel_cell])]
tallies_file.append(fuel_therm_abs_rate)

# Instantiate a tally mesh
mesh = openmc.RegularMesh(mesh_id=1)
mesh.dimension = [1, 1, 1]
mesh.lower_left = [-1,-1, -3]
mesh.upper_right= [1,1,3]
#mesh.width = [1.26, 1.26, 200.]
meshsurface_filter = openmc.MeshSurfaceFilter(mesh)

# Instantiate thermal, fast, and total leakage tallies
leak = openmc.Tally(name='leakage')
leak.filters = [meshsurface_filter]
leak.scores = ['current']
tallies_file.append(leak)

thermal_leak = openmc.Tally(name='thermal leakage')
thermal_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0., 0.625])]
thermal_leak.scores = ['current']
tallies_file.append(thermal_leak)

fast_leak = openmc.Tally(name='fast leakage')
fast_leak.filters = [meshsurface_filter, openmc.EnergyFilter([0.625, 20.0e6])]
fast_leak.scores = ['current']
tallies_file.append(fast_leak)

# K-Eigenvalue (infinity) tallies
fiss_rate = openmc.Tally(name='fiss. rate')
abs_rate = openmc.Tally(name='abs. rate')
fiss_rate.scores = ['nu-fission']
abs_rate.scores = ['absorption']
tallies_file += (fiss_rate, abs_rate)

# Beta-Delayed Neutrons vs Total Neutrons
# ---------------------------------------

# Beta-Delayed neutron fission
delayed_nu = openmc.Tally(name="del.nu.fis.")
delayed_nu.scores = ['delayed-nu-fission']
tallies_file.append(delayed_nu)

# Total fission
total_nu = openmc.Tally(name="tot.nu.fis.")
total_nu.scores = ['nu-fission']
tallies_file.append(total_nu)

tallies_file.export_to_xml()

##############################################################################################################################################################################################################################################################################################################################################################################################
# RUNNING OpenMC 
##############################################################################################################################################################################################################################################################################################################################################################################################

if __name__ == "__main__":
    openmc.run()  #threads = 4                                                                                                                                               # Threads set for heplab primarily (to limit utilized cores on heplab server, this can be removed when run locally)

##############################################################################################################################################################################################################################################################################################################################################################################################
# END
##############################################################################################################################################################################################################################################################################################################################################################################################
