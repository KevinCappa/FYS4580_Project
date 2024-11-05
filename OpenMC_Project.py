##############################################################################################################################################################################################################################################################################################################################################################################################
# FYS4580 PROJECT - PART 2                                                                                                                                                       # COMMENTS
##############################################################################################################################################################################################################################################################################################################################################################################################

import openmc # type: ignore
import openmc.deplete # type: ignore
import openmc.model # type: ignore
import matplotlib.pyplot as plt
import numpy as np
import csv

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

borcar = openmc.Material(name= "Boron Carbide")
borcar.add_element("B", 4)
borcar.add_element("C", 1)
borcar.set_density("g/cm3", 2.5)

stainless_steel = openmc.Material(name= "Stainless Steel")
stainless_steel.add_element("Fe", 0.8686)
stainless_steel.add_element("C", 0.0214)
stainless_steel.add_element("Cr", 0.1100)
stainless_steel.set_density("g/cm3", 8.05)

gadolinia = openmc.Material(name = "Gadolinium Oxide Gd2O3 Fuel")           # By mistake I run the whole depletion model only on Gd152 and
gadolinia.add_nuclide("Gd152", 0.02*2)                                      # optimized k-effective for this case. I have afterwards updated
gadolinia.add_nuclide("Gd154", 0.0218*2)                                    # the gadolinia material with the right Gd isotopes and get a 
gadolinia.add_nuclide("Gd155", 0.148*2)                                     # k-effective a factor of ca 0.1 lower. I want to cry.
gadolinia.add_nuclide("Gd156", 0.2047*2)
gadolinia.add_nuclide("Gd157", 0.1565*2)
gadolinia.add_nuclide("Gd158", 0.2484*2)
gadolinia.add_nuclide("Gd160", 0.2189*2)
gadolinia.add_element("O", 3)  
gadolinia.set_density("g/cm3", 7.07)
gadolinia.temperature = 1474

light_water = openmc.Material(name = "Water")                                                                                                                               # Defining Water material
light_water.add_element("H", 2)
light_water.add_element("O", 1)
light_water.set_density("g/cm3", 1)
light_water.temperature = 902

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

air = openmc.Material(name= "air")
air.add_element("N", 0.79, "ao")
air.add_element("O", 0.21, "ao")
air.set_density("g/cm3", 0.001293)

mox = openmc.Material.mix_materials([urox, gadolinia], [0.9, 0.1], 'wo')
mox.add_s_alpha_beta("c_U_in_UO2")

# Saving material temperatures to data file
# temperatures = [urox.temperature, gadolinia.temperature, light_water.temperature, boric_acid.temperature]
# csvfile = open('temperatures.csv', 'a') 
# csvwriter = csv.writer(csvfile)
# csvwriter.writerow(temperatures)
# csvfile.close()

##############################################################################################################################################################################################################################################################################################################################################################################################
# GEOMETRY
##############################################################################################################################################################################################################################################################################################################################################################################################

# Fuel Rod Diameter
frd = 0.95                                                                                                                                                                  # All of the reactor geometry dimentions will be decided as a function of the fuel rod diamater

# Define Boundary Type
boundary = "transmission"

# Top and Bottom Planes
z_min = openmc.ZPlane(-1.83, boundary_type = boundary)                                                                                                                         # Defining a minimum amd a maximum in the z-plane which will determine the height of the cylinders and applying boundary conditions: transmissive, vacuum or reflektive (for neutrons)
z_max = openmc.ZPlane(1.83, boundary_type = boundary)

# Fuel Rods
fuel_pin = -openmc.ZCylinder(r = frd/2 - 0.05 )                                                                                                                                       # Defining cylinders with infinte z axis and r (always given in cm)
helium_gap = -openmc.ZCylinder(r = frd/2 - 0.03)                                                                                                                                     # - sign defines interior of cylinder, + sign defines exterior. Helium_gap is a thin cylindrical shell 2mm thick, so r = 0.40 + 0.02
cladding_region = -openmc.ZCylinder(r = frd/2)                                                                                                                                 # r = 0.40 + 0.02 + 0.03

fuel_pin = fuel_pin & +z_min & -z_max                                                                                                                                       # Defining the cylinders' height / + sign means upwards while - sign means downwards
helium_gap = helium_gap & +z_min & -z_max
cladding_region = cladding_region & +z_min & -z_max

helium_gap = helium_gap & ~fuel_pin                                                                                                                                         # Making the cylindres hollow, ~ sign removes the following volume
cladding_region = cladding_region & ~helium_gap & ~fuel_pin

# Burnable Fuel Rods
bfuel_pin = -openmc.ZCylinder(r = frd/2 - 0.05 )                                                                                                                                       # Defining cylinders with infinte z axis and r (always given in cm)
bhelium_gap = -openmc.ZCylinder(r = frd/2 - 0.03)                                                                                                                                     # - sign defines interior of cylinder, + sign defines exterior. Helium_gap is a thin cylindrical shell 2mm thick, so r = 0.40 + 0.02
bcladding_region = -openmc.ZCylinder(r = frd/2)                                                                                                                                 # r = 0.40 + 0.02 + 0.03

bfuel_pin = bfuel_pin & +z_min & -z_max                                                                                                                                       # Defining the cylinders' height / + sign means upwards while - sign means downwards
bhelium_gap = bhelium_gap & +z_min & -z_max
bcladding_region = bcladding_region & +z_min & -z_max

bhelium_gap = bhelium_gap & ~bfuel_pin                                                                                                                                         # Making the cylindres hollow, ~ sign removes the following volume
bcladding_region = bcladding_region & ~bhelium_gap & ~bfuel_pin

# Control Rods
control_pin = -openmc.ZCylinder(r = frd/2 - 0.05)
gass_gap = -openmc.ZCylinder(r = frd/2 - 0.03) 
control_cladding = -openmc.ZCylinder(r = frd/2) 
control_thimble = -openmc.ZCylinder(r = frd/2 + 0.03) 

control_pin = control_pin & +z_min & -z_max 
gass_gap = gass_gap & +z_min & -z_max
control_cladding = control_cladding & +z_min & -z_max
control_thimble = control_thimble & +z_min & -z_max

gass_gap = gass_gap & ~control_pin
control_cladding = control_cladding & ~gass_gap & ~control_pin
control_thimble = control_thimble & ~control_cladding & ~gass_gap & ~control_pin

# Instrumental/Control Guide Thimbles
air_gap = -openmc.ZCylinder(r = frd/2) 
thimble_cladding = -openmc.ZCylinder(r = frd/2 + 0.03)  

air_gap = air_gap & +z_min & -z_max
thimble_cladding = thimble_cladding & +z_min & -z_max

thimble_cladding = thimble_cladding & ~air_gap

# Outer box enviroment
box = -openmc.model.RectangularPrism(width = 2*frd, height = 2*frd, boundary_type = boundary)                                                                                   # Defining an external box
box = box & +z_min & -z_max

fuel_box = box & ~cladding_region & ~helium_gap & ~fuel_pin                                                                                                                 # None of the geometry elements is overlapping
bfuel_box = box & ~bcladding_region & ~bhelium_gap & ~bfuel_pin 
control_box = box & ~control_thimble &~control_cladding & ~gass_gap & ~control_pin
thimble_box = box & ~thimble_cladding & ~air_gap
water_box = box     

##############################################################################################################################################################################################################################################################################################################################################################################################
# MATERIALS VOLUMES
##############################################################################################################################################################################################################################################################################################################################################################################################

urox.volume = (17*17-8*4-69)*(17*17-28-25)*(np.pi*(frd/2 - 0.05)**2)*366
mox.volume = (17*17-8*4-69)*(28)*(np.pi*(frd/2 - 0.05)**2)*366
helium.volume = (17*17-8*4)*(17*17-1)*(np.pi*((frd/2 - 0.03)**2 - (frd/2 - 0.05)**2))*366
zircaloy4.volume = (17*17-8*4)*(17*17)*(np.pi*((frd/2)**2 - (frd/2 - 0.03)**2))*366
borcar.volume = ((69)*(17*17)*(np.pi*(frd/2 - 0.05)**2) + (17*17-8*4-69)*(25)*(np.pi*(frd/2 - 0.05)**2))*366*0.255
stainless_steel.volume = ((69)*(17*17)*(np.pi*((frd/2 + 0.03)**2 - (frd/2)**2)) + (17*17-8*4-69)*(25)*(np.pi*((frd/2 + 0.03)**2 - (frd/2)**2)))*366
borated_water.volume = (17*17-8*4)*(21**2)*366 - (urox.volume - mox.volume - helium.volume - zircaloy4.volume - borcar.volume - stainless_steel.volume)

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

# Burnable Fuel Rods Cells
bfuel_cell = openmc.Cell(fill = mox, region = bfuel_pin)                                                                                                                     # Filling the geometry elements with fuel material
bhelium_cell = openmc.Cell(fill = helium, region = bhelium_gap)
bcladding_cell = openmc.Cell(fill = zircaloy4 , region = bcladding_region)
bwater_cell = openmc.Cell(fill = borated_water, region = bfuel_box)

cell_list_fuel = [bfuel_cell, bhelium_cell, bcladding_cell, bwater_cell] 

bfuel_universe = openmc.Universe(cells = cell_list_fuel) 

# Control Rods Cells
control_cell = openmc.Cell(fill = borcar, region = control_pin)                                                                                                             # Filling the geometry elements with fuel material
gass_cell = openmc.Cell(fill = helium, region = gass_gap)
control_cladding_cell = openmc.Cell(fill = zircaloy4 , region = control_cladding)
guide_cell = openmc.Cell(fill = stainless_steel, region = control_thimble)
moderator_cell = openmc.Cell(fill = borated_water, region = control_box)

cell_list_control = [control_cell, gass_cell, control_cladding_cell, guide_cell, moderator_cell] 

controlpin_universe = openmc.Universe(cells = cell_list_control)  

# Instrumental/Guide Thimbles Cells                                                                                                               
air_cell = openmc.Cell(fill = air, region = air_gap)
thimble_cell = openmc.Cell(fill = stainless_steel, region = thimble_cladding)
water_cell = openmc.Cell(fill = borated_water, region = thimble_box)

cell_list_thimble = [air_cell, thimble_cell, water_cell] 

thimble_universe = openmc.Universe(cells = cell_list_thimble) 

# Water Box (Box without Control Rods) Cells
control_water_cell = openmc.Cell(fill = borated_water, region = water_box)
water_universe = openmc.Universe(cells = [control_water_cell])

##############################################################################################################################################################################################################################################################################################################################################################################################
# FUEL ASSEMBLY
##############################################################################################################################################################################################################################################################################################################################################################################################

f = fuel_universe
b = bfuel_universe
t = thimble_universe
c = controlpin_universe
w = water_universe

# Gap Between Pincells
gap = 0.255

# Pin Size
pinx = frd + gap
piny = pinx
pinz = 3.66

# Assembly Size
assembly_x = 17
assembly_y = 17
assembly_z = 100

# Making an empty universe for the fuel/control rods assembly
universii = np.zeros((assembly_z,assembly_x,assembly_y), dtype = openmc.Universe)

# Control Rods Insertion Percentage (000/100)
crip = 25

# Making array for universes
for z in range(0,assembly_z):
    for i in range (0,assembly_x):
        for j in range(0,assembly_y):
            if i in [3,6,10,13] and j in [2,8,14]:                                                                                                                               # Decide placement of controlpins
                universii[z][i][j] = b
            elif i in [2,8,14] and j in [3,6,10,13]:
                universii[z][i][j] = b
            elif i in [4,12] and j in [4,12]:
                universii[z][i][j] = b
            elif i in [5,8,11] and j in [2,5,8,11,14] and z < 100-crip:                                                                                                                              # Decide placement of controlpins
                universii[z][i][j] = t
            elif i in [2,5,8,11,14] and j in [5,8,11] and z < 100-crip:                                                                                                                              # Decide placement of controlpins
                universii[z][i][j] = t
            elif i in [3,13] and j in [3,13] and z < 100-crip:                                                                                                                              # Decide placement of controlpins
                universii[z][i][j] = t
            elif i in [5,8,11] and j in [2,5,8,11,14] and z >= 100-crip:                                                                                                                               # Decide placement of controlpins
                universii[z][i][j] = c
            elif i in [2,5,8,11,14] and j in [5,8,11] and z >= 100-crip:                                                                                                                               # Decide placement of controlpins
                universii[z][i][j] = c
            elif i in [3,13] and j in [3,13] and z >= 100-crip:                                                                                                                               # Decide placement of controlpins
                universii[z][i][j] = c
            elif i == 8 and j == 8:
                universii[z][i][j] = t                                                                                                                                      # Place controlpins
            else:
                universii[z][i][j] = f                                                                                                                                      # Place fuel pins

# Making the Fuel Assembly Lattice
f_assembly = openmc.RectLattice(name = "fuel_assembly")
f_assembly.lower_left = (-assembly_x*abs(pinx)/2,-assembly_y*abs(piny)/2,-assembly_z*abs(pinz)/2)
f_assembly.universes = universii
f_assembly.pitch = (pinx,piny,pinz)

# Making Outer Universe
outer_water_cell = openmc.Cell(fill = borated_water)
outer_uni = openmc.Universe(cells = [outer_water_cell])
f_assembly.outer = outer_uni

# Putting the Lattice into the Assembly Universe
outedge = "transmission"                                                                                                                                                      # Infinite reactor or singular assembly
outcap= openmc.ZPlane(pinz*assembly_z/2, boundary_type = "vacuum")
outbot = openmc.ZPlane(-pinz*assembly_z/2, boundary_type = "vacuum")
outcell = openmc.model.RectangularPrism(width = pinx*assembly_x, height = piny*assembly_y, boundary_type = outedge)        
f_cell = openmc.Cell(fill = f_assembly, region = -outcell &+outbot &-outcap)                                                                                                # Making cell for assembly
cell_list = [f_cell]
f_assembly_uni = openmc.Universe(cells = cell_list)                                                                                                                           # Putting the assembly into a universe

##############################################################################################################################################################################################################################################################################################################################################################################################
# CONTROL ASSEMBLY
##############################################################################################################################################################################################################################################################################################################################################################################################

f = fuel_universe
b = bfuel_universe
t = thimble_universe
c = controlpin_universe
w = water_universe

# Gap Between Pincells
gap = 0.255

# Pin Size
pinx = frd + gap
piny = pinx
pinz = 3.66

# Assembly Size
assembly_x = 17
assembly_y = 17
assembly_z = 100

# Making an empty universe for the fuel/control rods assembly
c_universii = np.zeros((assembly_z,assembly_x,assembly_y), dtype = openmc.Universe)

# Control Rods Insertion Percentage (000/100)
crip = 25

# Making array for universes
for z in range(0,assembly_z):
    for i in range (0,assembly_x):
        for j in range(0,assembly_y):
            if i == 8 and j == 8:
                c_universii[z][i][j] = t
            elif z < 100-crip:
                c_universii[z][i][j] = t                                                                                                                                  # Place controlpins
            else:
                c_universii[z][i][j] = c                                                                                                                                     # Place fuel pins

# Making the Control Assembly Lattice
c_assembly = openmc.RectLattice(name = "control_assembly")
c_assembly.lower_left = (-assembly_x*abs(pinx)/2,-assembly_y*abs(piny)/2,-assembly_z*abs(pinz)/2)
c_assembly.universes = c_universii
c_assembly.pitch = (pinx,piny,pinz)

# Making Outer Universe
outer_water_cell = openmc.Cell(fill = borated_water)
outer_uni = openmc.Universe(cells = [outer_water_cell])
c_assembly.outer = outer_uni

# Putting the Lattice into a Control Assembly Universe
outedge = "transmission"                                                                                                                                                      # Infinite reactor or singular assembly
outcap= openmc.ZPlane(pinz*assembly_z/2, boundary_type = "vacuum")
outbot = openmc.ZPlane(-pinz*assembly_z/2, boundary_type = "vacuum")
outcell = openmc.model.RectangularPrism(width = pinx*assembly_x, height = piny*assembly_y, boundary_type = outedge)        
control_cell = openmc.Cell(fill = c_assembly, region = -outcell &+outbot &-outcap)                                                                                                # Making cell for assembly
cell_list = [control_cell]
c_assembly_uni = openmc.Universe(cells = cell_list)                                                                                                                           # Putting the assembly into a universe

##############################################################################################################################################################################################################################################################################################################################################################################################
# REACTOR CORE
##############################################################################################################################################################################################################################################################################################################################################################################################

# Core Size
core_x = 17
core_y = 17

# Making an empty universe for the core assembly
universiii = np.zeros((core_x,core_y), dtype = openmc.Universe)

# Making array for universes
for i in range(0,core_x):
    for j in range (0,core_y): 
        if i in [4,6,8,10,12] and j in [4,6,8,10,12]:
            universiii[i][j] = c_assembly_uni
        elif i in [3,5,7,9,11,13] and j in [1,15]:
            universiii[i][j] = c_assembly_uni
        elif i in [1,15] and j in [3,5,7,9,11,13]:
            universiii[i][j] = c_assembly_uni
        elif i in [2,6,8,10,14] and j in [2,14]:
            universiii[i][j] = c_assembly_uni
        elif i in [2,14] and j in [6,8,10]:
            universiii[i][j] = c_assembly_uni
        elif i in [3,13] and j in [3,13]:
            universiii[i][j] = c_assembly_uni
        else:
            universiii[i][j] = f_assembly_uni

# Making Outer Neutron Reflector Universe
outer_reflector_cell = openmc.Cell(fill = graphite)
outer_reflector_uni = openmc.Universe(cells = [outer_reflector_cell])

# Filling upper left corner with neutron reflector
x=0               
for i in range(0,4):
    for j in range (0,4-x):
        if j == 2 and i == 1:
            universiii[i][j] = f_assembly_uni
        
        elif j == 1 and i == 2:
            universiii[i][j] = f_assembly_uni
        
        else:
            universiii[i][j] = outer_reflector_uni
    x += 1

# Filling lower left corner with neutron reflector
x=3
for i in range(core_x-4,core_x):
    for j in range (0,4-x): 
        if j == 2 and i == core_x-2:
            universiii[i][j] = f_assembly_uni
        
        elif j == 1 and i == core_x-3:
            universiii[i][j] = f_assembly_uni
        
        else:
            universiii[i][j] = outer_reflector_uni                               
    x -= 1
        
# Filling upper right corner with neutron reflector
x=4
for i in range(0,4):
    for j in range (core_x-x,core_y):
        if j == core_y-3 and i == 1:
            universiii[i][j] = f_assembly_uni
        
        elif j == core_y-2 and i == 2:
            universiii[i][j] = f_assembly_uni
        
        else:
            universiii[i][j] = outer_reflector_uni                              
    x -= 1
        
# Filling lower right corner with neutron reflector
x=0
for i in range(core_x-(5-x),core_x):
    for j in range (core_y-x,core_y):
        if j == core_y-3 and i == core_x-2:
            universiii[i][j] = f_assembly_uni
        
        elif j == core_y-2 and i == core_x-3:
            universiii[i][j] = f_assembly_uni
        
        else:
            universiii[i][j] = outer_reflector_uni
    x += 1

# Making the Core Lattice
c_assembly = openmc.RectLattice(name = "core_assembly")
c_assembly.lower_left = (-core_x*assembly_x*abs(pinx)/2,-core_y*assembly_y*abs(piny)/2)
c_assembly.universes = universiii
c_assembly.pitch = (pinx*assembly_x,piny*assembly_y)

# Making Outer Neutron Reflector Universe
c_outer_reflector_cell = openmc.Cell(fill = graphite)
c_outer_uni = openmc.Universe(cells = [c_outer_reflector_cell])
c_assembly.outer = c_outer_uni

# Putting the Lattice into the Core Universe
c_outedge = "transmission"                                                                                                                                                    # Infinite reactor or singular assembly
c_outcap= openmc.ZPlane(pinz*assembly_z/2, boundary_type = "vacuum")
c_outbot = openmc.ZPlane(-pinz*assembly_z/2, boundary_type = "vacuum")
c_outcell = openmc.model.RectangularPrism(width = pinx*assembly_x*core_x*2, height = piny*assembly_y*core_y*2, boundary_type = outedge)        
c_cell = openmc.Cell(fill = c_assembly, region = -c_outcell &+c_outbot &-c_outcap)                                                              
cell_list = [c_cell]
core_uni = openmc.Universe(cells = cell_list) 

##############################################################################################################################################################################################################################################################################################################################################################################################
# REACTOR PRESSURE VESSEL
##############################################################################################################################################################################################################################################################################################################################################################################################

# Neutron Reflector
neutron_reflector = -openmc.ZCylinder(r = (pinx*assembly_x*(core_x+4))/2 , boundary_type = outedge)        
neutron_reflector_cell = openmc.Cell(fill = core_uni, region = neutron_reflector & +c_outbot & -c_outcap)                                                
graphite.volume = (np.pi * ((pinx*assembly_x*(core_x+4))/2)**2 - ((core_x*core_y) - (4*8))*(pinx*assembly_x)**2)*366

# Tank
r_t = 10 # Tank thickness
tank= -openmc.ZCylinder(r = (pinx*assembly_x*(core_x+4))/2 + r_t , boundary_type = c_outedge)
tank_thickness = tank & ~neutron_reflector & +c_outbot & -c_outcap
tank_cell = openmc.Cell(fill = stainless_steel, region = tank_thickness)

# Tank Air Gap
r_g = 30  
tank_gap = -openmc.ZCylinder(r = (pinx*assembly_x*(core_x+4))/2 + r_g, boundary_type = c_outedge)
tank_gap = tank_gap & ~tank_thickness & +c_outbot & -c_outcap
tank_gap_cell = openmc.Cell(fill = air, region = tank_gap)

# Reactor Pressure Vessel
r_rpv = 50  
rpv = -openmc.ZCylinder(r = (pinx*assembly_x*(core_x+4))/2 + r_rpv, boundary_type = "vacuum")
rpv = rpv & ~tank_gap & +c_outbot & -c_outcap
rpv_cell = openmc.Cell(fill = stainless_steel, region = rpv)

# Complete Reactor
cell_list = [neutron_reflector_cell, tank_cell, tank_gap_cell, rpv_cell]
rpv_uni = openmc.Universe(cells = cell_list)    

##############################################################################################################################################################################################################################################################################################################################################################################################
# SIMULATION SETTINGS 
##############################################################################################################################################################################################################################################################################################################################################################################################

materials = openmc.Materials([urox, helium, zircaloy4, light_water, mox, borcar, stainless_steel, graphite, air, boric_acid, borated_water])
materials.export_to_xml()                                                                                                                                                   # OpenMc does not run a python script so everything is saved into xml files. (The smallest set of xml files neede is: materials.xml, geometry.xml og settings.xml)

geometry = openmc.Geometry()
geometry.root_universe = rpv_uni
geometry.export_to_xml()

settings = openmc.Settings()
settings.batches = 100                                                                                                                                                      # How many iterations/batches of the same neutron generation to simulate
settings.inactive = 10                                                                                                                                                      # Number of inactive batches = How many iterations are used to determine sample spots (MC simulation)
settings.particles = 10000                                                                                                                                                  # Number of neutrons in each generation
settings.temperature = {"method":"interpolation", "tolerance":100}
#settings.generations_per_batch = (int)                                                                                                                                     # Number of generations per batch (not used here)

##############################################################################################################################################################################################################################################################################################################################################################################################
# OpenMc PARTICLE SAMPLE SITES 
##############################################################################################################################################################################################################################################################################################################################################################################################

bounds = [-(pinx*assembly_x*core_x+4)/2 + r_rpv/2,-(piny*assembly_y*core_y+4)/2+ r_rpv/2,-(pinz*assembly_z)/2,(pinx*assembly_x*core_x+4)/2+ r_rpv/2,(piny*assembly_y*core_y+4)/2+ r_rpv/2,(pinz*assembly_z)/2]           # Coordinates of lower left corner (first 3 entries) and of upper right corner (last 3 entries). Bound is 1.4 cm wide and 2 cm high
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:])                                                                                         
settings.source = openmc.IndependentSource(space=uniform_dist, constraints={'fissionable': True})                                                                           # Generating neutrons only in fissible material (e.g. not in water but in urox)            
settings.export_to_xml()

##############################################################################################################################################################################################################################################################################################################################################################################################
# PLOTS AND PICTURES
##############################################################################################################################################################################################################################################################################################################################################################################################
'''
fuel_cell.plot(color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})                                                                                        
plt.savefig("Pictures/Fuelpin")

helium_cell.plot(color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Heliumgap")

fuel_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Fuelrods_plot_xy")

fuel_universe.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Fuelrods_plot_yz")

bfuel_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Burnablefuelrods_plot_xy")

bfuel_universe.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Burnablefuelrods_plot_yz")

controlpin_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Controlrods_plot_xy")

controlpin_universe.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Controlrods_plot_yz")

thimble_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Thimbleguides_plot_xy")

thimble_universe.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Thimbleguides_plot_yz")

water_universe.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Water_Universe_plot_xy")

f_assembly_uni.plot(basis = "xy", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/f_Assembly_Universe_plot_xy")

f_assembly_uni.plot(basis = "yz", pixels = [800,800], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/f_Assembly_Universe_plot_yz")

core_uni.plot(basis = "xy", pixels = [1600,1600], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/Core_Universe_plot_yz")

rpv_uni.plot(basis = "xy", pixels = [5000,5000], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"})
plt.savefig("Pictures/RPV_Universe_plot_xy")

rpv_uni.plot(basis = "yz", pixels = [5000,5000], color_by = "material",
                  colors = {urox: "lime", helium: "gold", air: "white", mox: "fuchsia", borated_water: "blue",
                  zircaloy4: "lightgray", light_water : "blue", borcar : "red", stainless_steel : "gray", graphite : "midnightblue"}, origin=(0,0,0))
plt.savefig("Pictures/RPV_Universe_plot_yz")

'''
##############################################################################################################################################################################################################################################################################################################################################################################################
# TALLY ARITHMETICS 
##############################################################################################################################################################################################################################################################################################################################################################################################

# Defining Tallies File
tallies_file = openmc.Tallies()                                                                                                                                             # For our tallies to be included we have to put them into an inp

# Defining a Spatial Mesh for Fuel Pin Cell
pin_mesh = openmc.RegularMesh()                                                                                                                                                 # Making a mesh instance
pin_mesh.dimension = [100,100]                                                                                                                                                  # Setting the bins of the mesh to 100 x 100 bins
pin_mesh.lower_left = [-18*(pinx) - pinx/2,- piny/2]                                                                                                
pin_mesh.upper_right = [-18*(pinx) + pinx/2, piny/2]                                                                                                  
pin_mesh_filter = openmc.MeshFilter(pin_mesh)

# Spatial Distribution of the Prompt Fission Sites and the Spatial Neutron Flux Distribution for Fuel Pin Cell
pin_tally = openmc.Tally(name = "Pin Tally")                                                                                                                                        # Making a tally instance
pin_tally.filters = [pin_mesh_filter]                                                                                                                                               # Designatingen the mesh for the tally
pin_tally.scores = ["flux", 'prompt-nu-fission']                                                                                                                                # Designating the objects to be tallied
tallies_file.append(pin_tally)  

# Defining a Spatial Mesh for Reactor Pressure Vessel
rpv_mesh = openmc.RegularMesh()                                                                                                                                                 # Making a mesh instance
rpv_mesh.dimension = [100,100]                                                                                                                                                  # Setting the bins of the mesh to 100 x 100 bins
rpv_mesh.lower_left = [-(pinx*assembly_x*core_x+100)/2,-(piny*assembly_y*core_y+100)/2]                                                                                                 # Setting lower left coordinates
rpv_mesh.upper_right = [(pinx*assembly_x*core_x+100)/2,(piny*assembly_y*core_y+100)/2]                                                                                                  # Setting upper right coordinates
mesh_filter = openmc.MeshFilter(rpv_mesh)                                                                                                                                       # Making a meshfilter

# Spatial Distribution of the Prompt Fission Sites and the Spatial Neutron Flux Distribution for Reactor Pressure Vessel
tally = openmc.Tally(name = "Tally")                                                                                                                                        # Making a tally instance
tally.filters = [mesh_filter]                                                                                                                                               # Designatingen the mesh for the tally
tally.scores = ["flux", 'prompt-nu-fission']                                                                                                                                # Designating the objects to be tallied
tallies_file.append(tally)                                                                                                                                                  # Appending the tally to our tallies file

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
mesh.lower_left = [-(pinx*assembly_x*core_x+100)/2,-(piny*assembly_y*core_y+100)/2, -366/2]
mesh.upper_right= [(pinx*assembly_x*core_x+100)/2,(piny*assembly_y*core_y+100)/2, 366/2]
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
# DEPLETION MODEL
##############################################################################################################################################################################################################################################################################################################################################################################################
'''
reactor_model = openmc.Model()
reactor_model.geometry = geometry
reactor_model.materials = materials
reactor_model.settings = settings
reactor_model.tallies = tallies_file

# Designating the depletion chain being used for our depletion calculation
chain = openmc.deplete.Chain.from_xml("./chain_endfb80_pwr.xml")

operator = openmc.deplete.CoupledOperator(reactor_model, "./chain_endfb80_pwr.xml")

# Designating power settings in watts
power = 1200*10**6 

# Number of timesteps and their length: in this case 56 steps of one week length in seconds, resulting in a 13 month cycle.
time_steps = [7* 24 * 60 * 60] * 56

#Running the depletion calculation
integrator = openmc.deplete.PredictorIntegrator(operator, time_steps, power)
integrator.integrate()
'''

##############################################################################################################################################################################################################################################################################################################################################################################################
# RUNNING OpenMC 
##############################################################################################################################################################################################################################################################################################################################################################################################

if __name__ == "__main__":
    openmc.run()  #threads = 4                                                                                                                                               # Threads set for heplab primarily (to limit utilized cores on heplab server, this can be removed when run locally)

##############################################################################################################################################################################################################################################################################################################################################################################################
# END
##############################################################################################################################################################################################################################################################################################################################################################################################

exec(open("Post_processing.py").read())