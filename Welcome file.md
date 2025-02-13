# exercise 01 :: hop on
Aims of the exercise is to introduce linux software platform used for the exercises, get familiar with the linux command line environment, basics of bash scripting and python programming languge, and atomistic simulation environment (ASE).
- for linux commands line [see cheat sheets](https://webdisk.science.upjs.sk/~martin_gmitra/Atomistic%20Computer%20Modeling%20of%20Materials/linux%20cheat%20sheets)
- see [what is Python?](https://wiki.fysik.dtu.dk/ase/python.html) for more details
- more about ASE see [homepage](https://wiki.fysik.dtu.dk/ase/about.html)


## TASK 01

1) open command line terminal and create folder ASE, change working directory to the ASE folder, create and edit a file ASE_ase.py within the aluminium folder
>    mkdir ASE
>  cd ASE
>  nano ASE_ase.py

Save the file content with the "Hello world" program should look like:
```sh
#!/usr/bin/python
print('Hello, world!')
```
2) Let's explore the basic elements of the ASE package. The most fundamnetal structural unit is ```Atom``` object (see [Atom](https://wiki.fysik.dtu.dk/ase/ase/atom.html) for more details). Try to create the aluminium atom at position r=(0,0,0), with aluminium atomic mass and charge (see [ptable](https://ptable.com/#Properties)) using 	```python_apmm```.
>from ase import Atom
>al = Atom('Al', position = (0, 0, 0), charge = ?, mass = ?)

Inspect the properties of the ```al``` object with 
>al.mass (charge/momentum/...)

In order to crate more interesting systems consisting of several atoms ASE provides us with ```Atoms``` object. Check the parameters of this object at  [Atoms](https://wiki.fysik.dtu.dk/ase/ase/atoms.html). The Atoms object can represent an isolated molecule, or a periodically repeated structure. It has a unit cell and there may be periodic boundary conditions along any of the three unit cell axes. Information about the atoms (atomic numbers and position) is stored in ndarrays. Optionally, there can be information about tags, momenta, masses, magnetic moments and charges.

Create the simple ```N3``` molecule:
>from  ase  import Atoms
>atoms = Atoms("N3", [(0, 0, 0), (1, 0, 0), (0, 0, 1)])

Let's now check the positions on which were the individual atoms set and move them to new ones and print them out
>atoms.get_positions()
>atoms.set_positions([(2, 0, 0), (0, 2, 2), (2, 2, 0)])
>atoms.get_positions()

Atoms object can be assigned a ```cell``` determined by 3 cell parameters and angles between them or by an array of cell vectors. Set the atoms a cubic cell with lattice parameters ```2 ang```.
>atoms.set_cell(2*np.identity(3))
>atoms.get_cell()

Alternative way is to create a cell object with the same properties as
>from ase.cell import Cell
>cellpar = [2, 2, 2, 90, 90, 90]         #[a, b, c, alpha, beta, gamma]
>cell = Cell.fromcellpar(cellpar) 
>atoms.set_cell(cell)

Visualization of the Atoms object can be performed using [view](https://wiki.fysik.dtu.dk/ase/ase/visualize/visualize.html) function
>from ase.visualize import view
>view(atoms) 

<center>
  <img src="https://github.com/timmos1/APMM_2025/blob/master/N3.png?raw=true" width="200" height="200">
</center>

If there is a need to change the original cell (squeeze or stretch, bend it) together with the atomic structure inside, we just set the new cell with 
>atoms.set_cell(3*cell,scale_atoms=True)
>atoms.get_positions()
>view(atoms)





## TASK 02

1)Create a new directory silllicon 
> cd ..
> mkdir silicon
> cd silicon

2) create the file Po_ase.py of the following content
```sh
#!/usr/bin/python3

import ase
import ase.io as io
from ase import Atoms, Atom
from ase.lattice.cubic import SimpleCubic
from ase.visualize import view

polonium = SimpleCubic(directions=[[1,0,0], [0,1,0], [0,0,1]], symbol='Po', latticeconstant = 3.359)
view(polonium)

#write output for crystal visualization (use e.g. VESTA)
struct = polonium
ase.io.write("Po_ase.cif", struct, "cif")
```
 Run the script
> python3 Po_ase.py

which uses ASE to generate a [simple cubic crystal](https://wiki.fysik.dtu.dk/ase/ase/lattice.html) of polonium,  and writes cif structure file for further visualization, e.g., using VESTA program
> VESTA Po_ase.cif


 
3) create slab extending in <111> direction made of GaAs in zincblende structure

```sh
#!/usr/bin/python3

import numpy as np
import ase.io as io
from ase.spacegroup import crystal
from ase.visualize import view

a = 5.6533
X = 'Ga'
Y = 'As'

# basis is defined in terms of the unit vectors
system = crystal([X, Y], basis = [(0.0, 0.0, 0.0), (0.25, 0.25,0.25)], spacegroup=216, cellpar=[a, a, a, 90, 90, 90], size = (1,1,4))

view(system)

```
<center>
  <img src="
https://github.com/timmos1/APMM_2025/blob/master/Screenshot%20from%202025-02-13%2019-28-40.png?raw=true" width=400>
</center>

4) create unit cell of InSb in wurtzite crystal structure
```sh
#!/usr/bin/python3

from ase import *
from math import *
from ase.build import bulk
from ase.visualize import view
import ase.io as io

a = 4.5712
c = 7.5221
ca = c/a
epsilon = -0.00097
ua=3/8.+epsilon

atoms = Atoms([Atom('Sb', (2/3., 1/3.,    0)),
               Atom('Sb', (1/3., 2/3., 1/2.)),
               Atom('In', (2/3., 1/3., ua)),
               Atom('In', (1/3., 2/3., 1/2.+ua))])
cell = [(a*sqrt(3)/2, -a/2, 0),
        (          0,    a, 0),
        (          0,    0, c)]
atoms.set_cell(cell, scale_atoms=True)

view(atoms)

io.write("InSb_ase.cif", atoms, "cif")
```

**Note:** for more python info see [reference card](https://perso.limsi.fr/pointal/_media/python:pqrc:pqrc-2.4-a4-latest.pdf) or [python cheat sheet](https://www.pythoncheatsheet.org/)


## TASK 03
For the InSb structure from previous task create an input into [Quantum Espresso](https://www.quantum-espresso.org/)  [pw.x](https://www.quantum-espresso.org/Doc/INPUT_PW.html) code using ASE ```write_espresso_in``` function. 

```sh
#!/usr/bin/python3

import numpy as np
import ase
from math import *
from ase.build import bulk
from ase.visualize import view
import ase.io as io
from ase.io.espresso import write_espresso_in

atoms = io.read("InSb_ase.cube")

view(atoms)

pseudopotentials={
    'In':'In_ONCV_PBE_sr.upf',
    'Sb':'Sb_ONCV_PBE'
        }

kpoints = np.ndarray([6,6,6])
input_data = {
    'calculation': 'relax',
    'restart_mode': 'from_scratch',
    'tprnfor': True,
    'etot_conv_thr': 1e-5,
    'forc_conv_thr': 1e-4,
    'ecutwfc': 60,
    'ecutrho': 480,
    'input_dft': 'rpbe',
    'vdw_corr': 'dft-d3',
    'occupations': 'smearing',
    'degauss': 0.01,
    'smearing': 'cold',
    'conv_thr': 1e-8,
    'mixing_mode': 'local-TF',
    'mixing_beta': 0.35,
    'diagonalization': 'david',
    'ion_dynamics': 'bfgs',
    'bfgs_ndim': 6,
    'startingwfc': 'random',
}  # This flat dictionary will be converted to a nested dictionary where, for example, "calculation" will be put into the "control" section


inp = open('pw.in','w')

ase.io.espresso.write_espresso_in(fd=inp, atoms = atoms, pseudopotentials = pseudopotentials,kpts=(4,4,4), input_data = input_data, format='espresso-in')

#write output for crystal visualization (use e.g. VESTA)
io.write("InSb.cif", atoms, "cif")
io.write("InSb.xsf", atoms, "xsf")
```

Visualize the outputs with ```VESTA``` and ```xcrysden``` programs.


## TASK 04
1) use ASE to generate structure of bulk Aluminium (space group #225, [FCC structure](http://lampx.tugraz.at/~hadley/ss1/crystalstructure/structures/fcc/fcc_jsmol.php))
2) use ASE to generate structure of bulk Silicon (space group #227, [diamond crystal structure](http://lampx.tugraz.at/~hadley/memm/materials/silicon/silicon.php))

**Note 1:** for lattice constants see also https://periodictable.com/Properties/A/LatticeConstants.html
**Note 2:** for more crystallographic data see [Bilbao Crystallographic Server](https://www.cryst.ehu.es/) 

---

#### Further links:
- [ASE-workshop 2019 @ Chalmers](https://ase-workshop.materialsmodeling.org/)
- [The Structure of Materials](http://som.web.cmu.edu/frames2.html) (book)
- [The Material Project](https://materialsproject.org/) web-based access to computed material properties
---
Martin Gmitra :: martin.gmitra@upjs.sk
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTc1MzUwNjA2NSwtMTg1NzE5MDY4MiwtMT
Y4ODQ2NzUyNSwxMTI0MzcxNTQ2LDE3ODA3Njg0NzRdfQ==
-->
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0OTUzNzkxMzQsMTY2Mzc1MjIxNiwtMT
A5MjMzMjI4OCwxMTA1MjAyMTQ1LDE5NzkyMDMwMjcsLTIxMjk4
MzYxNjVdfQ==
-->