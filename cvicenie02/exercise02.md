
# exercise 02 :: Self-consistent field DFT calculations of Aluminium

Aliminium has 13 electrons. Its crystalline form is cubic face-centered with space group Fm-3m (number 225) with one atom in unit cell in Wyckhoff position 4a, see [Bilbao Crystallographic Server](https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-wp-list). Lattice constant is 4.0495 A.


## TASK 01

Construct structure of Al crystal in VESTA program (video tutorial is [here](https://webdisk.science.upjs.sk/~martin_gmitra/Atomistic%20Computer%20Modeling%20of%20Materials/Exercises/exercise%2002/building_Al_by_vesta.mp4)). Compare with `Al_cubic_SG225.vesta` on  [webdisk](https://webdisk.science.upjs.sk/~martin_gmitra/Atomistic%20Computer%20Modeling%20of%20Materials/Exercises/exercise%2002/). 


## TASK 02

Create the input file for the self-consistent field (SCF) calculations within Quantum Espresso pw.x program. For input variables refer to [Input File Description](https://www.quantum-espresso.org/Doc/INPUT_PW.html).


1. Save the following input file with the name `pw-scf_Al.in`, or download it from [webdisk](https://webdisk.science.upjs.sk/~martin_gmitra/Atomistic%20Computer%20Modeling%20of%20Materials/Exercises/exercise%2002/)

        &control
            calculation = 'scf'
            restart_mode = 'from_scratch',
            pseudo_dir = './',
            outdir = './tmp',
            prefix = 'al'
        /
        &system
            ibrav = 2, celldm(1) = 7.652446, nat = 1, ntyp = 1, 
            ecutwfc = 15.0, occupations = 'smearing',
            smearing = 'marzari-vanderbilt', degauss = 0.05
        /
        &electrons
            diagonalization = 'david'
            mixing_beta = 0.7
            conv_thr = 1.D-6
        /
        ATOMIC_SPECIES
            Al  26.98 Al.pz-vbc.UPF
        ATOMIC_POSITIONS
            Al 0.00 0.00 0.00
        K_POINTS (automatic)
            8 8 8  0 0 0

2. Visualize the  `pw-scf_Al.in` input file with the `XCrySDen` program, e.g., using command line `xcrysden --pwi pw-scf_Al.in`.
3. Download the pseudopotential file [`Al.pz-vbc.UPF`](https://webdisk.science.upjs.sk/~martin_gmitra/Atomistic%20Computer%20Modeling%20of%20Materials/Exercises/exercise%2002/Al.pz-vbc.UPF) with LDA approximation to exchange-correlation functional.

4. Run the SCF calculations `pw.x -i pw-scf_Al.in > pw-scf_Al.out`

5. Analyze the output file `pw-scf_Al.out` and seek for the total energy and contributions to the exchange-correlation, Hartree, and single electron energies. What is the Fermi energy?

6. How many valence electrons were used in the calculation?

---
Martin Gmitra :: martin.gmitra@upjs.sk
<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/80x15.png" /></a> This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.
<!--stackedit_data:
eyJoaXN0b3J5IjpbNDQ3MzU1MTc5LDExMjkzNjI3OTQsLTkwNz
UyODc5MCwxMjcxOTE4MjUwLDE0NzE3MjA1MjQsMTUxMjExODk4
NCwtNDM1OTc3OTg2LDExNzE5NjcyNjgsLTEyNjIxMTE0OTMsOD
QzNTI4NzkwLDE2NDQ0Mjc0MDksNjg3MDg2ODg5LDE4MDU2NTI5
MjYsLTk1ODI4MDUxNF19
-->

> Written with [StackEdit](https://stackedit.io/).
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTE0ODgwMjY0OTksNzMwOTk4MTE2XX0=
-->