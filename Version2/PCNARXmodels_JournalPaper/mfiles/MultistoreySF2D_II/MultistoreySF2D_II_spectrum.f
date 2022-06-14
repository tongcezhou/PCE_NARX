/CLEAR
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7

! Material 1: Vertical beams 
ET,1,BEAM23
R,1,0.04,1.33333E-04,0.20				! Real Const: area,I,height
R,2,0.030625,7.81576E-05,0.175			! Real Const: area,I,height
R,3,0.0225,4.21875E-05,0.15				! Real Const: area,I,height
R,4,0.015625,2.03451E-05,0.125			! Real Const: area,I,height
R,5,0.01,8.333333E-06,0.10			! Real Const: area,I,height
MP,EX,1,200E9
MP,PRXY,1,0.29
MP,DENS,1,7850
KEYOPT,1,2,0							! No shear deflection
KEYOPT,1,4,1							! Print out member forces and moments in the element coordinate system
KEYOPT,1,6,0							! No shear deflection

! Material 2: Horizontal beams
R,6,0.0225,4.21875E-05,0.15				! Real Const: area,I,height
MP,EX,2,200E11
MP,PRXY,2,0.29
MP,DENS,2,25972


! Keypoints (front profile)
! Ground floor
K,1, 0.0, 0.0
K,2, 4.0, 0.0
K,3, 8.0, 0.0
K,4, 12.0, 0.0
! 1st floor
K,5, 0.0, 3.0
K,6, 4.0, 3.0
K,7, 8.0, 3.0
K,8, 12.0, 3.0
! 2nd floor
K,9, 0.0, 6.0
K,10, 4.0, 6.0
K,11, 8.0, 6.0
K,12, 12.0, 6.0
! 3rd floor
K,13, 0.0, 9.0
K,14, 4.0, 9.0
K,15, 8.0, 9.0
K,16, 12.0, 9.0
! 4th floor
K,17, 0.0, 12.0
K,18, 4.0, 12.0
! 5th floor
K,19, 0.0, 15.0
K,20, 4.0, 15.0

! Lines
! Vertical Beams b/w grounf-1st floor
TYPE,1
L,1,5	! Element 1
L,2,6	! Element 2
L,3,7	! Element 3
L,4,8	! Element 4
! Vertical Beams b/w 1st floor-2nd floor
L,5,9	! Element 5
L,6,10	! Element 6
L,7,11	! Element 7
L,8,12	! Element 8
! Vertical Beams b/w 2nd floor-3rd floor
L,9,13	! Element 9
L,10,14	! Element 10
L,11,15	! Element 11
L,12,16	! Element 12
! Vertical Beams b/w 3rd floor-4th floor
L,13,17	! Element 13
L,14,18	! Element 14
! Vertical Beams b/w 4th floor-5th floor
L,17,19	! Element 15
L,18,20	! Element 16

! Horizontal Beams
L,5,6	! Element 17
L,6,7	! Element 18
L,7,8	! Element 19
L,9,10	! Element 20
L,11,12	! Element 21
L,13,14	! Element 22
L,15,16	! Element 23
L,17,18	! Element 24
L,19,20	! Element 25


! Meshing 
LESIZE,ALL,,,1,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
! Verical beams (ground-1st floor) 
LATT,1,1,1								! Mat, Real, Type
LSEL,S,LINE,,1,4,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (1st-2nd floor) 
LATT,1,2,1
LSEL,S,LINE,,5,8,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (2nd-3rd floor) 
LATT,1,3,1
LSEL,S,LINE,,9,12,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (3rd-4th floor) 
LATT,1,4,1
LSEL,S,LINE,,13,14,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (4th-5th floor) 
LATT,1,5,1
LSEL,S,LINE,,15,16,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Horizontal beams 
LATT,2,6,1
LSEL,S,LINE,,17,25,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL

FINISH

! =======================================================
! SPECTRUM SOLUTION 
! =======================================================
/SOLU

! Rayleigh Damping
! ================
ALPHAD,4.0
BETAD,0.0002
	
! Boundary conditions
! ===============
NSEL,S,LOC,Y,0
D,ALL,ALL
ALLSEL,ALL
	
ANTYPE,3			! Harmonic analysis
F,20,FX,5000
F,16,FX,5000

HARFRQ,0,40,   	! Frequency range
NSUBST,200, 		! Number of frequency steps
KBC,1				! Stepped loads

SOLVE
FINISH

/POST26
NSOL,2,2,U,X,UX2	! Get y-deflection data
PRVAR,2 			! Print data
PLVAR,2				! Plot data