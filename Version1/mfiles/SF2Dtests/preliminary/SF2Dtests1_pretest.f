/CLEAR
/CWD 'D:\SpaceFrame\'
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
*SET,NoS,1000
*SET,DT,0.05

! Material 1: Vertical beams 
ET,1,BEAM23
R,1,0.0225,4.21875E-05,0.15				! Real Const: area,I,height
MP,EX,1,200E9
MP,PRXY,1,0.29
MP,DENS,1,7850
! Nonlinear material properties
TB,BISO,1								! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,200E6,10E9
KEYOPT,1,2,0							! No shear deflection
KEYOPT,1,4,1							! Print out member forces and moments in the element coordinate system
KEYOPT,1,6,0							! No shear deflection

! Material 2: Horizontal beams
ET,2,BEAM23
R,2,0.0225,4.21875E-05,0.15				! Real Const: area,I,height
MP,EX,2,200E11
MP,PRXY,2,0.29
MP,DENS,2,33946
KEYOPT,2,2,0							! No shear deflection
KEYOPT,2,4,1							! Print out member forces and moments in the element coordinate system
KEYOPT,2,6,0							! No shear deflection

! Keypoints (front profile)
! Ground floor
K,1, 0.0, 0.0
K,2, 4.0, 0.0
! 1st floor
K,3, 0.0, 3.0
K,4, 4.0, 3.0
! 2nd floor
K,5, 0.0, 6.0
K,6, 4.0, 6.0
! 3rd floor
K,7, 0.0, 9.0
K,8, 4.0, 9.0
! 4th floor
K,9, 0.0, 12.0
K,10, 4.0, 12.0
! 5th floor
K,11, 0.0, 15.0
K,12, 4.0, 15.0

! Lines
! TYPE,1
! Vertical Beams b/w grounf-1st floor
TYPE,1
L,1,3	! Element 1
L,2,4	! Element 2
! Vertical Beams b/w 1st floor-2nd floor
L,3,5	! Element 3
L,4,6	! Element 4
! Vertical Beams b/w 2nd floor-3rd floor
L,5,7	! Element 5
L,6,8	! Element 6
! Vertical Beams b/w 3rd floor-4th floor
L,7,9	! Element 7
L,8,10	! Element 8
! Vertical Beams b/w 4th floor-5th floor
L,9,11	! Element 9
L,10,12	! Element 10

! TYPE,2
! Horizontal Beams 1st floor (outer beams)
L,3,4	! Element 11
L,5,6	! Element 12
L,7,8	! Element 13
L,9,10	! Element 14
L,11,12	! Element 15

! Meshing 
LESIZE,ALL,,,1,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
TYPE,1 									! Element type                                                                                                                                                                                                                                                                             
MAT,1 
REAL,1
LSEL,S,LINE,,1,10,1
LMESH,ALL
ALLSEL,ALL  

TYPE,2 									! Element type                                                                                                                                                                                                                                                                             
MAT,2 
REAL,2
LSEL,S,LINE,,11,15,1
LMESH,ALL								! Generates nodes and line elements along lines


! List nodes
! NLIST,ALL,,,,NODE
! List sections and their properties
! SLIST,ALL,,,FULL
! List Elements and its properties
! ELIST,ALL,,,0,1

FINISH

! =========================================================
! SOLUTION 
! =========================================================
/SOLU
SOLCONTROL,ON
! Rayleigh Damping
! ================
ALPHAD,4.265
BETAD,0.0006
	
! Transient Analysis Options
! ==========================
ANTYPE,STATIC			! Specifies the analysis type and restart status.

! Controls the solution data written to the database.                                                                                                                  
OUTRES,ALL,NONE 
OUTRES,ESOL,LAST
OUTRES,NSOL,LAST
OUTRES,NLOAD,LAST                                                                                                                                                       
OUTRES,RSOL,LAST      
OUTRES,A,LAST 
OUTRES,V,LAST 
	
! DOF constraints
! ===============
NSEL,S,LOC,Y,0
D,ALL,ALL				! Constraint all dofs for fixed supports
ALLSEL,ALL

F,2,FX,1000,,2,1		! 1st floor
F,5,FX,1000*2,,5,1	! 2nd floor
F,7,FX,1000*3,,7,1	! 3rd floor
F,9,FX,1000*4,,9,1	! 4th floor
F,11,FX,1000*5,,11,1	! 5th floor
NSEL,ALL

SOLVE
FINISH
