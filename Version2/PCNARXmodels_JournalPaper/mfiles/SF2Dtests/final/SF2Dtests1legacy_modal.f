/CLEAR
/CWD 'D:\SpaceFrame\'
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
! Nonlinear material properties
TB,BISO,1,1								! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,200E6,10E9
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
! Verical beams (ground-1st floor) 
LATT,1,1,1								! Mat, Real, Type
LSEL,S,LINE,,1,2,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (1st-2nd floor) 
LATT,1,2,1
LSEL,S,LINE,,3,4,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (2nd-3rd floor) 
LATT,1,3,1
LSEL,S,LINE,,5,6,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (3rd-4th floor) 
LATT,1,4,1
LSEL,S,LINE,,7,8,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (4th-5th floor) 
LATT,1,5,1
LSEL,S,LINE,,9,10,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Horizontal inner beams 
LATT,2,6,1
LSEL,S,LINE,,11,15,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL


! List nodes
NLIST,ALL,,,,NODE
! List sections and their properties
SLIST,ALL,,,FULL
! List Elements and its properties
ELIST,ALL,,,0,1

FINISH

! =========================================================
! SOLUTION 
! =========================================================
/SOLU
ANTYPE,MODAL 					! Specifies the analysis type
MODOPT,LANB,20,0,100			! Specifies modal analysis options                                                                                                                                                        
EQSLV,SPAR          			! Specifies the type of equation solver                                                                                                                                                   
MXPAND,0,,,0      			! Number of modes to expand                                                                                                                                                   
LUMPM,0 						! Lumped mass matrix formulation (on/off)                                                                                                                                                            
PSTRES,0            			! Prestress effects included (on/off)                                                                                                                                                   
	
! Rayleigh Damping
! ================
ALPHAD,4.5
BETAD,0.0002
	
! Boundary conditions
! ===============
NSEL,S,LOC,Y,0
D,ALL,ALL
ALLSEL,ALL
	
SOLVE       		! Starts a solution                                                                                                                                                           
FINISH