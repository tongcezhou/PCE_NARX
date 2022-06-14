/CLEAR
/CWD 'D:\SpaceFrame\'
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
*SET,NoS,1000
*SET,DT,0.05

! Material 1
ET,1,BEAM188						! BEAM188 elements 
MP,EX,1,200E9						! Young modulus
MP,PRXY,1,0.29						! Poisson ratio
MP,DENS,1,7850						! Density
KEYOPT,1,3,3						! Cubic shape functions
KEYOPT,1,15,1						! Results given on section's integration points
! Nonlinear material properties
TB,BISO,1							! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,200E6,10E9					! Yield stress = 200 MPa, Region of plasticity strength = 10GPa

! Material 1: Beams 
ET,2,BEAM188
MP,EX,2,200E11
MP,PRXY,2,0.29
MP,DENS,2,7850
KEYOPT,2,3,3						! Cubic shape functions
KEYOPT,2,15,1						! Results given on section's integration points
! Nonlinear material properties
! TB,BISO,2							! Bilinear isotropic hardening                                                                                                                      
! TBDATA,1,200E6,10E9					! Yield stress = 200 MPa, Region of plasticity strength = 10GPa

! Cross-section 1 (vertical elements): Square cross section w = h = 0.15 m, 
SECTYPE,1,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 2 (vertical elements)
SECTYPE,2,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 3 (vertical elements)
SECTYPE,3,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 4 (vertical elements)
SECTYPE,4,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 5 (vertical elements)
SECTYPE,5,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 6 (horizontal outer elements)
SECTYPE,6,BEAM,RECT
SECDATA,0.15,0.15,3,3
SECCONTROLS,,,,407.75

! Keypoints (front profile)
! Ground floor
K,1, 0.0, 0.0, 0.0
K,2, 4.0, 0.0, 0.0
! 1st floor
K,3, 0.0, 0.0, 3.0
K,4, 4.0, 0.0, 3.0
! 2nd floor
K,5, 0.0, 0.0, 6.0
K,6, 4.0, 0.0, 6.0
! 3rd floor
K,7, 0.0, 0.0, 9.0
K,8, 4.0, 0.0, 9.0
! 4th floor
K,9, 0.0, 0.0, 12.0
K,10, 4.0, 0.0, 12.0
! 5th floor
K,11, 0.0, 0.0, 15.0
K,12, 4.0, 0.0, 15.0

! Lines
! Vertical Beams b/w grounf-1st floor
!TYPE,1
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

!TYPE,2
! Horizontal Beams 1st floor (outer beams)
L,3,4	! Element 11
L,5,6	! Element 12
L,7,8	! Element 13
L,9,10	! Element 14
L,11,12	! Element 15

! Meshing 
LESIZE,ALL,,,1,1,1  					! Specifies the divisions (3) and spacing ratio on unmeshed lines (no division, no spacing ratio)
! Verical beams (ground-1st floor) 
LATT,1,1,1,,,,1							! Associates element attributes with the selected, unmeshed lines
LSEL,S,LINE,,1,2,1,1					! Select ground floor columns
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  							! Select all
! Verical beams (1st-2nd floor) 
LATT,1,1,1,,,,2
LSEL,S,LINE,,3,4,1,1
LMESH,ALL
ALLSEL,ALL  
! Verical beams (2nd-3rd floor) 
LATT,1,1,1,,,,3
LSEL,S,LINE,,5,6,1,1
LMESH,ALL
ALLSEL,ALL  
! Verical beams (3rd-4th floor) 
LATT,1,1,1,,,,4
LSEL,S,LINE,,7,8,1,1
LMESH,ALL
ALLSEL,ALL  
! Verical beams (4th-5th floor) 
LATT,1,1,1,,,,5
LSEL,S,LINE,,9,10,1,1
LMESH,ALL
ALLSEL,ALL  

! Horizontal inner beams 
LATT,2,2,2,,,,6
LSEL,S,LINE,,11,15,1,1
LESIZE,ALL,,,1,1,1  					! Specifies the divisions (1) and spacing ratio on unmeshed lines (no division, no spacing ratio)
LMESH,ALL								
ALLSEL,ALL

! List nodes
NLIST,ALL,,,,NODE,X,Z
! List sections and their properties
SLIST,ALL,,,FULL
! List Elements and its properties
ELIST,ALL,,,0,1

! Read the Acceleration Input file                                                                                                                                      
! ================================
*DIM,EXCIT,ARRAY,NoS,1
*CREATE,ansuitmp                                                                                                                    
*VREAD,EXCIT,'testIN','txt',,IJK,NoS,1
(E16.7)                                                                                                                                                                
*END                                                                                                                                                                   
/INPUT,ansuitmp
FINISH

! =========================================================
! SOLUTION 
! =========================================================
/SOLU
SOLCONTROL,ON
! Rayleigh Damping
! ================
ALPHAD,4.00
BETAD,0.0002
	
! Transient Analysis Options
! ==========================
ANTYPE,TRANS			! Specifies the analysis type and restart status.
TRNOPT,FULL				! Specifies transient analysis options.
NLGEOM,ON				! Includes large-deflection effects in a static or full transient analysis.                                                                                                                                                           
NROPT,FULL				! Specifies the Newton-Raphson options in a static or full transient analysis.                                                                                                                                                                                                                                                      
EQSLV,SPARSE			! Specifies the type of equation solver.    
AUTOTS,ON				! Automatic time and/or load stepping

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
NSEL,ALL
D,ALL,UY				! Constraint all dofs for fixed supports
D,ALL,ROTX				! Constraint all dofs for fixed supports
D,ALL,ROTZ				! Constraint all dofs for fixed supports
ALLSEL,ALL

NSEL,S,LOC,Z,0
D,ALL,ALL				! Constraint all dofs for fixed supports
ALLSEL,ALL


STRd = 'SF2Dtest1_dspl'
STRv = 'SF2Dtest1_vel'
STRa = 'SF2Dtest1_acc'

TIMINT,OFF,STRUC 		! Static load step (dead weight)
KBC,1 
NSUBST,100
ACEL,0,0,9.81			! Acceleration of gravity 
TIME,1e-6 
SOLVE

TIMINT,ON				! Transient Analysis
AUTOTS,ON				! Automatic time and/or load stepping

*DO,I,1,NoS				! Time index DO LOOP
	! Time index
	Tindx = DT*(I)
	TIME,Tindx
	! X-axis input acceleration
	Gacc = EXCIT(I)
	! Application of the X-axis & gravity acceleration 
	ACEL,Gacc,0,9.81
	KBC,0					! Ramped load steps
	NSUBST,100				! Specifies the number of substeps to be taken this load step 
	SOLVE

		
	*DO,ElemI,1,1				! Time index DO LOOP
		STRx = STRCAT('SF2Dtest1_Sx_',CHRVAL(ElemI))
		STRy = STRCAT('SF2Dtest1_Sy_',CHRVAL(ElemI))
		STRz = STRCAT('SF2Dtest1_Sz_',CHRVAL(ElemI))
		STRxy = STRCAT('SF2Dtest1_Sxy_',CHRVAL(ElemI))
		STRyz = STRCAT('SF2Dtest1_Syz_',CHRVAL(ElemI))
		STRxz = STRCAT('SF2Dtest1_Sxz_',CHRVAL(ElemI))
		
		EPSx = STRCAT('SF2Dtest1_Ex_',CHRVAL(ElemI))
		EPSy = STRCAT('SF2Dtest1_Ey_',CHRVAL(ElemI))
		EPSz = STRCAT('SF2Dtest1_Ez_',CHRVAL(ElemI))
		EPSxy = STRCAT('SF2Dtest1_Exy_',CHRVAL(ElemI))
		EPSyz = STRCAT('SF2Dtest1_Eyz_',CHRVAL(ElemI))
		EPSxz = STRCAT('SF2Dtest1_Exz_',CHRVAL(ElemI))
			
		! X stresses at element 1, node I at 8 integration points  
		*GET,SxI5,SECR,ElemI,S,X,IVAL,5 	
		*GET,SxI6,SECR,ElemI,S,X,IVAL,6
		*GET,SxI32,SECR,ElemI,S,X,IVAL,32
		*GET,SxI31,SECR,ElemI,S,X,IVAL,31
		*GET,SxI13,SECR,ElemI,S,X,IVAL,13 	
		*GET,SxI16,SECR,ElemI,S,X,IVAL,16
		*GET,SxI22,SECR,ElemI,S,X,IVAL,22
		*GET,SxI23,SECR,ElemI,S,X,IVAL,23
		! X stresses at element 1, node J at 8 integration points  
		*GET,SxJ5,SECR,ElemI,S,X,JVAL,5 	
		*GET,SxJ6,SECR,ElemI,S,X,JVAL,6
		*GET,SxJ32,SECR,ElemI,S,X,JVAL,32
		*GET,SxJ31,SECR,ElemI,S,X,JVAL,31
		*GET,SxJ13,SECR,ElemI,S,X,JVAL,13 	
		*GET,SxJ16,SECR,ElemI,S,X,JVAL,16
		*GET,SxJ22,SECR,ElemI,S,X,JVAL,22
		*GET,SxJ23,SECR,ElemI,S,X,JVAL,23
		
	
		! Y stresses at element 1, node I at 8 integration points  
		*GET,SyI5,SECR,ElemI,S,Y,IVAL,5 	
		*GET,SyI6,SECR,ElemI,S,Y,IVAL,6
		*GET,SyI32,SECR,ElemI,S,Y,IVAL,32
		*GET,SyI31,SECR,ElemI,S,Y,IVAL,31
		*GET,SyI13,SECR,ElemI,S,Y,IVAL,13
		*GET,SyI16,SECR,ElemI,S,Y,IVAL,16
		*GET,SyI22,SECR,ElemI,S,Y,IVAL,22
		*GET,SyI23,SECR,ElemI,S,Y,IVAL,23
		! Y stresses at element 1, node J at 8 integration points  
		*GET,SyJ5,SECR,ElemI,S,Y,JVAL,5 	
		*GET,SyJ6,SECR,ElemI,S,Y,JVAL,6
		*GET,SyJ32,SECR,ElemI,S,Y,JVAL,32
		*GET,SyJ31,SECR,ElemI,S,Y,JVAL,31
		*GET,SyJ13,SECR,ElemI,S,Y,JVAL,13 	
		*GET,SyJ16,SECR,ElemI,S,Y,JVAL,16
		*GET,SyJ22,SECR,ElemI,S,Y,JVAL,22
		*GET,SyJ23,SECR,ElemI,S,Y,JVAL,23
	
		
		! Z stresses at element 1, node I at 8 integration points  
		*GET,SzI5,SECR,ElemI,S,Z,IVAL,5 	
		*GET,SzI6,SECR,ElemI,S,Z,IVAL,6
		*GET,SzI32,SECR,ElemI,S,Z,IVAL,32
		*GET,SzI31,SECR,ElemI,S,Z,IVAL,31
		*GET,SzI13,SECR,ElemI,S,Z,IVAL,13
		*GET,SzI16,SECR,ElemI,S,Z,IVAL,16
		*GET,SzI22,SECR,ElemI,S,Z,IVAL,22
		*GET,SzI23,SECR,ElemI,S,Z,IVAL,23
		! Z stresses at element 1, node J at 8 integration points  
		*GET,SzJ5,SECR,ElemI,S,Z,JVAL,5 	
		*GET,SzJ6,SECR,ElemI,S,Z,JVAL,6
		*GET,SzJ32,SECR,ElemI,S,Z,JVAL,32
		*GET,SzJ31,SECR,ElemI,S,Z,JVAL,31
		*GET,SzJ13,SECR,ElemI,S,Z,JVAL,13 	
		*GET,SzJ16,SECR,ElemI,S,Z,JVAL,16
		*GET,SzJ22,SECR,ElemI,S,Z,JVAL,22
		*GET,SzJ23,SECR,ElemI,S,Z,JVAL,23
		
		! XY shear-stresses at element 1, node I at 8 integration points  
		*GET,SxyI5,SECR,ElemI,S,XY,IVAL,5 	
		*GET,SxyI6,SECR,ElemI,S,XY,IVAL,6
		*GET,SxyI32,SECR,ElemI,S,XY,IVAL,32
		*GET,SxyI31,SECR,ElemI,S,XY,IVAL,31
		*GET,SxyI13,SECR,ElemI,S,XY,IVAL,13
		*GET,SxyI16,SECR,ElemI,S,XY,IVAL,16
		*GET,SxyI22,SECR,ElemI,S,XY,IVAL,22
		*GET,SxyI23,SECR,ElemI,S,XY,IVAL,23
		! XY shear-stresses at element 1, node J at 8 integration points  
		*GET,SxyJ5,SECR,ElemI,S,XY,JVAL,5 		
		*GET,SxyJ6,SECR,ElemI,S,XY,JVAL,6
		*GET,SxyJ32,SECR,ElemI,S,XY,JVAL,32
		*GET,SxyJ31,SECR,ElemI,S,XY,JVAL,31
		*GET,SxyJ13,SECR,ElemI,S,XY,JVAL,13
		*GET,SxyJ16,SECR,ElemI,S,XY,JVAL,16
		*GET,SxyJ22,SECR,ElemI,S,XY,JVAL,22
		*GET,SxyJ23,SECR,ElemI,S,XY,JVAL,23
					
		! YZ shear-stresses at element 1, node I at 8 integration points  		
		*GET,SyzI5,SECR,ElemI,S,YZ,IVAL,5		
		*GET,SyzI6,SECR,ElemI,S,YZ,IVAL,6
		*GET,SyzI32,SECR,ElemI,S,YZ,IVAL,32
		*GET,SyzI31,SECR,ElemI,S,YZ,IVAL,31
		*GET,SyzI13,SECR,ElemI,S,YZ,IVAL,13		
		*GET,SyzI16,SECR,ElemI,S,YZ,IVAL,16
		*GET,SyzI22,SECR,ElemI,S,YZ,IVAL,22
		*GET,SyzI23,SECR,ElemI,S,YZ,IVAL,23
		! YZ shear-stresses at element 1, node J at 8 integration points  
		*GET,SyzJ5,SECR,ElemI,S,YZ,JVAL,5 
		*GET,SyzJ6,SECR,ElemI,S,YZ,JVAL,6
		*GET,SyzJ32,SECR,ElemI,S,YZ,JVAL,32
		*GET,SyzJ31,SECR,ElemI,S,YZ,JVAL,31
		*GET,SyzJ13,SECR,ElemI,S,YZ,JVAL,13
		*GET,SyzJ16,SECR,ElemI,S,YZ,JVAL,16
		*GET,SyzJ22,SECR,ElemI,S,YZ,JVAL,22
		*GET,SyzJ23,SECR,ElemI,S,YZ,JVAL,23
		
		! XZ shear-stresses at element 1, node I at 8 integration points  
		*GET,SxzI5,SECR,ElemI,S,XZ,IVAL,5
		*GET,SxzI6,SECR,ElemI,S,XZ,IVAL,6
		*GET,SxzI32,SECR,ElemI,S,XZ,IVAL,32
		*GET,SxzI31,SECR,ElemI,S,XZ,IVAL,31
		*GET,SxzI13,SECR,ElemI,S,XZ,IVAL,13
		*GET,SxzI16,SECR,ElemI,S,XZ,IVAL,16
		*GET,SxzI22,SECR,ElemI,S,XZ,IVAL,22
		*GET,SxzI23,SECR,ElemI,S,XZ,IVAL,23
		! XZ shear-stresses at element 1, node J at 8 integration points  
		*GET,SxzJ5,SECR,ElemI,S,XZ,JVAL,5
		*GET,SxzJ6,SECR,ElemI,S,XZ,JVAL,6
		*GET,SxzJ32,SECR,ElemI,S,XZ,JVAL,32
		*GET,SxzJ31,SECR,ElemI,S,XZ,JVAL,31
		*GET,SxzJ13,SECR,ElemI,S,XZ,JVAL,13 	! X stress at element 1		
		*GET,SxzJ16,SECR,ElemI,S,XZ,JVAL,16
		*GET,SxzJ22,SECR,ElemI,S,XZ,JVAL,22
		*GET,SxzJ23,SECR,ElemI,S,XZ,JVAL,23
		
		
		
		! X strains at element 1, node I at 8 integration points  
		*GET,ExI5,SECR,ElemI,EPTO,X,IVAL,5 	
		*GET,ExI6,SECR,ElemI,EPTO,X,IVAL,6
		*GET,ExI32,SECR,ElemI,EPTO,X,IVAL,32
		*GET,ExI31,SECR,ElemI,EPTO,X,IVAL,31
		*GET,ExI13,SECR,ElemI,EPTO,X,IVAL,13 	
		*GET,ExI16,SECR,ElemI,EPTO,X,IVAL,16
		*GET,ExI22,SECR,ElemI,EPTO,X,IVAL,22
		*GET,ExI23,SECR,ElemI,EPTO,X,IVAL,23
		! X stresses at element 1, node J at 8 integration points  
		*GET,ExJ5,SECR,ElemI,EPTO,X,JVAL,5 	
		*GET,ExJ6,SECR,ElemI,EPTO,X,JVAL,6
		*GET,ExJ32,SECR,ElemI,EPTO,X,JVAL,32
		*GET,ExJ31,SECR,ElemI,EPTO,X,JVAL,31
		*GET,ExJ13,SECR,ElemI,EPTO,X,JVAL,13 	
		*GET,ExJ16,SECR,ElemI,EPTO,X,JVAL,16
		*GET,ExJ22,SECR,ElemI,EPTO,X,JVAL,22
		*GET,ExJ23,SECR,ElemI,EPTO,X,JVAL,23
		
		! Y strains at element 1, node I at 8 integration points  
		*GET,EyI5,SECR,ElemI,EPTO,Y,IVAL,5 	
		*GET,EyI6,SECR,ElemI,EPTO,Y,IVAL,6
		*GET,EyI32,SECR,ElemI,EPTO,Y,IVAL,32
		*GET,EyI31,SECR,ElemI,EPTO,Y,IVAL,31
		*GET,EyI13,SECR,ElemI,EPTO,Y,IVAL,13
		*GET,EyI16,SECR,ElemI,EPTO,Y,IVAL,16
		*GET,EyI22,SECR,ElemI,EPTO,Y,IVAL,22
		*GET,EyI23,SECR,ElemI,EPTO,Y,IVAL,23
		! Y strains at element 1, node J at 8 integration points  
		*GET,EyJ5,SECR,ElemI,EPTO,Y,JVAL,5 	
		*GET,EyJ6,SECR,ElemI,EPTO,Y,JVAL,6
		*GET,EyJ32,SECR,ElemI,EPTO,Y,JVAL,32
		*GET,EyJ31,SECR,ElemI,EPTO,Y,JVAL,31
		*GET,EyJ13,SECR,ElemI,EPTO,Y,JVAL,13 	
		*GET,EyJ16,SECR,ElemI,EPTO,Y,JVAL,16
		*GET,EyJ22,SECR,ElemI,EPTO,Y,JVAL,22
		*GET,EyJ23,SECR,ElemI,EPTO,Y,JVAL,23
		
		! Z strains at element 1, node I at 8 integration points  
		*GET,EzI5,SECR,ElemI,EPTO,Z,IVAL,5 	
		*GET,EzI6,SECR,ElemI,EPTO,Z,IVAL,6
		*GET,EzI32,SECR,ElemI,EPTO,Z,IVAL,32
		*GET,EzI31,SECR,ElemI,EPTO,Z,IVAL,31
		*GET,EzI13,SECR,ElemI,EPTO,Z,IVAL,13
		*GET,EzI16,SECR,ElemI,EPTO,Z,IVAL,16
		*GET,EzI22,SECR,ElemI,EPTO,Z,IVAL,22
		*GET,EzI23,SECR,ElemI,EPTO,Z,IVAL,23
		! Z strains at element 1, node J at 8 integration points  
		*GET,EzJ5,SECR,ElemI,EPTO,Z,JVAL,5 	
		*GET,EzJ6,SECR,ElemI,EPTO,Z,JVAL,6
		*GET,EzJ32,SECR,ElemI,EPTO,Z,JVAL,32
		*GET,EzJ31,SECR,ElemI,EPTO,Z,JVAL,31
		*GET,EzJ13,SECR,ElemI,EPTO,Z,JVAL,13 	
		*GET,EzJ16,SECR,ElemI,EPTO,Z,JVAL,16
		*GET,EzJ22,SECR,ElemI,EPTO,Z,JVAL,22
		*GET,EzJ23,SECR,ElemI,EPTO,Z,JVAL,23
		
		! XY shear-strains at element 1, node I at 8 integration points  
		*GET,ExyI5,SECR,ElemI,EPTO,XY,IVAL,5 	
		*GET,ExyI6,SECR,ElemI,EPTO,XY,IVAL,6
		*GET,ExyI32,SECR,ElemI,EPTO,XY,IVAL,32
		*GET,ExyI31,SECR,ElemI,EPTO,XY,IVAL,31
		*GET,ExyI13,SECR,ElemI,EPTO,XY,IVAL,13
		*GET,ExyI16,SECR,ElemI,EPTO,XY,IVAL,16
		*GET,ExyI22,SECR,ElemI,EPTO,XY,IVAL,22
		*GET,ExyI23,SECR,ElemI,EPTO,XY,IVAL,23
		! XY shear-strains at element 1, node J at 8 integration points  
		*GET,ExyJ5,SECR,ElemI,EPTO,XY,JVAL,5 		
		*GET,ExyJ6,SECR,ElemI,EPTO,XY,JVAL,6
		*GET,ExyJ32,SECR,ElemI,EPTO,XY,JVAL,32
		*GET,ExyJ31,SECR,ElemI,EPTO,XY,JVAL,31
		*GET,ExyJ13,SECR,ElemI,EPTO,XY,JVAL,13
		*GET,ExyJ16,SECR,ElemI,EPTO,XY,JVAL,16
		*GET,ExyJ22,SECR,ElemI,EPTO,XY,JVAL,22
		*GET,ExyJ23,SECR,ElemI,EPTO,XY,JVAL,23
				
		! YZ shear-strains at element 1, node I at 8 integration points  		
		*GET,EyzI5,SECR,ElemI,EPTO,YZ,IVAL,5		
		*GET,EyzI6,SECR,ElemI,EPTO,YZ,IVAL,6
		*GET,EyzI32,SECR,ElemI,EPTO,YZ,IVAL,32
		*GET,EyzI31,SECR,ElemI,EPTO,YZ,IVAL,31
		*GET,EyzI13,SECR,ElemI,EPTO,YZ,IVAL,13		
		*GET,EyzI16,SECR,ElemI,EPTO,YZ,IVAL,16
		*GET,EyzI22,SECR,ElemI,EPTO,YZ,IVAL,22
		*GET,EyzI23,SECR,ElemI,EPTO,YZ,IVAL,23
		! YZ shear-strains at element 1, node J at 8 integration points  
		*GET,EyzJ5,SECR,ElemI,EPTO,YZ,JVAL,5 
		*GET,EyzJ6,SECR,ElemI,EPTO,YZ,JVAL,6
		*GET,EyzJ32,SECR,ElemI,EPTO,YZ,JVAL,32
		*GET,EyzJ31,SECR,ElemI,EPTO,YZ,JVAL,31
		*GET,EyzJ13,SECR,ElemI,EPTO,YZ,JVAL,13
		*GET,EyzJ16,SECR,ElemI,EPTO,YZ,JVAL,16
		*GET,EyzJ22,SECR,ElemI,EPTO,YZ,JVAL,22
		*GET,EyzJ23,SECR,ElemI,EPTO,YZ,JVAL,23
		
		! XZ shear-strains at element 1, node I at 8 integration points  
		*GET,ExzI5,SECR,ElemI,EPTO,XZ,IVAL,5
		*GET,ExzI6,SECR,ElemI,EPTO,XZ,IVAL,6
		*GET,ExzI32,SECR,ElemI,EPTO,XZ,IVAL,32
		*GET,ExzI31,SECR,ElemI,EPTO,XZ,IVAL,31
		*GET,ExzI13,SECR,ElemI,EPTO,XZ,IVAL,13
		*GET,ExzI16,SECR,ElemI,EPTO,XZ,IVAL,16
		*GET,ExzI22,SECR,ElemI,EPTO,XZ,IVAL,22
		*GET,ExzI23,SECR,ElemI,EPTO,XZ,IVAL,23
		! XZ shear-strains at element 1, node J at 8 integration points  
		*GET,ExzJ5,SECR,ElemI,EPTO,XZ,JVAL,5
		*GET,ExzJ6,SECR,ElemI,EPTO,XZ,JVAL,6
		*GET,ExzJ32,SECR,ElemI,EPTO,XZ,JVAL,32
		*GET,ExzJ31,SECR,ElemI,EPTO,XZ,JVAL,31
		*GET,ExzJ13,SECR,ElemI,EPTO,XZ,JVAL,13 	! X stress at element 1		
		*GET,ExzJ16,SECR,ElemI,EPTO,XZ,JVAL,16
		*GET,ExzJ22,SECR,ElemI,EPTO,XZ,JVAL,22
		*GET,ExzJ23,SECR,ElemI,EPTO,XZ,JVAL,23
		

		*CFOPEN,STRx,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SxI5,SxI6,SxI32,SxI31,SxI13,SxI16,SxI22,SxI23,SxJ5,SxJ6,SxJ32,SxJ31,SxJ13,SxJ16,SxJ22,SxJ23
(18(E16.7,2X))
		*CFOPEN,STRy,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SyI5,SyI6,SyI32,SyI31,SyI13,SyI16,SyI22,SyI23,SyJ5,SyJ6,SyJ32,SyJ31,SyJ13,SyJ16,SyJ22,SyJ23
(18(E16.7,2X))
		*CFOPEN,STRz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SzI5,SzI6,SzI32,SzI31,SzI13,SzI16,SzI22,SzI23,SzJ5,SzJ6,SzJ32,SzJ31,SzJ13,SzJ16,SzJ22,SzJ23
(18(E16.7,2X))
		*CFOPEN,STRxy,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SxyI5,SxyI6,SxyI32,SxyI31,SxyI13,SxyI16,SxyI22,SxyI23,SxyJ5,SxyJ6,SxyJ32,SxyJ31,SxyJ13,SxyJ16,SxyJ22,SxyJ23
(18(E16.7,2X))
		*CFOPEN,STRyz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SyzI5,SyzI6,SyzI32,SyzI31,SyzI13,SyzI16,SyzI22,SyzI23,SyzJ5,SyzJ6,SyzJ32,SyzJ31,SyzJ13,SyzJ16,SyzJ22,SyzJ23
(18(E16.7,2X))
		*CFOPEN,STRxz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,SxzI5,SxzI6,SxzI32,SxzI31,SxzI13,SxzI16,SxzI22,SxzI23,SxzJ5,SxzJ6,SxzJ32,SxzJ31,SxzJ13,SxzJ16,SxzJ22,SxzJ23
(18(E16.7,2X))				! Close 'SFresI' txt file


		*CFOPEN,EPSx,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,ExI5,ExI6,ExI32,ExI31,ExI13,ExI16,ExI22,ExI23,ExJ5,ExJ6,ExJ32,ExJ31,ExJ13,ExJ16,ExJ22,ExJ23
(18(E16.7,2X))
		*CFOPEN,EPSy,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,EyI5,EyI6,EyI32,EyI31,EyI13,EyI16,EyI22,EyI23,EyJ5,EyJ6,EyJ32,EyJ31,EyJ13,EyJ16,EyJ22,EyJ23
(18(E16.7,2X))
		*CFOPEN,EPSz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,EzI5,EzI6,EzI32,EzI31,EzI13,EzI16,EzI22,EzI23,EzJ5,EzJ6,EzJ32,EzJ31,EzJ13,EzJ16,EzJ22,EzJ23
(18(E16.7,2X))
		*CFOPEN,EPSxy,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,ExyI5,ExyI6,ExyI32,ExyI31,ExyI13,ExyI16,ExyI22,ExyI23,ExyJ5,ExyJ6,ExyJ32,ExyJ31,ExyJ13,ExyJ16,ExyJ22,ExyJ23
(18(E16.7,2X))
		*CFOPEN,EPSyz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,EyzI5,EyzI6,EyzI32,EyzI31,EyzI13,EyzI16,EyzI22,EyzI23,EyzJ5,EyzJ6,EyzJ32,EyzJ31,EyzJ13,EyzJ16,EyzJ22,EyzJ23
(18(E16.7,2X))
		*CFOPEN,EPSxz,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,ExzI5,ExzI6,ExzI32,ExzI31,ExzI13,ExzI16,ExzI22,ExzI23,ExzJ5,ExzJ6,ExzJ32,ExzJ31,ExzJ13,ExzJ16,ExzJ22,ExzJ23
(18(E16.7,2X))

	
	*ENDDO
	
	*GET,Ux1r,NODE,2,U,X 		! X displacement at node 2
	*GET,Ux2r,NODE,5,U,X 		! X displacement at node 5		
	*GET,Ux3r,NODE,7,U,X 		! X displacement at node 7		
	*GET,Ux4r,NODE,9,U,X 		! X displacement at node 9		
	*GET,Ux5r,NODE,11,U,X 		! X displacement at node 11		
	*GET,Ux1l,NODE,4,U,X 		! X displacement at node 4
	*GET,Ux2l,NODE,6,U,X 		! X displacement at node 6		
	*GET,Ux3l,NODE,8,U,X 		! X displacement at node 8		
	*GET,Ux4l,NODE,10,U,X 		! X displacement at node 10		
	*GET,Ux5l,NODE,12,U,X 		! X displacement at node 12		
	
	*GET,Vx1r,NODE,2,V,X 		! X velocity at node 2
	*GET,Vx2r,NODE,5,V,X 		! X velocity at node 5		
	*GET,Vx3r,NODE,7,V,X 		! X velocity at node 7		
	*GET,Vx4r,NODE,9,V,X 		! X velocity at node 9		
	*GET,Vx5r,NODE,11,V,X 		! X velocity at node 11		
	*GET,Vx1l,NODE,4,V,X 		! X velocity at node 4
	*GET,Vx2l,NODE,6,V,X 		! X velocity at node 6		
	*GET,Vx3l,NODE,8,V,X 		! X velocity at node 8		
	*GET,Vx4l,NODE,10,V,X 		! X velocity at node 10		
	*GET,Vx5l,NODE,12,V,X 		! X velocity at node 12		
		
	*GET,Ax1r,NODE,2,A,X		! X acceleration at node 2
	*GET,Ax2r,NODE,5,A,X		! X acceleration at node 5
	*GET,Ax3r,NODE,7,A,X		! X acceleration at node 7
	*GET,Ax4r,NODE,9,A,X		! X acceleration at node 9
	*GET,Ax5r,NODE,11,A,X		! X acceleration at node 11
	*GET,Ax1l,NODE,4,A,X		! X acceleration at node 4
	*GET,Ax2l,NODE,6,A,X		! X acceleration at node 6
	*GET,Ax3l,NODE,8,A,X		! X acceleration at node 8
	*GET,Ax4l,NODE,10,A,X		! X acceleration at node 10
	*GET,Ax5l,NODE,12,A,X		! X acceleration at node 12
		
	*GET,ReacFx1,NODE,1,RF,FX 	! X Reaction Force at node 1
	*GET,ReacFx2,NODE,3,RF,FX 	! X Reaction Force at node 2
		
	*CFOPEN,STRd,txt,,APPEND		! Open txt file
	*VWRITE,Tindx,Gacc,ReacFx1+ReacFx2,Ux1r,Ux2r,Ux3r,Ux4r,Ux5r,Ux1l,Ux2l,Ux3l,Ux4l,Ux5l
(13(E16.7,2X))
	*CFCLOS					! Close 'SFresI' txt file	
	*CFOPEN,STRv,txt,,APPEND		! Open txt file
	*VWRITE,Tindx,Gacc,ReacFx1+ReacFx2,Vx1r,Vx2r,Vx3r,Vx4r,Vx5r,Vx1l,Vx2l,Vx3l,Vx4l,Vx5l
(13(E16.7,2X))
	*CFCLOS					! Close 'SFresI' txt file	
	*CFOPEN,STRa,txt,,APPEND		! Open txt file
	*VWRITE,Tindx,Gacc,ReacFx1+ReacFx2,Ax1r,Ax2r,Ax3r,Ax4r,Ax5r,Ax1l,Ax2l,Ax3l,Ax4l,Ax5l
(13(E16.7,2X))
	*CFCLOS					! Close 'SFresI' txt file
*ENDDO
FINISH