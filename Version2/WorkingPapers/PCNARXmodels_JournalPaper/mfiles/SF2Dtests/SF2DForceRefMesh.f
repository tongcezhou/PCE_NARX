/CLEAR
/CWD 'D:\SpaceFrame\'
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
*SET,NoS,2500
*SET,DT,0.1

! Material 1: Vertical beams 
ET,1,BEAM188
MP,EX,1,200E9
MP,PRXY,1,0.29
MP,DENS,1,7850
KEYOPT,1,3,3
! Nonlinear material properties
TB,BISO,1								! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,200E6,10E9

! Cross-section 1 (vertical elements)
SECTYPE,1,BEAM,RECT
SECDATA,0.20,0.20,3,3
! Cross-section 2 (vertical elements)
SECTYPE,2,BEAM,RECT
SECDATA,0.175,0.175,3,3
! Cross-section 3 (vertical elements)
SECTYPE,3,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 4 (vertical elements)
SECTYPE,4,BEAM,RECT
SECDATA,0.125,0.125,3,3
! Cross-section 5 (vertical elements)
SECTYPE,5,BEAM,RECT
SECDATA,0.10,0.10,3,3
! Cross-section 6 (horizontal outer elements)
SECTYPE,6,BEAM,RECT
SECDATA,0.15,0.15,3,3
SECCONTROLS,,,,6524

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
! Verical beams (ground-1st floor) 
LESIZE,ALL,,,3,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
LATT,1,1,1,,,,1
LSEL,S,LINE,,1,2,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (1st-2nd floor) 
LATT,1,1,1,,,,2
LSEL,S,LINE,,3,4,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (2nd-3rd floor) 
LATT,1,1,1,,,,3
LSEL,S,LINE,,5,6,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (3rd-4th floor) 
LATT,1,1,1,,,,4
LSEL,S,LINE,,7,8,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (4th-5th floor) 
LATT,1,1,1,,,,5
LSEL,S,LINE,,9,10,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Horizontal inner beams 
LATT,1,1,1,,,,6
LSEL,S,LINE,,11,15,1,1
LESIZE,ALL,,,1,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL

! Read the Acceleration Input file                                                                                                                                      
! ================================
*DIM,EXCIT,ARRAY,NoS,1
*CREATE,ansuitmp                                                                                                                    
*VREAD,EXCIT,'InputForce','txt',,IJK,NoS,1
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
ALPHAD,6.1
BETAD,0.00028

! Transient Analysis Options
! ==========================
ANTYPE,TRANS			! Specifies the analysis type and restart status.
TRNOPT,FULL				! Specifies transient analysis options.
EQSLV,SPARSE         	! Specifies the type of equation solver
NLGEOM,ON				! Includes large-deflection effects in a static or full transient analysis.                                                                                                                                                           
NROPT,FULL				! Specifies the Newton-Raphson options in a static or full transient analysis.                                                                                                                                                                                                                                                      
AUTOTS,ON				! Automatic time and/or load stepping

! Controls the solution data written to the database.                                                                                                                  
OUTRES,ALL,NONE 
OUTRES,ESOL,LAST
OUTRES,NSOL,LAST
OUTRES,NLOAD,LAST                                                                                                                                                       
OUTRES,RSOL,LAST      
OUTRES,A,LAST 

! DOF constraints
! ===============
NSEL,S,LOC,Z,0
D,ALL,ALL				! Constraint all dofs for fixed supports
ALLSEL,ALL

*CFOPEN,SF_ForceResponse,txt,		! Open txt file

TIMINT,OFF,STRUC 		! Static load step (dead weight)
TIME,1e-6 

! Input Force
!Finput = EXCIT(1)
!F,2,FX,Finput,,4,2		! 1st floor
!F,5,FX,Finput*2,,6,1	! 2nd floor
!F,7,FX,Finput*3,,8,1	! 3rd floor
!F,9,FX,Finput*4,,10,1	! 4th floor
!F,11,FX,Finput*5,,12,1	! 5th floor
!ALLSEL,ALL

ACEL,0,0,9.81
NSUBST,5
KBC,1
SOLVE

TIMINT,ON				! Transient Analysis
AUTOTS,ON				! Automatic time and/or load stepping

*DO,I,1,NoS				! Time index DO LOOP
		
	! DELTIM,0.001			! Specifies the time step sizes
	! Time index
	Tindx = DT*(I)
	TIME,Tindx
	
	! Input Force
	Finput = EXCIT(I)
	F,2,FX,Finput,,8,6		! 1st floor
	F,13,FX,Finput*2,,18,5	! 2nd floor
	F,23,FX,Finput*3,,28,5	! 3rd floor
	F,33,FX,Finput*4,,38,5	! 4th floor
	F,43,FX,Finput*5,,48,5	! 5th floor
	NSEL,ALL
	
	ACEL,0,0,9.81
	ALLSEL
	KBC,0					! Ramped load steps
	! DELTIME,DT,DT/10,DT	
	NSUBST,100				! Specifies the number of substeps to be taken this load step 
	SOLVE

	*GET,Sx1,SECR,1,S,X,MAX 	! X displacement at element 1		
	*GET,Sx2,SECR,7,S,X,MAX 	! X displacement at element 3
	*GET,Sx3,SECR,13,S,X,MAX 	! X displacement at element 5
	*GET,Sx4,SECR,19,S,X,MAX 	! X displacement at element 7
	*GET,Sx5,SECR,25,S,X,MAX 	! X displacement at element 9

	*GET,Sl1,SECR,1,S,EQV,MAX 	! 1st principal stress at element 1		
	*GET,Sl2,SECR,7,S,EQV,MAX 	! 1st principal stress at element 2		
	*GET,Sl3,SECR,13,S,EQV,MAX 	! 1st principal stress at element 3		
	*GET,Sl4,SECR,19,S,EQV,MAX 	! 1st principal stress at element 4		
	*GET,Sl5,SECR,25,S,EQV,MAX 	! 1st principal stress at element 5		
	
	*GET,Ux1,NODE,2,U,X 		! X displacement at node 4
	*GET,Ux2,NODE,9,U,X 		! X displacement at node 6		
	*GET,Ux3,NODE,15,U,X 		! X displacement at node 8		
	*GET,Ux4,NODE,21,U,X 		! X displacement at node 10		
	*GET,Ux5,NODE,27,U,X 		! X displacement at node 12		
		
	*GET,ReacFx1,NODE,1,RF,FX 	! X Reaction Force at node 1
	*GET,ReacFx2,NODE,5,RF,FX 	! X Reaction Force at node 2
		
	*VWRITE,Tindx,Finput,ReacFx1+ReacFx2,Sl1,Sl2,Sl3,Sl4,Sl5,Sx1,Sx2,Sx3,Sx4,Sx5,Ux1,Ux2,Ux3,Ux4,Ux5
(18(E16.7,2X))

*ENDDO
*CFCLOS					! Close 'SpaceFrameForceResponse' txt file

FINISH