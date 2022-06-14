/CLEAR
/CWD 'D:\SpaceFrame\'
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
*SET,NoS,1000
*SET,NoExp,50
*SET,DT,0.025
*SET,UNCERT,1				! Material uncertainty flag 

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
LESIZE,ALL,,,3,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
! Verical beams (ground-1st floor) 
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

! If material uncertainty flag == 1
! =======================================                                                                                                                                      
*IF,UNCERT,EQ,1,THEN
	! Read the Material Properties Input file
	! =======================================                                                                                                                                      
	*DIM,BeamProps,ARRAY,NoExp,1 
	*VREAD,BeamProps(1),'BeamProps','txt',,IJK,NoExp,1
(E16.7)                                                                                                                                                                
*ENDIF

! Read the Acceleration Input file                                                                                                                                      
! ================================
*DIM,EXCIT,ARRAY,NoExp*NoS,1
*CREATE,ansuitmp                                                                                                                    
*VREAD,EXCIT,'AccEQ','txt',,IJK,NoExp*NoS,1
(E16.7)                                                                                                                                                                
*END                                                                                                                                                                   
/INPUT,ansuitmp

FINISH

! =========================================================
! SOLUTION 
! =========================================================
*DO,J,1,NoExp			! Experiments DO LOOP
	
	! PREP7 
	! =====================================================
	*IF,UNCERT,EQ,1,THEN
	
	/DELETE,'MProps',txt,
	/PREP7
	MPDELE,ALL,1						! Delete material properties
	TBDELE,ALL,1						! Delete nonlinear material properties
	
	! Material 1: Vertical beams 
	MP,EX,1,200E9*BeamProps(J,1)        ! N/m^2	
	MP,PRXY,1,0.29
	MP,DENS,1,7850
	TB,BISO,1	 						! Redefine bilinear isotropic hardening  
	TBDATA,1,200e6,1E10
	KEYOPT,1,3,3
	
	FINISH

	! Write material properties for the current experiment
	! =====================================================
	*CFOPEN,'MProps',txt,,APPEND			! Open txt file
	*GET,MatEx,EX,1							! Material 1 Young moduli
	*GET,YieldStress,BISO,1,TEMP,0,CONST,1	! Material 1 yiels stress
	*GET,BEAMheight1,SECP,1,DATA,1			! Beam 1 cross section area
	*GET,BEAMheight2,SECP,2,DATA,1			! Beam 2 cross section area
	*GET,BEAMheight3,SECP,3,DATA,1			! Beam 1 cross section area
	*GET,BEAMheight4,SECP,4,DATA,1			! Beam 2 cross section area
	*GET,BEAMheight5,SECP,5,DATA,1			! Beam 1 cross section area
	*GET,BEAMheight6,SECP,6,DATA,1			! Beam 2 cross section area
	*VWRITE,J,MatEx,YieldStress,BEAMheight1,BEAMheight2,BEAMheight3,BEAMheight4,BEAMheight5,BEAMheight6
(9(E16.7,2X))
	*CFCLOS								! Close 'MaterialProps' txt file

	*ENDIF
	
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
	NSEL,S,LOC,Z,0
	D,ALL,ALL				! Constraint all dofs for fixed supports
	ALLSEL,ALL
	
	STRd = STRCAT('SF2D_EQ_dspl',CHRVAL(J))
	STRv = STRCAT('SF2D_EQ_vel',CHRVAL(J))
	STRa = STRCAT('SF2D_EQ_acc',CHRVAL(J))
	STRs = STRCAT('SF2D_EQ_stress',CHRVAL(J))
	
	
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
		Gacc = EXCIT((J-1)*NoS+I)
		! Application of the X-axis & gravity acceleration 
		ACEL,Gacc,0,9.81

		KBC,0					! Ramped load steps
		NSUBST,100				! Specifies the number of substeps to be taken this load step 
		SOLVE
	
		
		*GET,Ux1r,NODE,2,U,X 		! X displacement at node 2
		*GET,Ux2r,NODE,9,U,X 		! X displacement at node 5		
		*GET,Ux3r,NODE,15,U,X 		! X displacement at node 7		
		*GET,Ux4r,NODE,21,U,X 		! X displacement at node 9		
		*GET,Ux5r,NODE,27,U,X 		! X displacement at node 11		
		*GET,Ux1l,NODE,6,U,X 		! X displacement at node 4
		*GET,Ux2l,NODE,12,U,X 		! X displacement at node 6		
		*GET,Ux3l,NODE,18,U,X 		! X displacement at node 8		
		*GET,Ux4l,NODE,24,U,X 		! X displacement at node 10		
		*GET,Ux5l,NODE,30,U,X 		! X displacement at node 12		
		
		*GET,Vx1r,NODE,2,V,X 		! X velocity at node 2
		*GET,Vx2r,NODE,9,V,X 		! X velocity at node 5		
		*GET,Vx3r,NODE,15,V,X 		! X velocity at node 7		
		*GET,Vx4r,NODE,21,V,X 		! X velocity at node 9		
		*GET,Vx5r,NODE,27,V,X 		! X velocity at node 11		
		*GET,Vx1l,NODE,6,V,X 		! X velocity at node 4
		*GET,Vx2l,NODE,12,V,X 		! X velocity at node 6		
		*GET,Vx3l,NODE,18,V,X 		! X velocity at node 8		
		*GET,Vx4l,NODE,24,V,X 		! X velocity at node 10		
		*GET,Vx5l,NODE,30,V,X 		! X velocity at node 12		
		
		*GET,Ax1r,NODE,2,A,X		! X acceleration at node 2
		*GET,Ax2r,NODE,9,A,X		! X acceleration at node 5
		*GET,Ax3r,NODE,15,A,X		! X acceleration at node 7
		*GET,Ax4r,NODE,21,A,X		! X acceleration at node 9
		*GET,Ax5r,NODE,27,A,X		! X acceleration at node 11
		*GET,Ax1l,NODE,6,A,X		! X acceleration at node 4
		*GET,Ax2l,NODE,12,A,X		! X acceleration at node 6
		*GET,Ax3l,NODE,18,A,X		! X acceleration at node 8
		*GET,Ax4l,NODE,24,A,X		! X acceleration at node 10
		*GET,Ax5l,NODE,30,A,X		! X acceleration at node 12
		
		*GET,Sx1,SECR,1,S,EQV,MAX 	! X stress at element 1		
		*GET,Sx2,SECR,4,S,EQV,MAX 	! X stress at element 2
		*GET,Sx3,SECR,7,S,EQV,MAX 	! X stress at element 3
		*GET,Sx4,SECR,10,S,EQV,MAX 	! X stress at element 4
		*GET,Sx5,SECR,13,S,EQV,MAX 	! X stress at element 5		
		*GET,Sx6,SECR,16,S,EQV,MAX 	! X stress at element 6		
		*GET,Sx7,SECR,19,S,EQV,MAX 	! X stress at element 7
		*GET,Sx8,SECR,22,S,EQV,MAX 	! X stress at element 8
		*GET,Sx9,SECR,25,S,EQV,MAX 	! X stress at element 9
		*GET,Sx10,SECR,28,S,EQV,MAX ! X stress at element 10		
		
		*GET,ReacFx1,NODE,1,RF,FX 	! X Reaction Force at node 1
		*GET,ReacFx2,NODE,5,RF,FX 	! X Reaction Force at node 2
		
		*CFOPEN,STRs,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,ReacFx1+ReacFx2,Sx1,Sx2,Sx3,Sx4,Sx5,Sx6,Sx7,Sx8,Sx9,Sx10
(13(E16.7,2X))
		*CFCLOS					! Close 'SFresI' txt file
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
*ENDDO