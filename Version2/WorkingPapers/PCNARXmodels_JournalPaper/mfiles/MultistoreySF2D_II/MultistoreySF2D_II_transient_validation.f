/CLEAR
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
*SET,NoS,1000
*SET,DT,0.02
*SET,NoExp,49
*SET,UNCERT,1	

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
TBDATA,1,200E6,20E9
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


! List nodes
NLIST,ALL,,,,NODE
! List Elements and its properties
ELIST,ALL,,,0,1


! If material uncertainty flag == 1
! =======================================                                                                                                                                      
*IF,UNCERT,EQ,1,THEN
	! Read the Material Properties Input file
	! =======================================                                                                                                                                      
	*DIM,BeamProps,ARRAY,NoExp,1 
	*VREAD,BeamProps(1),'BeamPropsVal','txt',,IJK,NoExp,1
(E16.7)                                                                                                                                                                
*ENDIF

! Read the Acceleration Input file                                                                                                                                      
! ================================
*DIM,EXCIT,ARRAY,NoExp*NoS,1
*CREATE,ansuitmp                                                                                                                    
*VREAD,EXCIT,'AccRandomVal','txt',,IJK,NoExp*NoS,1
(E16.7)                                                                                                                                                                
*END                                                                                                                                                                   
/INPUT,ansuitmp

FINISH


/DELETE,'MPropsVal',txt,
! =========================================================
! SOLUTION 
! =========================================================
*DO,J,1,NoExp			! Experiments DO LOOP
	
	STRin = STRCAT('MultiIIval_',CHRVAL(J))
	STRd = STRCAT('MultiIIval_dspl_',CHRVAL(J))
	STRv = STRCAT('MultiIIval_vel_',CHRVAL(J))
	STRa = STRCAT('MultiIIval_acc_',CHRVAL(J))

	! PREP7 
	! =====================================================
	*IF,UNCERT,EQ,1,THEN
	
		/PREP7
		MPDELE,ALL,1						! Delete material properties
		TBDELE,ALL,1						! Delete nonlinear material properties
	
		! Material 1: Vertical beams 
		MP,EX,1,200E9*BeamProps(J,1)        ! N/m^2	
		MP,PRXY,1,0.29
		MP,DENS,1,7850
		! Nonlinear material properties
		TB,BISO,1,1	 						! Redefine bilinear isotropic hardening  
		TBDATA,1,200E6,20E9
		KEYOPT,1,2,0						! No shear deflection
		KEYOPT,1,4,1						! Print out member forces and moments in the element coordinate system
		KEYOPT,1,6,0						! No shear deflection
		FINISH

		! Write material properties for the current experiment
		! =====================================================
		*CFOPEN,'MPropsVal',txt,,APPEND			! Open txt file
		*GET,MatEx,EX,1							! Material 1 Young moduli
		*GET,YieldStress,BISO,1,TEMP,0,CONST,1	! Material 1 yiels stress
		*GET,BEAMheight1,RCON,1,CONST,3			! Beam 1 cross section area
		*GET,BEAMheight2,RCON,2,CONST,3			! Beam 2 cross section area
		*GET,BEAMheight3,RCON,3,CONST,3			! Beam 1 cross section area
		*GET,BEAMheight4,RCON,4,CONST,3			! Beam 2 cross section area
		*GET,BEAMheight5,RCON,5,CONST,3			! Beam 1 cross section area
		*GET,BEAMheight6,RCON,6,CONST,3			! Beam 2 cross section area
		*VWRITE,J,MatEx,YieldStress,BEAMheight1,BEAMheight2,BEAMheight3,BEAMheight4,BEAMheight5,BEAMheight6
(9(E16.7,2X))
		*CFCLOS								! Close 'MaterialProps' txt file

	*ENDIF

	! =========================================================
	! SOLUTION 
	! =========================================================
	/SOLU
	SOLCONTROL,ON
	! Rayleigh Damping
	! ================
	ALPHAD,4.0
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
	NSEL,S,LOC,Y,0
	D,ALL,ALL				! Constraint all dofs for fixed supports
	ALLSEL,ALL

	! Initial conditions
	! ==================
	TIMINT,OFF,STRUC 		! Static load step (dead weight)
	KBC,1 
	NSUBST,100
	ACEL,0,9.81				! Acceleration of gravity 
	TIME,1e-6 
	SOLVE

	! Transient analysis
	! ==================
	TIMINT,ON				! Transient Analysis
	AUTOTS,ON				! Automatic time and/or load stepping

	*DO,I,1,NoS				! Time index DO LOOP
	
		! Time index
		Tindx = DT*(I)
		TIME,Tindx
		! X-axis input acceleration
		Gacc = EXCIT((J-1)*NoS+I)
		! Application of the X-axis & gravity acceleration 
		ACEL,Gacc,9.81

		KBC,0					! Ramped load steps
		NSUBST,10				! Specifies the number of substeps to be taken this load step 
		SOLVE
		
		*GET,Ux2,NODE,2,U,X 		! X displacement at node 2		
		*GET,Ux4,NODE,4,U,X 		! X displacement at node 4		
		*GET,Ux6,NODE,6,U,X 		! X displacement at node 6
		*GET,Ux8,NODE,8,U,X 		! X displacement at node 8		
		*GET,Ux9,NODE,9,U,X 		! X displacement at node 9		
		*GET,Ux10,NODE,10,U,X 		! X displacement at node 10		
		*GET,Ux11,NODE,11,U,X 		! X displacement at node 11		
		*GET,Ux12,NODE,12,U,X 		! X displacement at node 12		
		*GET,Ux13,NODE,13,U,X 		! X displacement at node 13		
		*GET,Ux14,NODE,14,U,X 		! X displacement at node 14		
		*GET,Ux15,NODE,15,U,X 		! X displacement at node 15		
		*GET,Ux16,NODE,16,U,X 		! X displacement at node 16		
		*GET,Ux17,NODE,17,U,X 		! X displacement at node 17		
		*GET,Ux18,NODE,18,U,X 		! X displacement at node 18		
		*GET,Ux19,NODE,19,U,X 		! X displacement at node 19		
		*GET,Ux20,NODE,20,U,X 		! X displacement at node 20		

		
		*GET,Vx2,NODE,2,V,X 		! X velocity at node 2		
		*GET,Vx4,NODE,4,V,X 		! X velocity at node 4		
		*GET,Vx6,NODE,6,V,X 		! X velocity at node 6
		*GET,Vx8,NODE,8,V,X 		! X velocity at node 8		
		*GET,Vx9,NODE,9,V,X 		! X velocity at node 9		
		*GET,Vx10,NODE,10,V,X 		! X velocity at node 10		
		*GET,Vx11,NODE,11,V,X 		! X velocity at node 11
		*GET,Vx12,NODE,12,V,X 		! X velocity at node 12	
		*GET,Vx13,NODE,13,V,X 		! X velocity at node 13	
		*GET,Vx14,NODE,14,V,X 		! X velocity at node 14	
		*GET,Vx15,NODE,15,V,X 		! X velocity at node 15		
		*GET,Vx16,NODE,16,V,X 		! X velocity at node 16
		*GET,Vx17,NODE,17,V,X 		! X velocity at node 17	
		*GET,Vx18,NODE,18,V,X 		! X velocity at node 18	
		*GET,Vx19,NODE,19,V,X 		! X velocity at node 19
		*GET,Vx20,NODE,20,V,X 		! X velocity at node 20	
		
		
		*GET,Ax2,NODE,2,A,X			! X acceleration at node 2
		*GET,Ax4,NODE,4,A,X			! X acceleration at node 4
		*GET,Ax6,NODE,6,A,X			! X acceleration at node 6
		*GET,Ax8,NODE,8,A,X			! X acceleration at node 8
		*GET,Ax9,NODE,9,A,X			! X acceleration at node 9
		*GET,Ax10,NODE,10,A,X		! X acceleration at node 10
		*GET,Ax11,NODE,11,A,X		! X acceleration at node 11
		*GET,Ax12,NODE,12,A,X		! X acceleration at node 12
		*GET,Ax13,NODE,13,A,X		! X acceleration at node 13
		*GET,Ax14,NODE,14,A,X		! X acceleration at node 14
		*GET,Ax15,NODE,15,A,X		! X acceleration at node 15
		*GET,Ax16,NODE,16,A,X		! X acceleration at node 16
		*GET,Ax17,NODE,17,A,X		! X acceleration at node 17
		*GET,Ax18,NODE,18,A,X		! X acceleration at node 18
		*GET,Ax19,NODE,19,A,X		! X acceleration at node 19
		*GET,Ax20,NODE,20,A,X		! X acceleration at node 20
	
		*GET,ReacFx1,NODE,1,RF,FX 	! X Reaction Force at node 1
		*GET,ReacFx3,NODE,3,RF,FX 	! X Reaction Force at node 3
		*GET,ReacFx5,NODE,5,RF,FX 	! X Reaction Force at node 5
		*GET,ReacFx7,NODE,7,RF,FX 	! X Reaction Force at node 7
	
		*CFOPEN,STRin,txt,,APPEND		! Open txt file
		*VWRITE,Tindx,Gacc,ReacFx1,ReacFx3,ReacFx5,ReacFx7,ReacFx1+ReacFx3+ReacFx5+ReacFx7
(7(E16.7,2X))
		*CFCLOS					
	
		*CFOPEN,STRd,txt,,APPEND		! Open txt file
		*VWRITE,Ux2,Ux4,Ux6,Ux8,Ux9,Ux10,Ux11,Ux12,Ux13,Ux14,Ux15,Ux16,Ux17,Ux18,Ux19,Ux20
(16(E16.7,2X))
		*CFCLOS
	
		*CFOPEN,STRv,txt,,APPEND		! Open txt file
		*VWRITE,Vx2,Vx4,Vx6,Vx8,Vx9,Vx10,Vx11,Vx12,Vx13,Vx14,Vx15,Vx16,Vx17,Vx18,Vx19,Vx20
(16(E16.7,2X))
		*CFCLOS					
		
		*CFOPEN,STRa,txt,,APPEND		! Open txt file
		*VWRITE,Ax2,Ax4,Ax6,Ax8,Ax9,Ax10,Ax11,Ax12,Ax13,Ax14,Ax15,Ax16,Ax17,Ax18,Ax19,Ax20
(16(E16.7,2X))
		*CFCLOS					
	*ENDDO
	FINISH
*ENDDO