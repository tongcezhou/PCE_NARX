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
*SET,UNCERT,1							! Material uncertainty flag 

! Material 1: Vertical beams 
ET,1,BEAM188
MP,EX,1,200E9
MP,PRXY,1,0.29
MP,DENS,1,7850
KEYOPT,1,3,3
! Nonlinear material properties
TB,BISO,1								! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,250E6,10E9

! Cross-section 1 (vertical elements)
SECTYPE,1,BEAM,RECT
SECDATA,0.30,0.30,3,3
! Cross-section 2 (vertical elements)
SECTYPE,2,BEAM,RECT
SECDATA,0.25,0.25,3,3
! Cross-section 3 (vertical elements)
SECTYPE,3,BEAM,RECT
SECDATA,0.20,0.20,3,3
! Cross-section 4 (vertical elements)
SECTYPE,4,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 5 (vertical elements)
SECTYPE,5,BEAM,RECT
SECDATA,0.15,0.15,3,3
! Cross-section 6 (horizontal outer elements)
SECTYPE,6,BEAM,RECT
SECDATA,0.15,0.15,3,3
SECCONTROLS,,,,1631
! Cross-section 7 (horizontal inner elements)
SECTYPE,7,BEAM,RECT
SECDATA,0.15,0.15,3,3
SECCONTROLS,,,,3262


! Keypoints (front profile)
! Ground floor
K,1, 0.0, 0.0, 0.0
K,2, 4.0, 0.0, 0.0
K,3, 8.0, 0.0, 0.0
K,4, 0.0, 4.0, 0.0
K,5, 4.0, 4.0, 0.0
K,6, 8.0, 4.0, 0.0
K,7, 0.0, 8.0, 0.0
K,8, 4.0, 8.0, 0.0
K,9, 8.0, 8.0, 0.0
! 1st floor
K,10, 0.0, 0.0, 3.0
K,11, 4.0, 0.0, 3.0
K,12, 8.0, 0.0, 3.0
K,13, 0.0, 4.0, 3.0
K,14, 4.0, 4.0, 3.0
K,15, 8.0, 4.0, 3.0
K,16, 0.0, 8.0, 3.0
K,17, 4.0, 8.0, 3.0
K,18, 8.0, 8.0, 3.0
! 2nd floor
K,19, 0.0, 0.0, 6.0
K,20, 4.0, 0.0, 6.0
K,21, 8.0, 0.0, 6.0
K,22, 0.0, 4.0, 6.0
K,23, 4.0, 4.0, 6.0
K,24, 8.0, 4.0, 6.0
K,25, 0.0, 8.0, 6.0
K,26, 4.0, 8.0, 6.0
K,27, 8.0, 8.0, 6.0
! 3rd floor
K,28, 0.0, 0.0, 9.0
K,29, 4.0, 0.0, 9.0
K,30, 8.0, 0.0, 9.0
K,31, 0.0, 4.0, 9.0
K,32, 4.0, 4.0, 9.0
K,33, 8.0, 4.0, 9.0
K,34, 0.0, 8.0, 9.0
K,35, 4.0, 8.0, 9.0
K,36, 8.0, 8.0, 9.0
! 4th floor
K,37, 0.0, 0.0, 12.0
K,38, 4.0, 0.0, 12.0
K,39, 8.0, 0.0, 12.0
K,40, 0.0, 4.0, 12.0
K,41, 4.0, 4.0, 12.0
K,42, 8.0, 4.0, 12.0
K,43, 0.0, 8.0, 12.0
K,44, 4.0, 8.0, 12.0
K,45, 8.0, 8.0, 12.0
! 5th floor
K,46, 0.0, 0.0, 15.0
K,47, 4.0, 0.0, 15.0
K,48, 8.0, 0.0, 15.0
K,49, 0.0, 4.0, 15.0
K,50, 4.0, 4.0, 15.0
K,51, 8.0, 4.0, 15.0
K,52, 0.0, 8.0, 15.0
K,53, 4.0, 8.0, 15.0
K,54, 8.0, 8.0, 15.0

! Lines
! TYPE,1
! Vertical Beams b/w grounf-1st floor
TYPE,1
L,1,10	! Element 1
L,2,11	! Element 2
L,3,12	! Element 3
L,4,13	! Element 4
L,5,14	! Element 5
L,6,15	! Element 6
L,7,16	! Element 7
L,8,17	! Element 8
L,9,18	! Element 9
! Vertical Beams b/w 1st floor-2nd floor
L,10,19	! Element 10
L,11,20	! Element 11
L,12,21	! Element 12
L,13,22	! Element 13
L,14,23	! Element 14
L,15,24	! Element 15
L,16,25	! Element 16
L,17,26	! Element 17
L,18,27	! Element 18
! Vertical Beams b/w 2nd floor-3rd floor
L,19,28	! Element 19
L,20,29	! Element 20
L,21,30	! Element 21
L,22,31	! Element 22
L,23,32	! Element 23
L,24,33	! Element 24
L,25,34	! Element 25
L,26,35	! Element 26
L,27,36	! Element 27
! Vertical Beams b/w 3rd floor-4th floor
L,28,37	! Element 28
L,29,38	! Element 29
L,30,39	! Element 30
L,31,40	! Element 31
L,32,41	! Element 32
L,33,42	! Element 33
L,34,43	! Element 34
L,35,44	! Element 35
L,36,45	! Element 36
! Vertical Beams b/w 4th floor-5th floor
L,37,46	! Element 37
L,38,47	! Element 38
L,39,48	! Element 39
L,40,49	! Element 40
L,41,50	! Element 41
L,42,51	! Element 42
L,43,52	! Element 43
L,44,53	! Element 44
L,45,54	! Element 45

! TYPE,2
! Horizontal Beams 1st floor (outer beams)
L,10,11	! Element 46
L,11,12	! Element 47
L,10,13	! Element 48
L,12,15	! Element 49
L,13,16	! Element 50
L,15,18	! Element 51
L,16,17	! Element 52
L,17,18	! Element 53
! Horizontal Beams 2nd floor (outer beams)
L,19,20	! Element 54
L,20,21	! Element 55
L,19,22	! Element 56
L,21,24	! Element 57
L,22,25	! Element 58
L,24,27	! Element 59
L,25,26	! Element 60
L,26,27	! Element 61
! Horizontal Beams 3rd floor (outer beams)
L,28,29	! Element 62
L,29,30	! Element 63
L,28,31	! Element 64
L,30,33	! Element 65
L,31,34	! Element 66
L,33,36	! Element 67
L,34,35	! Element 68
L,35,36	! Element 69
! Horizontal Beams 4th floor (outer beams)
L,37,38	! Element 70
L,38,39	! Element 71
L,37,40	! Element 72
L,39,42	! Element 73
L,40,43	! Element 74
L,42,45	! Element 75
L,43,44	! Element 76
L,44,45	! Element 77
! Horizontal Beams 5th floor (outer beams)
L,46,47	! Element 78
L,47,48	! Element 79
L,46,49	! Element 80
L,48,51	! Element 81
L,49,52	! Element 82
L,51,54	! Element 83
L,52,53	! Element 84
L,53,54	! Element 85

! TYPE,3
! Horizontal Beams 1st floor (inner beams)
L,11,14	! Element 86
L,13,14	! Element 87
L,14,15	! Element 88
L,14,17	! Element 89
! Horizontal Beams 2nd floor (inner beams)
L,20,23	! Element 90
L,22,23	! Element 91
L,23,24	! Element 92
L,23,26	! Element 93
! Horizontal Beams 3rd floor (inner beams)
L,29,32	! Element 94
L,31,32	! Element 95
L,32,33	! Element 96
L,32,35	! Element 97
! Horizontal Beams 4th floor (inner beams)
L,38,41	! Element 98
L,40,41	! Element 99
L,41,42	! Element 100
L,41,44	! Element 101
! Horizontal Beams 5th floor (inner beams)
L,47,50	! Element 102
L,49,50	! Element 103
L,50,51	! Element 104
L,50,53	! Element 105

! Meshing 
LESIZE,ALL,,,1,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
! Verical beams (ground-1st floor) 
LATT,1,1,1,,,,1
LSEL,S,LINE,,1,9,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (1st-2nd floor) 
LATT,1,1,1,,,,2
LSEL,S,LINE,,10,18,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (2nd-3rd floor) 
LATT,1,1,1,,,,3
LSEL,S,LINE,,19,27,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (3rd-4th floor) 
LATT,1,1,1,,,,4
LSEL,S,LINE,,28,36,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Verical beams (4th-5th floor) 
LATT,1,1,1,,,,5
LSEL,S,LINE,,37,45,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Horizontal outer beams 
LATT,1,1,1,,,,6
LSEL,S,LINE,,46,85,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL  
! Horizontal inner beams 
LATT,1,1,1,,,,7
LSEL,S,LINE,,86,105,1,1
LMESH,ALL								! Generates nodes and line elements along lines
ALLSEL,ALL

! If material uncertainty flag == 1
! =======================================                                                                                                                                      
*IF,UNCERT,EQ,1,THEN
	! Read the Material Properties Input file
	! =======================================                                                                                                                                      
	*DIM,BeamProps,ARRAY,NoExp,3 
	*VREAD,BeamProps(1),'BeamProps','txt',,IJK,NoExp,3
(E16.7)                                                                                                                                                                
*ENDIF

! Read the Acceleration Input file                                                                                                                                      
! ================================
*DIM,EXCIT,ARRAY,NoExp*NoS,1
*CREATE,ansuitmp                                                                                                                    
*VREAD,EXCIT,'AccRandom','txt',,IJK,NoExp*NoS,1
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
	TBDATA,1,250e6*BeamProps(J,2),1E10
	KEYOPT,1,3,3
	
	! Cross-section 1 (vertical elements)
	SECTYPE,1,BEAM,RECT
	SECDATA,0.30*BeamProps(J,3),0.30*BeamProps(J,3),3,3
	! Cross-section 2 (vertical elements)
	SECTYPE,2,BEAM,RECT
	SECDATA,0.25*BeamProps(J,3),0.25*BeamProps(J,3),3,3
	! Cross-section 3 (vertical elements)
	SECTYPE,3,BEAM,RECT
	SECDATA,0.20*BeamProps(J,3),0.20*BeamProps(J,3),3,3
	! Cross-section 4 (vertical elements)
	SECTYPE,4,BEAM,RECT
	SECDATA,0.15*BeamProps(J,3),0.15*BeamProps(J,3),3,3
	! Cross-section 5 (vertical elements)
	SECTYPE,5,BEAM,RECT
	SECDATA,0.15*BeamProps(J,3),0.15*BeamProps(J,3),3,3

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
	*GET,BEAMheight7,SECP,7,DATA,1			! Beam 1 cross section area
	*VWRITE,J,MatEx,YieldStress,BEAMheight1,BEAMheight2,BEAMheight3,BEAMheight4,BEAMheight5,BEAMheight6,BEAMheight7
(10(E16.7,2X))
	*CFCLOS								! Close 'MaterialProps' txt file

	*ENDIF
	
	/SOLU
	SOLCONTROL,ON
	! Rayleigh Damping
	! ================
	ALPHAD,0.5
	BETAD,0.0001
	
	! Transient Analysis Options
	! ==========================
	ANTYPE,TRANS			! Specifies the analysis type and restart status.
	TRNOPT,FULL				! Specifies transient analysis options.
	NLGEOM,OFF				! Includes large-deflection effects in a static or full transient analysis.                                                                                                                                                           
	NROPT,FULL				! Specifies the Newton-Raphson options in a static or full transient analysis.                                                                                                                                                                                                                                                      
	EQSLV,SPARSE			! Specifies the type of equation solver.    
	AUTOTS,ON				! Automatic time and/or load stepping

	! Controls the solution data written to the database.                                                                                                                  
	OUTRES,ALL,NONE 
	OUTRES,NSOL,LAST
	OUTRES,NLOAD,LAST                                                                                                                                                       
	OUTRES,RSOL,LAST      
	OUTRES,A,LAST 
	
	! DOF constraints
	! ===============
	NSEL,S,LOC,Z,0
	D,ALL,ALL				! Constraint all dofs for fixed supports
	ALLSEL,ALL
	
	STR1 = STRCAT('SFres',CHRVAL(J))
	*CFOPEN,STR1,txt,		! Open txt file
	
	TIMINT,OFF,STRUC 		! Static load step (dead weight)
	KBC,1 
	NSUBST,5
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
		NSUBST,5				! Specifies the number of substeps to be taken this load step 
		SOLVE
	
		*GET,Ax18,NODE,18,A,X		! X acceleration at node 18
		*GET,Ax27,NODE,27,A,X		! X acceleration at node 27
		*GET,Ax36,NODE,36,A,X		! X acceleration at node 36
		*GET,Ax45,NODE,45,A,X		! X acceleration at node 45
		*GET,Ax54,NODE,54,A,X		! X acceleration at node 54
		*GET,Ux18,NODE,18,U,X 		! X displacement at node 18
		*GET,Ux27,NODE,27,U,X 		! X displacement at node 27		
		*GET,Ux36,NODE,36,U,X 		! X displacement at node 36		
		*GET,Ux45,NODE,45,U,X 		! X displacement at node 45		
		*GET,Ux54,NODE,54,U,X 		! X displacement at node 54		
		
		*GET,ReacFx1,NODE,1,RF,FX 	! X Reaction Force at node 1
		*GET,ReacFx2,NODE,3,RF,FX 	! X Reaction Force at node 2
		*GET,ReacFx3,NODE,5,RF,FX 	! X Reaction Force at node 3
		*GET,ReacFx4,NODE,7,RF,FX 	! X Reaction Force at node 4
		*GET,ReacFx5,NODE,9,RF,FX 	! X Reaction Force at node 5
		*GET,ReacFx6,NODE,11,RF,FX 	! X Reaction Force at node 6
		*GET,ReacFx7,NODE,13,RF,FX 	! X Reaction Force at node 7
		*GET,ReacFx8,NODE,15,RF,FX 	! X Reaction Force at node 8
		*GET,ReacFx9,NODE,17,RF,FX 	! X Reaction Force at node 9
		
		*VWRITE,Tindx,Gacc,ReacFx1+ReacFx2+ReacFx3+ReacFx4+ReacFx5+ReacFx6+ReacFx7+ReacFx8+ReacFx9,Ux18,Ux27,Ux36,Ux45,Ux54,Ax18+Gacc,Ax27+Gacc,Ax36+Gacc,Ax45+Gacc,Ax54+Gacc
(13(E16.7,2X))

	*ENDDO
	*CFCLOS					! Close 'SFresI' txt file

	FINISH
	
*ENDDO
