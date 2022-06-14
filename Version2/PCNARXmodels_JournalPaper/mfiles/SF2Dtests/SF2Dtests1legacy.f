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
! Nonlinear material properties
TB,BISO,2	 							! Bilinear isotropic hardening                                                                                                                      
TBDATA,1,200E9,10E11
KEYOPT,1,2,0							! No shear deflection
KEYOPT,1,4,1							! Print out member forces and moments in the element coordinate system
KEYOPT,1,6,0							! No shear deflection

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
NLIST,ALL,,,,NODE
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
ALPHAD,4.265
BETAD,0.0006
	
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

STRd = 'SF2Dtest1legacy_dspl'
STRv = 'SF2Dtest1legacy_vel'
STRa = 'SF2Dtest1legacy_acc'

TIMINT,OFF,STRUC 		! Static load step (dead weight)
KBC,1 
NSUBST,100
ACEL,0,9.81				! Acceleration of gravity 
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
		ACEL,Gacc,9.81

		KBC,0					! Ramped load steps
		NSUBST,100				! Specifies the number of substeps to be taken this load step 
		SOLVE


		*DO,ElemI,1,35				! Time index DO LOOP
		
			! ESEL,S,ELEM,,ElemI,ElemI,1
			STRx = STRCAT('SF2Dtest1legacy_Sx_',CHRVAL(ElemI))
			EPSex = STRCAT('SF2Dtest1legacy_Eex_',CHRVAL(ElemI))
			EPSpx = STRCAT('SF2Dtest1legacy_Epx_',CHRVAL(ElemI))
			
			! X stresses at element 1, node I at 8 integration points  
			*GET,SxI5,ELEM,ElemI,ETAB,
			
			ESEL,S,ELEM,ElemI,ElemI,1
			ETABLE,SIm50,LS,1
			ETABLE,SIm30,LS,2
			ETABLE,SI0,LS,3
			ETABLE,SIp30,LS,4
			ETABLE,SIp50,LS,5
			ETABLE,SMm50,LS,6
			ETABLE,SMm30,LS,7
			ETABLE,SM0,LS,8
			ETABLE,SMp30,LS,9
			ETABLE,SMp50,LS,10
			ETABLE,SJm50,LS,11
			ETABLE,SJm30,LS,12
			ETABLE,SJ0,LS,13
			ETABLE,SJp30,LS,14
			ETABLE,SJp50,LS,15
			
			ETABLE,EeIm50,LEPEL,1
			ETABLE,EeIm30,LEPEL,2
			ETABLE,EeI0,LEPEL,3
			ETABLE,EeIp30,LEPEL,4
			ETABLE,EeIp50,LEPEL,5
			ETABLE,EeMm50,LEPEL,6
			ETABLE,EeMm30,LEPEL,7
			ETABLE,EeM0,LEPEL,8
			ETABLE,EeMp30,LEPEL,9
			ETABLE,EeMp50,LEPEL,10
			ETABLE,EeJm50,LEPEL,11
			ETABLE,EeJm30,LEPEL,12
			ETABLE,EeJ0,LEPEL,13
			ETABLE,EeJp30,LEPEL,14
			ETABLE,EeJp50,LEPEL,15
			
			ETABLE,EpIm50,LEPPL,1
			ETABLE,EpIm30,LEPPL,2
			ETABLE,EpI0,LEPPL,3
			ETABLE,EpIp30,LEPPL,4
			ETABLE,EpIp50,LEPPL,5
			ETABLE,EpMm50,LEPPL,6
			ETABLE,EpMm30,LEPPL,7
			ETABLE,EpM0,LEPPL,8
			ETABLE,EpMp30,LEPPL,9
			ETABLE,EpMp50,LEPPL,10
			ETABLE,EpJm50,LEPPL,11
			ETABLE,EpJm30,LEPPL,12
			ETABLE,EpJ0,LEPPL,13
			ETABLE,EpJp30,LEPPL,14
			ETABLE,EpJp50,LEPPL,15
			
						
			*CFOPEN,STRx,txt,,APPEND		! Open txt file
			*VWRITE,Tindx,Gacc,SIm50,SIm30,SI0,SIp30,SIp50,SMm50,SMm30,SM0,SMp30,SMp50,SJm50,SJm30,SJ0,SJp30,SJp50
(17(E16.7,2X))

			*CFOPEN,EPSex,txt,,APPEND		! Open txt file
			*VWRITE,Tindx,Gacc,EeIm50,EeIm30,EeI0,EeIp30,EeIp50,EeMm50,EeMm30,EeM0,EeMp30,EeMp50,EeJm50,EeJm30,EeJ0,EeJp30,EeJp50
(17(E16.7,2X))

			*CFOPEN,EPSpx,txt,,APPEND		! Open txt file
			*VWRITE,Tindx,Gacc,EpIm50,EpIm30,EpI0,EpIp30,EpIp50,EpMm50,EpMm30,EpM0,EpMp30,EpMp50,EpJm50,EpJm30,EpJ0,EpJp30,EpJp50
(17(E16.7,2X))

		
		*ENDDO
		ESEL,ALL
		
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
*ENDDO