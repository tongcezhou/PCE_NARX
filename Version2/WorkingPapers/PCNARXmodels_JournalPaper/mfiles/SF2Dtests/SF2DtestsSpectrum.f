/CLEAR
!/CWD 'C:\Users\sminas\Dropbox\MATLAB\EEASHM\PreliminaryAnalysis\'
/CONFIG,NRES,10000

! =======================================================
! PREPROCESSOR
! =======================================================
/PREP7
! Material 1: Vertical beams 
ET,1,BEAM188
MP,EX,1,200E9
MP,PRXY,1,0.29
MP,DENS,1,7850

! Cross-section 1 (vertical elements)
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
SECCONTROLS,,,,6524

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
K,10, 4.0,12.0
! 5th floor
K,11, 0.0,15.0
K,12, 4.0,15.0

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
LESIZE,ALL,,,10,1,1  					! Specifies the divisions and spacing ratio on unmeshed lines (no division, no spacing ratio)
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

FINISH


! =======================================================
! SPECTRUM SOLUTION 
! =======================================================
/SOLU

! Rayleigh Damping
! ================
ALPHAD,4.265
BETAD,0.0006

! DOF constraints
! ===============
NSEL,S,LOC,Y,0
D,ALL,ALL				! Constraint all dofs for fixed supports
ALLSEL,ALL
	
ANTYPE,3			! Harmonic analysis

F,2,FX,10000

HARFRQ,0,10,    	! Frequency range
NSUBST,500, 		! Number of frequency steps
KBC,1				! Stepped loads

SOLVE
FINISH

/POST26
NSOL,2,2,U,X,UX2	! Get y-deflection data
PRVAR,2  			! Print data
PLVAR,2				! Plot data
FINISH