Q1
Answer can be found in the ChemE 7770_Samavi_Prelim_1.pdf (pg 2-4)

Q2
Instuction for running the code
	Download "Q2.m" "model.m" "getparam.m" from Q2\~ in the same directory
	Run Q2.m from MATLAB
	-It will generate graphs for steady state protein production without inducer and also for Phase 1 and Phase 2 as 		
	described in the question.
	Additional
		Explanation and assumptions for the code is given in ChemE 7770_Samavi_Prelim_1.pdf (pg 5-6)
		U1, U2, U3 (ranks) are attached as Q2\U.txt 
		Time averaged sensitivity arrays are generated as the MATLAB variables S1,S2 and S3 for Phase 1 and early/late Phase 2.
Q3
Instuction for running the code
	Download "Flux.jl" "data.jl" and "solve.jl" files in the same directory.
	Set Julia working directory as the one containing the files.
	In the julia console, run 
			include("solve.jl")
	The code will generate a plot of protein level  vs inducer concentration on semi-log scale for a =n =2
	The code will generate a “flux.csv” file with optimized fluxes with varying levels of inducer. 
	Please see ChemE 7770_Samavi_Prelim_1.pdf (pg 7 and onwards)
