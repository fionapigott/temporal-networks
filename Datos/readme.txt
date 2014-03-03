---------------------------------------------------------------
Fiona Pigott
Date: Jan 23 2014

Network data for cellular social network
Monthly weighted (by total minute duration of calls during a month)
directed (caller -> receiver) adjacency matrices and edge lists.

DL: original data, including monthly edge lists and internet 
	adoption data. Edge list files are either labeled 
	monthYear.csv (feb05.csv) or monthNumber.csv (2.csv)
	All files are ready to be used as "edge lists" in Gephi.

GephiMat: Adjacency matrices exported from Gephi
		Note that matrices include phone number labels
		as the first row and column, and are unordered.
		To order phone numbers from lowest to highest to
		match the ordering used in this project, import to
		MATLAB and use command:
		A = sortrows(sortrows(A)')';
		All files are stored as comma-separated .csv files
		named as YearMonth.csv (2005B.csv) where months are 
		letters (A, B, C, D, E, F, G, H, I , J, K, L) 
		corresponding to month number 
		(A = Jan, B = Feb ... L = Dec)
	DO NOT change the format of the data titles: the alphabetical order
	and name length are used to read in the data in MATLAB

	For the interested: to change a folder of semicolon-delimited files 
	to a folder of comma-delimited files, I used this grep/sed command 
	on Mac OSX terminal from the file where the data is stored:
		grep -lF ';' * | xargs sed -i "" 's/;/,/g'

	FEBRUARY: Because of an error in the data, months 2004D and 2007J 
	          have been intentionally removed from the folder
		  and added to “badData”

GephiProjects: Project files used to convert edge lists into the
			adjacency matrices. The command "File->Export->Graph" 
			will create a .csv file as found in GephiMat.
			Files are YearMonthProj.gephi (2005BProj.gephi).

The version of Gephi used to format this data is gephi-0.8.1-beta 
on MacOSX v 10.9.1.
Gephi is an open-source and multiplatform software for creating 
graphs distributed under the dual license CDDL 1.0 and 
GNU General Public License v3.
-------------------------------------------------------------------

