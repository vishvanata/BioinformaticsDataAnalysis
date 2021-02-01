# BioinformaticsDataAnalysis
This is the file directory which contains all the scripts using Python pandas and matplotlib to analyze variant pathogenicity data. 

In many files, I define functions which are then imported and called by other files. This is so that I don't have just one big file with a tons of lines. Instead, I've decoupled steps in my analysis pipeline according to the exact operation performed on a given input. 

Alot of the operations these functions perform is reading data from one file, calling some other functions which organize or run statistical tests on the data, and then call a final plotting function I define to output a chart or another ".txt" file with analytics. 

The files I am analyzing and visualizing are machine and deep learning Bioinformatics algorithms developed by researchers at the Center for the Study of Systems Biology at Georgia Tech. The algorithms I had used are ENTPRISE, ENTPRISE-X, and MEDICASCY. Please check "https://sites.gatech.edu/cssb/entprise/" for more information on these specific algorithms. 
