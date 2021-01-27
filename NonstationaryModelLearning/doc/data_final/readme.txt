Naming convention for some learning algorithms in the .csv files is different than the one used in the paper and in the julia code.
While in the julia code we updated the names in order to match the ones finally used in the paper, we did not do so for the .csv result files, in order to avoid the introduction of errors and in order to keep them "as is".

More specifically:

smileextended in csv files ---> VarSmile in paper and julia code
gsor1 in csv files ---> MP1 for Gaussian task in paper and MPN with 1 particle in julia code
gsor20 in csv files ---> MP20 for Gaussian task in paper and MPN with 20 particle in julia code
gsorOriginal1 in csv files ---> SOR1 for Gaussian task in paper and SOR with 1 particle in  julia code
gsorOriginal20 in csv files ---> SOR20 for Gaussian task in paper and SOR with 20 particle in julia code
tsor1 in csv files ---> MP1 for Categorical task in paper and MPN with 1 particle in  julia code
tsor20 in csv files ---> MP20 for Categorical task in paper and MPN with 20 particle in julia code
tsorOriginal1 in csv files ---> SOR1 for Categorical task in paper and SOR with 1 particle in julia code
tsorOriginal20 in csv files ---> SOR20 for Categorical task in paper and SOR with 20 particle in julia code
gnassarJN in csv files ---> Nas10* for Gaussian task in paper and Nas10 in julia code
gnassarNatN in csv files ---> Nas12* for Gaussian task in paper and Nas12 in julia code
