table(PA_11_14$Source_ID)
nrow(PA_11_14)

sort(100*table(PA_11_14$Source_ID)/nrow(PA_11_14))

  DL1_2009__Woinarski    AD1_2008__Billeter 
           5.27823858            7.05775901 

#check
100*nrow(subset(PA_11_14, Source_ID == "DL1_2009__Woinarski"))/nrow(PA_11_14)
