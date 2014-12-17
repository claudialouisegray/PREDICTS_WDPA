# why is distance to boundary changing?

# what are the studies that are in when fungi were included, that are not here

old.data <- read.csv("matched.landuse_11_2014_with fungi.csv")
new.data <- read.csv("matched.landuse_11_2014.csv")

new.SS <- unique(new.data$SS)
old.SS <- unique(old.data$SS)

#make names match
new.SS <- gsub("Plants", "Plants and Fungi", new.SS)

# which studies were dropped.
to.check <- old.SS[which(old.SS %in% new.SS == F)]
to.check <- gsub(" Plants and Fungi", "", to.check)


diversity <- readRDS('N:/Documents/PREDICTS/data from Lawrence/diversity-2014-11-13-03-40-23.rds')


dropped <- subset(diversity, SS %in% to.check) # these are all fungi. 
unique(dropped$Kingdom)



# any studies new?

to.check2 <- new.SS[which(new.SS %in% old.SS == F)]
#yes - these walker ones, that for some reason were not in the first nov cut. 
# why?
# only thing that has changed is merge sites, but that would only combine sites with the same coordinates?