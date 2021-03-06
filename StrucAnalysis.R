# Dissertation Chapter 3 codes
# Structure in R


# http://membres-timc.imag.fr/Olivier.Francois/tutoRstructure.pdf
#############################################################################
# Installing packages

install.packages(c("fields","RColorBrewer","mapplots"))
source("http://bioconductor.org/biocLite.R")
biocLite("LEA")

install.packages("LEA_1.4.0_tar.gz", repos = NULL, type ="source")

source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

# So, to my understanding, this is a function that takes a stucture file, 
# and wrangles or reformats it somehow so it's usable in this package?
# Or maybe, input.file is just an excel file that is converted to the ./genotype.geno format
# I need to review the documentation for this package

struct2geno(file = input.file, TESS = FALSE, diploid = TRUE, FORMAT = 2,
extra.row = 0, extra.col = 0, output = "./genotype.geno")

input.file = "http://membres-timc.imag.fr/Olivier.Francois/secondary_contact.str"
struct2geno(file = input.file, TESS = TRUE, diploid = TRUE, FORMAT = 2,
extra.row = 0, extra.col = 0, output = "secondary_contact.geno")

library(LEA)
obj.snmf = snmf("secondary_contact.geno", K = 3, alpha = 100, project = "new")
qmatrix = Q(obj.snmf, K = 3)

barplot(t(qmatrix), col = c("orange","violet","lightgreen"), border = NA, space = 0,
xlab = "Individuals", ylab = "Admixture coefficients")

coord = read.table("coordinates.coord")
pop = rep(1:60, each = 10)

K = 3
Npop = length(unique(pop))
qpop = matrix(NA, ncol = K, nrow = Npop)
coord.pop = matrix(NA, ncol = 2, nrow = Npop)
for (i in unique(pop)){
qpop[i,] = apply(qmatrix[pop == i,], 2, mean)
coord.pop[i,] = apply(coord[pop == i,], 2, mean)}

library(mapplots)
plot(coord, xlab = "Longitude", ylab = "Latitude", type = "n")
map(add = T, col = "grey90", fill = TRUE)
for (i in 1:Npop){
add.pie(z = qpop[i,], x = coord.pop[i,1], y = coord.pop[i,2], labels = "",
col = c("orange","violet","lightgreen"))}
#maybe change the color scheme to viridis or magma

pop = scan("mypop.txt")

obj.snmf = snmf("secondary_contact.geno", K = 1:8, ploidy = 2, entropy = T,
alpha = 100, project = "new")
plot(obj.snmf, col = "blue4", cex = 1.4, pch = 19)

