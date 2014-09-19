#!sh

#Copyright 2007 Kevin Thornton, University of California, Irvine.

#This script is released under the terms of the GNU Public License.

#This code is provided without warranty, as-is, and is
#simply an example of how do use the reg program,
#but the truth is there are better ways to go for more complex problems.
#I do not support modifications to this code (although you are of course
#free to do so, under the license terms).  Basically, that means I will ignore 
#email about this script, unless it points out a bona-fide bug.

#The exercise in this shell script is to generate some random data, 
#and then estimate theta=4Nu under the standard Wright-Fisher diploid model,
#and the infinitely-many sites mutation model,
#using S (the number of segregating sites in the sample) as our summary statistic

#step 1: generate 1000 data sets with theta ~U(1,50), and store the summary statistic
#in a file called "data"

R --no-save <<EOF
write(runif(1000,1,50),file="true_values",ncolumns=1)
EOF
ms 10 1000 -t tbs < true_values | grep segs | cut -d":" -f 2 > data

#step 2: Simulate from a uniform prior distribution on theta, which goes from 0 to 50.
#We will simulate 1 million draws from such a prior, and save the results in the file 
#"prior_temp"

R --no-save <<EOF
write(runif(1000000,0,50),file="prior_temp",ncolumns=1)
EOF

#step 3: simulate data under the model, using prior_temp to provide values of theta

ms 10 1000000 -t tbs < prior_temp | grep segs | cut -d":" -f 2 > prior_summaries

#step 4: generate a prior file of parameters + summaries

paste prior_temp prior_summaries > prior

#step 5: do the ABC step

../src/reg -p prior -d data -P 1 -S 1 -b output -T -t 0.001

#step 6: plot the distribution of the estimator
#This bit uses the locfit package in R to estimate the mode
#of each posterior distribution (the so-called "maximum a-posteriori" estimate,
#or MAP esstimate, of the parameter.  This is the value which maximises the 
#posterior probability of the parameters, given the data.
#The distribution of the MAP estimates are plotted, and a vertical line placed
#at the true value.  Hopefully, the peak of the distribution of the estimator
#is close to the vertical line.
R --no-save <<EOF
source("getmap.R")
estimates=array()
j=1
tv=scan("true_values")
for( i in 0:999 )
{
fn=paste("output.",i,".tangent.post.gz",sep="")
infile=gzfile(fn,"rb")
x=scan(fn,quiet=T)
estimates[j] = getmap(x)
j=j+1
}
write(estimates,file="estimates.txt",ncolumns=1)
estimates = estimates/tv
pdf("estimates.pdf",height=10,width=10,pointsize=24)
hist(estimates,xlab=expression(hat(theta)/theta),ylab="Density",main="",xlim=c(0,5),breaks="FD")
abline(v=1)
dev.off()
print(mean(estimates))
EOF
