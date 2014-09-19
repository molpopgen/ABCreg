library(locfit);
getmap <- function(x)
{
	temp<-as.data.frame(x);
	d<-locfit(~temp[,1],temp);
	map<-temp[,1][which.max(predict(d,newdata=temp))]
	map;
}

process2plot<-function(x,relative=TRUE)
{
        temp <- sort(x);
        fit <- locfit(~temp);
        ffit <- fitted(fit);
        if(relative)
          {
            ffit <- ffit/max(ffit);
          }
        return(cbind(temp,ffit))
}
