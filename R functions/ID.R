ID <- function(NT, mu){

	  NT 	<-	as.numeric(NT)
	  TT 	<-	length(NT)
	  s 	<- 	(1:TT)/mu
	  y		<-	log(NT)/s
	  
	  R0_ID	<-	exp(sum(y)/TT)

  	  return(list=c(Rhat=R0_ID, inputs=list(NT=NT, mu=mu)))
}
