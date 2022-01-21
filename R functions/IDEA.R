IDEA <- function(NT, mu){

	if(length(NT)<2) {
		print("Warning: length of NT should be at least two.")
	}
	else{
	  NT 		<-	as.numeric(NT)
	  TT 		<-	length(NT)
	  s 		<- 	(1:TT)/mu
	  
	  y1		<-	log(NT)/s
	  y2		<-	s^2
	  y3		<-	log(NT)
#	  IDEA1		<-  cumsum(y2)*cumsum(y1)-cumsum(s)*cumsum(y3)
#	  IDEA2		<-	(1:TT)*cumsum(y2)-(cumsum(s))^2	
#	  IDEA		<-	exp(IDEA1/IDEA2)
#	  Rhat		<-	tail(IDEA,1)
	  IDEA1		<-  sum(y2)*sum(y1)-sum(s)*sum(y3)
	  IDEA2		<-	TT*sum(y2)-(sum(s))^2	
	  IDEA		<-	exp(IDEA1/IDEA2)
	  

  	  return(list(Rhat=IDEA, inputs=list(NT=NT, mu=mu)))
  	  }
}
