library(data.table)
library(pbapply)

MontyHall <- function( N,CHANGE.DOOR=T,NUM.DOOR=3 ) {

 attempts.dt             <- data.table( winning.door=as.character(sample( x=seq(NUM.DOOR),size=N,replace=T )),chosen.door=as.character(sample( x=seq(NUM.DOOR),size=N,replace=T )) )
 attempts.dt[["result"]] <-
  pbapply(
   X      = attempts.dt,
   MARGIN = 1,
   FUN    = function(LINE) {
#=##>> works for classic; NUM.DOOR=3
#=#    MY.DOOR <- LINE["winning.door"] == LINE["chosen.door"]
#=#    return( c("WIN","LOST")[ as.integer(MY.DOOR==CHANGE.DOOR)+1 ] )
# ! check ! #cat( " >>>>>>>>>>\n",sep="" )
# ! check ! #cat( " > winning door : ",LINE["winning.door"],"\n",sep="" )
# ! check ! #cat( " > chosen door : ", LINE["chosen.door"], "\n",sep="" )
    DCLOSED <- as.character(seq(NUM.DOOR))
    DCLOSED <- DCLOSED[ DCLOSED != LINE["chosen.door"]  ]
    D.EMPTY <- DCLOSED[ DCLOSED != LINE["winning.door"] ]
# ! check ! #cat( " > set of closed doors : ",paste(DCLOSED,collapse=","),"\n",sep="" )
# ! check ! #print(DCLOSED)
# ! check ! #print(D.EMPTY)
    REMOVED <- sample( x=D.EMPTY,size=1,replace=F )
# ! check ! #cat( " > removed door : ",REMOVED,"\n",sep="" )
    DCLOSED <- DCLOSED[ DCLOSED != REMOVED ]
# ! check ! #cat( " > set of closed doors : ",paste(DCLOSED,collapse=","),"\n",sep="" )
# ! check ! #print(DCLOSED)
    if ( CHANGE.DOOR ) {
     MY.DOOR <- sample( x=DCLOSED,size=1,replace=F )
    } else             {
     MY.DOOR <- LINE["chosen.door"]
    }
# ! check ! #cat( " > do I change my door : ",c("no","yes")[ as.integer(CHANGE.DOOR)+1 ],"\n",sep="" )
# ! check ! #cat( " > my door : ",MY.DOOR, "\n",sep="" )
# ! check ! #cat( " > result : ",c("LOST","WIN")[ as.integer( MY.DOOR==LINE["winning.door"] )+1 ], "\n\n",sep="" )
    return( c("LOST","WIN")[ as.integer( MY.DOOR==LINE["winning.door"] )+1 ] )
   }
  )
 
 return(attempts.dt)

}

