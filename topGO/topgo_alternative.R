
#if(!isGeneric("GOFisher_alternative"))
setGeneric("GOFisher_greater", function(object) standardGeneric("GOFisher_greater"))
setGeneric("GOFisher_less", function(object) standardGeneric("GOFisher_less"))

setMethod("GOFisher_greater", c("classicCount"), 
          function(object) {
            
            contMat <- contTable(object)
            
            if(all(contMat == 0))
              p.value <- 1
            else
              p.value <- fisher.test(contMat, alternative = "greater")$p.value
            
            return(p.value)
          })

setMethod("GOFisher_less", c("classicCount"), 
          function(object) {
            
            contMat <- contTable(object)
            
            if(all(contMat == 0))
              p.value <- 1
            else
              p.value <- fisher.test(contMat, alternative = "less")$p.value
            
            return(p.value)
          })
