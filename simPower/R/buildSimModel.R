buildSimModel <-
function (EScond = .2, 
                           ESmod = .2,
                           ESint = .2,
                           bpath = c(.4,.3),
                           model = "growth") 
{
  
  modela1 <- paste0("mediator", " ~ " ,EScond,"*","condition" ,  collapse = " \n ") 
  modela2 <- paste0("mediator", " ~ " ,ESmod,"*","moderator" ,  collapse = " \n ") 
  modela3 <- paste0("mediator", " ~ " ,ESint,"*","interaction" ,  collapse = " \n ") 
  
  if (model == "growth") {
    modelb1 <- paste0("li", " ~ " ,bpath[1],"*","mediator" ,  collapse = " \n ") 
    modelb2 <- paste0("ls", " ~ " ,bpath[2],"*","mediator" ,  collapse = " \n ") 
    modelb <- paste0(modelb1,  "\n ", modelb2)
    
    modelc1 <- ("li =~ 1*y1 + 1*y2 + 1*y3 + 1*y4")
    modelc2 <- ("ls =~ 3*y1 + 2*y2 + 1*y3 + 0*y4")
  }
  else
  {
    modelb <- paste0("y1", " ~ " ,bpath[1],"*","mediator" ,  collapse = " \n ") 
    for (i in 2:length(bpath)){
      value <- paste0("y",i, " ~ " ,bpath[i],"*","mediator" ,  collapse = " \n ") 
      modelb <- paste0(modelb,  " ;  ", value," ; ", collapse = " \n ")
    }
    
    modelc1 <- "  "
    modelc2 <- "  "
    
  }
  
  modelind1 <- paste0("ind1 := ",EScond,"*",bpath[1],collapse = " \n ")
  modelind2 <- ifelse(model != "growth", "  ", paste0("ind2 := ",EScond,"*",bpath[2], collapse = " \n "))
  
  
  model <- paste0(modela1," \n ",modela2," \n ",modela3," \n ",
                  modelb," \n ",
                  modelc1," \n ", modelc2," \n ",
                  modelind1," \n ", modelind2 )
  return(model)
  
}
