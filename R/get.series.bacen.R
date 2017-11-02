#' @import RCurl httr

get.series.bacen<- function(x, from = "", to = ""){


  if (missing(x)){
    stop("Need to specify at least one serie.")
  }
  
  if (! is.numeric(x)){
    stop("Argument x must be numeric.")
  }
  
  if (from == ""){
    data_init = "01/01/1980"
  } else {data_init = from}
  
  if (to == ""){
    data_end = format(Sys.Date(), "%d/%m/%Y")
  } else {data_end = to}
  
  inputs = as.character(x)
  len = seq_along(inputs)
  serie = mapply(paste0, "serie_", inputs, USE.NAMES = FALSE)
  
  
  for (i in len){ 
    
    result = tryCatch({
      
      
      RCurl::getURL(paste0('http://api.bcb.gov.br/dados/serie/bcdata.sgs.',
                           inputs[i], 
                           '/dados?formato=csv&dataInicial=', data_init, '&dataFinal=',
                           data_end),
                    ssl.verifyhost=FALSE, ssl.verifypeer=FALSE, .opts = list(timeout = 1, maxredirs = 2))
    },
    error = function(e) {
      
      return(RCurl::getURL(paste0('http://api.bcb.gov.br/dados/serie/bcdata.sgs.',
                                  inputs[i], 
                                  '/dados?formato=csv&dataInicial=', data_init,
                                  '&dataFinal=',
                                  data_end),
                           ssl.verifyhost=FALSE, ssl.verifypeer=FALSE))
      
    })
    
    assign(serie[i], result) 
  }
  
  
  sinal = tryCatch({
    for (i in len){
      texto = utils::read.csv2(textConnection(eval(as.symbol(
        serie[i]))), header=T)
      texto$data = gsub(' .*$','', eval(texto$data))
      assign(serie[i], texto)
      
    }},
    error = function(e){
      return("error")
    })
  
  
  
  if(sinal == "error"){
    for(i in len){
      texto=tryCatch({
        k = paste0('http://api.bcb.gov.br/dados/serie/bcdata.sgs.',
                   inputs[i], 
                   '/dados?formato=csv&dataInicial=', data_init, '&dataFinal=',
                   data_end)
        dados = httr::GET(k)
        aux = content(dados,'raw')
        aux2=rawToChar(aux)
        
        DF <- data.frame(do.call(cbind, strsplit(aux2, "\r\n", fixed=TRUE)))
        names(DF) <-"mist"
        DF$mist <- as.character(DF$mist)
        DF$mist<- gsub(x = DF$mist,pattern = '"',replacement = "")
        DF$data <- gsub(x = DF$mist,pattern = ";.*",replacement = "")
        DF$valor <- gsub(x = DF$mist,pattern = ".*;",replacement = "")
        DF$valor <- gsub(x = DF$valor,pattern = ",",replacement = ".")
        DF <- DF[-1,-1]
      })}
    assign(serie[i], result)
  }
  
  if(sinal == "error"){
    for (i in len){
      texto$data = gsub(' .*$','', eval(texto$data))
      assign(serie[i], texto)
    }
  }
  
  
  if(ncol(texto) == 1){
    for (i in len){
      texto = utils::read.csv(textConnection(eval(as.symbol(
        serie[i]))), header=T)
      texto$data = gsub(' .*$','', eval(texto$data))
      assign(serie[i], texto)
    }
  }
  
  rm(texto) 
  
  lista = list()
  ls_df = ls()[grepl('data.frame', sapply(ls(), function(x) class(get(x))))]
  for ( obj in ls_df ) { lista[obj]=list(get(obj)) }
  
  # return(invisible(lista))
  return(lista)
  
}