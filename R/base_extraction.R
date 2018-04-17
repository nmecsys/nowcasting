#' @title Create a real time data base
#' @description Create a time series matrix \code{mts} extracting information from Bacen (Banco Central do Brasil) API.
#' @param series_code Vector with the series encoding following the Bacen (Banco Central do Brasil) standards.
#' @import xts
#' @importFrom stats ts
#' @import zoo
#' @examples 
#' # Extracting GDP serie at real-time from Central Bank of Brasil data base
#' \dontrun{
#' gdp<-base_extraction(22099)
#' # Industrial production (21859) serie at real-time from Central Bank of Brasil data base
#' ind_prod<-base_extraction(21859)
#' 
#' # Creating real time data base with the series: 
#' # Vehicles production (1373);
#' # Industrial production, general index (21859).
#' mybase<-base_extraction(c(1373,21859))
#' 
#' # Creating real time data base with the series: 
#' # Exchange rate - Free - United States dollar (1);
#' #  Interest rate - CDI (12).
#' mybase<-base_extraction(c(1,12))
#' 
#' # Creating real time data base with the series: 
#' # Vehicles production (1373);
#' # Credit Sales Index (1453);
#' # Retail sales (1455);
#' # Industrial production, general index (21859).
#' mybase<-base_extraction(c(1373,1453,1455,21859))}
#' 
#' @references Central Bank of Brazil
#' 
#' @export


base_extraction<-function(series_code){

# Seleção da série de dados
# Abaixo a lista com os códigos das variáveis tanto no Bacen
codigos<-series_code

datas<-seq(as.Date("1994-07-01"),Sys.Date(),by="days") # vetor de datas desde 1994-07-01 até a data atual
base<-data.frame(datas)
start.time<-Sys.time()
for (i in 1:length(codigos)){
serie <-{}
serie_aux<-{}
serie<-get.series.bacen(codigos[i])     # trago a série em formato de data.frame
  for (jdatas in 1:length(datas)){                # séries diárias, mensais e trimestrais no mesmo data.frame
    ind <- which(as.Date(serie$DF[,1],"%d/%m/%Y")==datas[jdatas])
    if (length(ind) == 0){
      serie_aux[jdatas] <- NA
    } else {
      serie_aux[jdatas]<- serie$DF[ind,2]
    }
  }
  base<-cbind(base,serie_aux)
  names(base)[i+1]<-paste0('serie',codigos[i])
  message(paste(i,'from',length(codigos),'series extracted'))
}

# Transformação da base para mensal ----

# Diária para mensal
# A média mensal representa a variável mensal
basexts <- xts(base[,-1],as.Date(base[,1]))
basemonth<-data.frame(apply.monthly(basexts,mean,na.rm=T))
basemonth[is.na(basemonth)]<-NA

year<-as.numeric(substr(row.names(basemonth)[1],1,4))
month<-as.numeric(substr(row.names(basemonth)[1],6,7))
mybase<-stats::ts(basemonth,start=c(year,month),freq=12)

return(mybase)

}



