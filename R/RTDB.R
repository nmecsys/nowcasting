<<<<<<< HEAD
# #' @title Create Real Time Data Base
# #' @description Create a time series matrix \code{mts} replicating the information available in a given date.
# #' @param series_code vector with the series encoding following the Bacen (Banco Central do Brasil) standards.
# #' @param vintage the vintage encoded by the day of the extraction
# #' @import DBI RMySQL xts zoo
# #' @importFrom stats ts
# #' @examples
# #' \dontrun{
# #' # Show series available:
# #' RTDB()
# #' 
# #' # Show vintages available for the series 1:
# #' RTDB(series_code = 1)
# #' 
# #' # Show series 1 data at vintage 2017-04-04:
# #' RTDB(series_code = 1, vintage = "2017-04-04")
# #' }
# #' @references Central Bank of Brazil
# #' @export
=======
# @title Create Real Time Data Base
# @description Create a time series matrix \code{mts} replicating the information available in a given date.
# @param series_code vector with the series encoding following the Bacen (Banco Central do Brasil) standards.
# @param vintage the vintage encoded by the day of the extraction
# @import DBI RMySQL xts zoo
# @importFrom stats ts
# @examples
# \dontrun{
# # Show series available:
# RTDB()
# 
# # Show vintages available for the series 1:
# RTDB(series_code = 1)
# 
# # Show series 1 data at vintage 2017-04-04:
# RTDB(series_code = 1, vintage = "2017-04-04")
# }
# @references Central Bank of Brazil
# @export
>>>>>>> 9a4dd10ce08de40411dbe3382053a487cd5c8451


RTDB <- function(series_code = NULL, vintage = NULL){

  # library(DBI)
  # library(RMySQL)
  # series_code<-22099
  # vintage<-as.Date('2017-01-01')
  # vintage<-'now'
  
  v_ind<-as.character(vintage)
  
  # if(v_ind == 'now'){
  # 
  # # Seleção da série de dados
  # # Abaixo a lista com os códigos das variáveis tanto no Bacen
  # codigos<-series_code
  # 
  # datas<-seq(as.Date("1994-07-01"),Sys.Date(),by="days") # vetor de datas desde 1994-07-01 até a data atual
  # base<-data.frame(datas)
  # start.time<-Sys.time()
  # for (i in 1:length(codigos)){
  #   serie <-{}
  #   serie_aux<-{}
  #   serie<-get.series.bacen(codigos[i])     # trago a série em formato de data.frame
  #   for (jdatas in 1:length(datas)){                # séries diárias, mensais e trimestrais no mesmo data.frame
  #     ind <- which(as.Date(serie$DF[,1],"%d/%m/%Y")==datas[jdatas])
  #     if (length(ind) == 0){
  #       serie_aux[jdatas] <- NA
  #     } else {
  #       serie_aux[jdatas]<- serie$DF[ind,2]
  #     }
  #   }
  #   base<-cbind(base,serie_aux)
  #   names(base)[i+1]<-paste0('serie',codigos[i])
  #   message(paste(i,'from',length(codigos),'series extracted'))
  # }
  # 
  # # Transformação da base para mensal ----
  # 
  # # Diária para mensal
  # # A média mensal representa a variável mensal
  # basexts <- xts(base[,-1],as.Date(base[,1]))
  # basemonth<-data.frame(apply.monthly(basexts,mean,na.rm=T))
  # basemonth[is.na(basemonth)]<-NA
  # 
  # year<-as.numeric(substr(row.names(basemonth)[1],1,4))
  # month<-as.numeric(substr(row.names(basemonth)[1],6,7))
  # mybase<-stats::ts(basemonth,start=c(year,month),freq=12)
  # 
  # return(mybase)
  # 
  # }else{

#### Não oferece nenhum argumento
  
  
  SQL<-paste("SELECT * FROM dbvintage")
  conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
  dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
  DBI::dbDisconnect(conn)
  temounaotem<-sapply(3:(dim(dados0)[2]-1),function(y) tapply(dados0[,y], dados0[,dim(dados0)[2]], function(x) sum(is.na(x))==length(x)))
  colnames(temounaotem)<-colnames(dados0)[3:(dim(dados0)[2]-1)]

    if(is.null(series_code) & is.null(vintage)){  # Caso não ofereça nenhum argumento
      SQL<-"SELECT `COLUMN_NAME` FROM `INFORMATION_SCHEMA`.`COLUMNS` WHERE `TABLE_SCHEMA`='pibnow' AND `TABLE_NAME`='dbvintage'"
      conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
      dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
      DBI::dbDisconnect(conn)
      dados1<-data.frame(series_code=substr(dados0[3:(dim(dados0)[1]-1),],6,50))
      return(dados1)
      
  #### Oferece apenas o argumento da vintage
      
    }else if(is.null(series_code)){ # Caso não ofereça o argumento series_code
      if(vintage>Sys.Date()){
        message('Are you interest in the nowcast? \nGreat! Try the function nowcast.')
      }else{
        SQL<-paste("SELECT * FROM dbvintage WHERE vintage_cod =",paste0("\'",v_ind,"\'"))
        conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
        dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
        DBI::dbDisconnect(conn)
        if(nrow(dados0)==0){ # Caso tal vintage não esteja disponível
          SQL<-paste("SELECT vintage_cod FROM dbvintage")
          conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
          dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
          DBI::dbDisconnect(conn)
          
          if(as.Date(vintage)<vintage_cod[1]){ # Se tal vintage é menor do que a primeira
            message('Sorry, this vintage is not available.',paste('\nTry vintage',vintage_cod[1]))
          }else{ # Caso essa vintage não seja a primeira
            vintage_cod<-as.Date(unique(dados0)[,1])
            ind_previous<-max(which(vintage_cod<as.Date(vintage)))
            message('Sorry, this vintage is not available.',paste('\nTry vintage',vintage_cod[ind_previous]),' or ',vintage_cod[ind_previous+1])
          }
        }else{ # Caso tal vintage esteja disponível
          
          if(as.Date(vintage)<as.Date(rownames(temounaotem)[dim(temounaotem)[1]])){ # se for uma vintage antes de 26/10/2017
            dados0<-colnames(temounaotem)[temounaotem[as.Date(rownames(temounaotem))==as.Date(vintage),]==FALSE]
            dados1<-data.frame(series_code=substr(dados0,6,50))
          }else{ # Se for uma vintage depois ou igual a 26/10/2017
            dados0<-colnames(temounaotem)
            dados1<-data.frame(series_code=substr(dados0,6,50))
          }
        return(dados1)
        }
      }
  #### Oferece apenas o argumento series_code
  
    }else if(is.null(vintage)){
      if(!(paste0('serie',series_code) %in% colnames(temounaotem))){   # Caso não tenha essa série
        message('Sorry, this serie(s) is(are) not available.')
      }else{
      dados0<-rownames(temounaotem)[temounaotem[,paste0('serie',series_code)]==FALSE]
      dados0<-c(as.Date(dados0),seq(as.Date(dados0[length(dados0)]),Sys.Date(),by = 'days')[-1])
      dados1<-data.frame(vintages=dados0)
      return(dados1)
      }
      
          
  #### Oferece os dois argumentos series_code e vintage    
      
    }else{
      if(vintage>Sys.Date()){
        message('Are you interest in the nowcast? \nGreat! Try the function nowcast.')
      }else{
        return_try<-tryCatch({
                    SQL<-paste("SELECT", paste0('X,',paste0("serie",series_code,collapse = ',')) ,"FROM dbvintage WHERE vintage_cod =",paste0("\'",v_ind,"\'"))
                    conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
                    dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
                    DBI::dbDisconnect(conn)},
                    error = function(err) {
                        return(FALSE)
          })
        if(return_try==FALSE){  # Se a série não está no banco de dados
        message('Sorry, this serie(s) is(are) not available.')
        }else{ # Se a série está no banco de dados
          if(nrow(dados0)==0){  # mas não há vintage disponível para ela
            # SQL<-paste("SELECT vintage_cod FROM dbvintage")
            # conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
            # dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
            # DBI::dbDisconnect(conn)
            # temounaotem
            # vintage_cod<-as.Date(rownames(temounaotem))
            vintage_cod<-as.Date(rownames(temounaotem)[!temounaotem[,colnames(temounaotem)==paste0('serie',series_code)]])
            if(as.Date(vintage)<vintage_cod[1]){  # Se a vintage é anterior à primeira vintage
              # v_ind<-vintage_cod[1]
              # SQL<-paste("SELECT", paste0('X,',paste0("serie",series_code,collapse = ',')) ,"FROM dbvintage WHERE vintage_cod =",paste0("\'",v_ind,"\'"))
              # conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
              # dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
              # DBI::dbDisconnect(conn)
              # dados1<-ts(dados0[,-1],start=as.numeric(c(substr(dados0[1,1],1,4),substr(dados0[1,1],6,7))),frequency=12)
              message('Sorry, this vintage is not available for this(ese) serie(s).',
                             paste('\nTry vintage',vintage_cod[1]))
              # invisible(dados1)
              
            }else{ # Se está entre uma das vintages.
              
              ind_previous<-max(which(vintage_cod<as.Date(vintage)))
              v_ind<-vintage_cod[ind_previous]
              SQL<-paste("SELECT", paste0('X,',paste0("serie",series_code,collapse = ',')) ,"FROM dbvintage WHERE vintage_cod =",paste0("\'",v_ind,"\'"))
              conn = dbConnect(MySQL(),db="pibnow",user="pibnow_user",password="123456",host="200.18.49.107",port=3306)
              dados0 <- DBI::dbGetQuery(conn = conn,statement = SQL)
              DBI::dbDisconnect(conn)
              dados1<-ts(dados0[,-1],start=as.numeric(c(substr(dados0[1,1],1,4),substr(dados0[1,1],6,7))),frequency=12)
              message('Sorry, this vintage is not available for this(ese) serie(s).',
                             paste('\nHowever, I return (invisible) the last vintage available for this serie:',vintage_cod[ind_previous],':)'))
              invisible(dados1)
              }
          }else{
              dados1<-ts(dados0[,-1],start=as.numeric(c(substr(dados0[1,1],1,4),substr(dados0[1,1],6,7))),frequency=12)
              return(dados1)
          }
        }
      }
    }
  }
# }








