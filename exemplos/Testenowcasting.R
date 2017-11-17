
# Importo a base de dados
database<-readRDS('exemplos/database.rds')

# Seleciono a base completa
base<-database$dados2

### Method 2sq
pib<-BRGDP[,8]
y<-month2qtr(diff(diff(pib,3),12))
x<-Bpanel(BRGDP[,-8],rep(4,dim(BRGDP)[2]),aggregate = T)
q<-1
r<-2
p<-1
now_2sq<-nowcast(y,x,q,r,p,method = '2sq')
now_2sq$main
nowcast.plot(now_2sq)

### Method 2sq
x1<-Bpanel(base,rep(1,dim(base)[2]),aggregate = T)
x2<-Bpanel(base,rep(2,dim(base)[2]),aggregate = T)
x3<-Bpanel(base,rep(3,dim(base)[2]),aggregate = T)
x4<-Bpanel(base,rep(4,dim(base)[2]),aggregate = T)
q<-1
r<-2
p<-1
now_2sq<-nowcast(y,x4,q,r,p,method = '2sq')
now_2sq$main
nowcast.plot(now_2sq)
nowcast.plot(now_2sq,type = 'factors')
nowcast.plot(now_2sq,type = 'eigenvalues')

### Method 2sm
x<-Bpanel(base,rep(3,dim(base)[2]),aggregate = F)
x1<-Bpanel(base,rep(4,dim(base)[2]),aggregate = F)
now_2sm<-nowcast(y,x1,q,r,p,method = '2sm')
now_2sm$main
nowcast.plot(now_2sm)
nowcast.plot(now_2sm,type = 'factors')
nowcast.plot(now_2sm,type = 'month_y')



i<-5
ts.plot(cbind(x[,i],x1[,i]),col=1:2)
cor(x[,i],x1[,i],use = "complete.obs")

### Method EM
now_em<-nowcast(y,x,q,r,p,'EM')
now_em$main
nowcast.plot(now_em)
