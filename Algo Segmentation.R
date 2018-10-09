# ALGO D'AUTOMATISATION DU PROCESSUS DE SEGMENTATION

# rm(list=ls())

echantillons_appr_test<-function(matrice,taux_appr,ech_fixe=c("T","F")){
  if (taux_appr<=0 || taux_appr>=1){
    stop("taux_appr doit etre compris entre 0 et 1")
  } else {
    if (ech_fixe=="T"){
      set.seed(123)
      a<-runif(nrow(matrice)) # renvoie un entier entre 0 et 1
      d1<-(1 :nrow(matrice))[a<=taux_appr]
      d2<-(1 :nrow(matrice))[a>taux_appr]
      apprentissage<-matrice[d1,]
      test<-matrice[d2,]
      return(list(apprentissage=apprentissage,test=test))
      
    } else {
      a<-runif(nrow(matrice)) # renvoie un entier entre 0 et 1
      d1<-(1 :nrow(matrice))[a>taux_appr]
      d2<-(1 :nrow(matrice))[a<=taux_appr]
      apprentissage<-matrice[d1,]
      test<-matrice[d2,]
      return(list(apprentissage=apprentissage,test=test))
    }
  }
}

Err.rate<-function(observe,predit){
  mc<-table(observe,predit)
  error <- 1-sum(diag(mc))/sum(mc)
}

giga_correlation<-function(x,y){
  if (is.factor(x)==T){
    if (is.factor(y)==T) {
      # V de Cramer
      require(lsr)
      if(levels(x)>1 & levels(y)>1){
        v<-cramersV(x,y)
      } else {
        v<-NA
      }
      return(v)
    } else {
      # Rapport de correlation lineaire
      require(BioStatR)
      r<-sqrt(eta2(y,x))
      return(r)
    }
  } else {
    if (is.factor(y)==T) {
      # Rapport de correlation lineaire
      require(BioStatR)
      r<-sqrt(eta2(x,y))
      return(r)
    } else {
      # Coefficient de correlation lineaire
      c<-cor(x,y)
      return(c)
    }
  }
}

matrix_correlation<-function(x) {
  
  mat<-matrix(data=x,nrow=ncol(x),ncol = ncol(x))
  colnames(mat)<-colnames(x)
  rownames(mat)<-colnames(x)
  
  for (i in 1:ncol(x)){
    for (j in 1:ncol(x)){
      val<-giga_correlation(x[,i],x[,j])
      mat[i,j]<-val
    }
  }
  print(mat)
  mat[lower.tri(mat)] <- NA
  return(mat)
  print("**********************************************************************")
  if (ncol(x)>8) print("Essayer de reduire le tableau pour faciliter la lisibilite")
}

findBestClustering<-function(dfMCA,dfPCA,dfFAMD,dfCART=dfFAMD,nb_classes=seq(5,8),nb_facteurs=seq(2,5),
                             xySOM_preClust=c(7,7),
                             Factoriel=c('PCA','MCA','FAMD'),Methodes=c('Kmeans','CAH','Kohonen'),
                             hclust.method="ward.D2", taux_train=0.7){
  
  # Packages
  library(rpart)
  library(FactoMineR)
  library(ROCR)
  library(pROC)
  library(ade4)
  library(kohonen)
  library(rpart.plot)
  library(factoextra)
  library(sqldf)
  library(formattable)
  library(DT)
  library(nnet)
  library(beepr)
  library(lsr)
  library(BioStatR)
  
  debut<-proc.time()
  
  evolution<-length(nb_facteurs)*length(Factoriel)*length(Methodes)*length(nb_classes)
  
  Resultats<-data.frame(Prediction=factor(evolution),Methode=factor(evolution),Donnees=factor(evolution),
                        Nb_axes=numeric(evolution),k=numeric(evolution),CP=numeric(evolution),
                        Nb_split=numeric(evolution),Erreur_train=numeric(evolution),
                        Erreur_test=numeric(evolution),AUC_test=numeric(evolution),Axe1=numeric(evolution),
                        Axe2=numeric(evolution),Axe3=numeric(evolution),Axe4=numeric(evolution))
  levels(Resultats$Methode)<-c('Kmeans','CAH','Kohonen')
  levels(Resultats$Donnees)<-c('PCA','MCA','FAMD')
  levels(Resultats$Prediction)<-c('CART')
  
  # ANALYSE FACTO
  if("MCA" %in% Factoriel){
    facto1<-MCA(dfMCA,
                method="Burt",ncp=ncol(dfMCA),graph = F)
  }
  if("FAMD" %in% Factoriel){
    facto2<-FAMD(dfFAMD,
                 ncp=ncol(dfFAMD),graph = F)
  }
  if("PCA" %in% Factoriel){
    facto3<-PCA(dfPCA,
                ncp=ncol(dfPCA),graph = F)
  }
  
  p<-xySOM_preClust[1]
  q<-xySOM_preClust[2]
  
  # DEMARRAGE ALGO
  pb <- winProgressBar("test progress bar", "Some information in %",0, 100, 50)
  Sys.sleep(3)
  index<-0 
  
  for(i in 1:length(nb_facteurs)){ # Pour chaque nb de facteurs
    for(j in 1:length(Factoriel)){ # Parcourir les methodes factorielles specifiees
      
      if(Factoriel[j]=="MCA"){
        facto<-facto1
      }else if(Factoriel[j]=="FAMD"){
        facto<-facto2
      }else{
        facto<-facto3
      }
      
      # On applique le clustering sur les facteurs
      
      coord<-data.frame(facto$ind$coord)
      percent<-data.frame(Percent=round(facto$eig[,3][1:4],2))# vecteur des %
      for(k in 1:length(Methodes)){ # Sur les methodes de segmentation
        if(k==1){# Kmeans
          for(l in 1:length(nb_classes)){ # boucle sur les nb de classes specifies
            index<-index+1
            
            set.seed(123)
            classif<-kmeans(coord,centers = nb_classes[l])
            cluster<-classif$cluster
            new<-data.frame(dfCART,cluster=cluster)
            
            # Recherche du meilleur arbre
            ech<-echantillons_appr_test(new,taux_train,"T")
            train<-ech$apprentissage
            test<-ech$test
            set.seed(123)
            cart<-rpart(cluster~.,data=train,method="class",parms = list(split="gini"),cp=0)
            
            # CART: on choisit l'arbre elague au "CP+xstd" qui contient entre 3 et nb_classes+4 feuilles 
            # et qui minimise l'erreur sur l'echantillon de test
            cp<-cart$cptable
            cp<-subset(cp,cp[,'nsplit']>3 & cp[,'nsplit']<nb_classes[l]+5)
            auc<-matrix(NA,nrow(cp),4)
            for(g in 1:nrow(cp)){
              cartp<-prune(cart,cp=cp[g,"CP"]+cp[g,"xstd"])
              predc<-predict(cartp,type="class",test)
              
              auc[g,1]<-cp[g,'CP']+cp[g,'xstd']
              auc[g,2]<-cp[g,'nsplit']+1
              auc[g,4]<-Err.rate(test$cluster,predc)
              
              # Erreur en aprentissage
              predc2<-predict(cartp,type="class",train)
              
              auc[g,3]<-Err.rate(train$cluster,predc2)
            }
            colnames(auc)<-c('CP','nsplit','err_appr','err_test')
            
            
            # Meilleur modele CART
            top<-which.min(auc[,'err_test'])[1] # Le 1er top a forcement moins de feuilles
            
            # Remplissage du tableau
            Resultats[index,1]<-'CART'
            Resultats[index,2]<-Methodes[k]
            Resultats[index,3]<-Factoriel[j]
            Resultats[index,4]<-nb_facteurs[i]
            Resultats[index,5]<-nb_classes[l] # nb classe
            Resultats[index,6]<-round(auc[top,1],8)
            Resultats[index,7]<-auc[top,2]
            Resultats[index,8]<-round(auc[top,3],4)
            Resultats[index,9]<-round(auc[top,4],4)
            Resultats[index,11]<-percent[1,1]
            Resultats[index,12]<-percent[2,1]
            Resultats[index,13]<-percent[3,1]
            Resultats[index,14]<-percent[4,1]
            print(Resultats[index,])
            
            info <- sprintf("%d%% done", round(((index)/evolution)*100))
            setWinProgressBar(pb, round((index/evolution)*100), sprintf("test (%s)", info), info)
          }
        } else if (k==2){# CAH
          for(l in 1:length(nb_classes)){
            index<-index+1
            
            # Pre-classification
            set.seed(123)
            carte <- som(scale(coord),grid=somgrid(p,q,"hexagonal")) # Carte topographique de Kohonen
            pre_cluster<-carte$unit.classif
            
            # Matrice des distances entre les noeuds (codebooks)
            dc <- dist(carte$codes[[1]])
            
            # Demarrage CAH
            set.seed(123)
            cah <- hclust(dc,method=hclust.method)
            #plot(cah)
            groupes <- cutree(cah,k=nb_classes[l]) # groupes de noeuds
            cluster <- groupes[pre_cluster] # cluster des indiv
            
            new<-data.frame(dfCART,cluster=cluster)
            
            # Recherche du meilleur arbre
            ech<-echantillons_appr_test(new,taux_train,"T")
            train<-ech$apprentissage
            test<-ech$test
            set.seed(123)
            cart<-rpart(cluster~.,data=train,method="class",parms = list(split="gini"),cp=0)
            
            # CART: on choisit l'arbre elague au "CP+xstd" qui contient entre 3 et nb_classes+4 feuilles 
            # et qui minimise l'erreur sur l'echantillon de test
            cp<-cart$cptable
            cp<-subset(cp,cp[,'nsplit']>3 & cp[,'nsplit']<nb_classes[l]+5)
            auc<-matrix(NA,nrow(cp),4)
            for(g in 1:nrow(cp)){
              cartp<-prune(cart,cp=cp[g,"CP"]+cp[g,"xstd"])
              # Elagage au min(xerror)+std (1SE) comme le suggere Breiman mais xerror+1SE
              predc<-predict(cartp,type="class",test)
              
              auc[g,1]<-cp[g,'CP']+cp[g,'xstd']
              auc[g,2]<-cp[g,'nsplit']+1
              auc[g,4]<-Err.rate(test$cluster,predc)
              
              # Erreur en aprentissage
              predc2<-predict(cartp,type="class",train)
              
              auc[g,3]<-Err.rate(train$cluster,predc2)
            }
            colnames(auc)<-c('CP','nsplit','err_appr','err_test')
            
            
            # Meilleur modele CART
            top<-which.min(auc[,'err_test'])[1] # Le 1er top a forcement moins de feuilles
            
            # Remplissage du tableau
            Resultats[index,1]<-'CART'
            Resultats[index,2]<-Methodes[k]
            Resultats[index,3]<-Factoriel[j]
            Resultats[index,4]<-nb_facteurs[i]
            Resultats[index,5]<-nb_classes[l] # nb classe
            Resultats[index,6]<-round(auc[top,1],8)
            Resultats[index,7]<-auc[top,2]
            Resultats[index,8]<-round(auc[top,3],4)
            Resultats[index,9]<-round(auc[top,4],4)
            Resultats[index,11]<-percent[1,1]
            Resultats[index,12]<-percent[2,1]
            Resultats[index,13]<-percent[3,1]
            Resultats[index,14]<-percent[4,1]
            print(Resultats[index,])
            
            info <- sprintf("%d%% done", round(((index)/evolution)*100))
            setWinProgressBar(pb, round((index/evolution)*100), sprintf("test (%s)", info), info)
          }
        } else { # Kohonen
          for(l in 1:length(nb_classes)){
            
            if(nb_classes[l]%%2==0 & nb_classes[l]<=16){ # On parcourt seulement les nb de classes pairs pour cette methode
              index<-index+1
              
              if(nb_classes[l]==4){x<-2;y<-2}
              if(nb_classes[l]==6){x<-3;y<-2}
              if(nb_classes[l]==8){x<-4;y<-2}
              if(nb_classes[l]==10){x<-5;y<-2}
              if(nb_classes[l]==12){x<-4;y<-3}
              if(nb_classes[l]==14){x<-2;y<-7}
              if(nb_classes[l]==16){x<-4;y<-4}
              
              set.seed(123)
              carte <- som(scale(coord),grid=somgrid(x,y,"hexagonal")) # Carte topographique de Kohonen
              cluster<-carte$unit.classif
              
              new<-data.frame(dfCART,cluster=cluster)
              
              # Recherche du meilleur arbre
              ech<-echantillons_appr_test(new,taux_train,"T")
              train<-ech$apprentissage
              test<-ech$test
              set.seed(123)
              cart<-rpart(cluster~.,data=train,method="class",parms = list(split="gini"),cp=0)
              
              # CART: on choisit l'arbre elague au "CP+xstd" qui contient entre 3 et nb_classes+4 feuilles 
              # et qui minimise l'erreur sur l'echantillon de test
              cp<-cart$cptable
              cp<-subset(cp,cp[,'nsplit']>3 & cp[,'nsplit']<nb_classes[l]+5)
              auc<-matrix(NA,nrow(cp),4)
              for(g in 1:nrow(cp)){
                cartp<-prune(cart,cp=cp[g,"CP"]+cp[g,"xstd"])
                # Elagage au min(xerror)+std (1SE) comme le suggere Breiman mais xerror+1SE
                predc<-predict(cartp,type="class",test)
                
                auc[g,1]<-cp[g,'CP']+cp[g,'xstd']
                auc[g,2]<-cp[g,'nsplit']+1
                auc[g,4]<-Err.rate(test$cluster,predc)
                
                # Erreur en aprentissage
                predc2<-predict(cartp,type="class",train)
                
                auc[g,3]<-Err.rate(train$cluster,predc2)
              }
              colnames(auc)<-c('CP','nsplit','err_appr','err_test')
              
              
              # Meilleur modele CART
              top<-which.min(auc[,'err_test'])[1] # Le 1er top a forcement moins de feuilles
              
              # Remplissage du tableau
              Resultats[index,1]<-'CART'
              Resultats[index,2]<-Methodes[k]
              Resultats[index,3]<-Factoriel[j]
              Resultats[index,4]<-nb_facteurs[i]
              Resultats[index,5]<-nb_classes[l] # nb classe
              Resultats[index,6]<-round(auc[top,1],8)
              Resultats[index,7]<-auc[top,2]
              Resultats[index,8]<-round(auc[top,3],4)
              Resultats[index,9]<-round(auc[top,4],4)
              Resultats[index,11]<-percent[1,1]
              Resultats[index,12]<-percent[2,1]
              Resultats[index,13]<-percent[3,1]
              Resultats[index,14]<-percent[4,1]
              print(Resultats[index,])
              
              info <- sprintf("%d%% done", round(((index)/evolution)*100))
              setWinProgressBar(pb, round((index/evolution)*100), sprintf("test (%s)", info), info)
            }
          }
        }
      }
      
    }
  }
  colnames(Resultats)<-c('Prediction','Methode','Donnees','Nb_axes','k','CP','Nb_split','Erreur_train','Erreur_test','AUC',
                         'Axe1','Axe2','Axe3','Axe4')
  Sys.sleep(2)
  close(pb)
  library(beepr)
  beep('fanfare')
  print(time2<-proc.time()-debut)
  
  Resultats<-subset(Resultats,k!=0)
  print(datatable(formattable(Resultats)))
  return(Resultats)
}

caracterisation<-function(tableForCaract,clusterName="cluster"){
  if(clusterName!="cluster" || !clusterName%in%colnames(tableForCaract)){
    stop("clusterName must be equal to 'cluster' and tableForCaract must contains a column named 'cluster'")
  }else{
    listQuanti<-c()
    listQuali<-c()
    for(i in 1:ncol(tableForCaract)){
      if(is.numeric(tableForCaract[,i])){
        listQuanti<-c(listQuanti,colnames(tableForCaract)[i])
      }else{
        listQuali<-c(listQuali,colnames(tableForCaract)[i])
      }
    }
    listQuali<-setdiff(listQuali,clusterName)
    
    # Summary sur var quanti
    
    caract1<-aggregate(x=tableForCaract[,colnames(tableForCaract)%in%c(setdiff(listQuanti,clusterName))],
                       by=list(tableForCaract[,c(clusterName)]),FUN=mean,na.rm=T)
    caract1$nb<-table(data.frame(tableForCaract[,c(clusterName)]))
    
    # summary sur var quali
    
    if(!is.null(listQuali)){
      mytemp<-subset(tableForCaract,select = listQuali)
      varTocreate<-c()
      levelsToParsed<-c()
      corresp<-data.frame()
      for(i in 1:ncol(mytemp)){
        if(colnames(mytemp)[i]!='Id'){
          if(length(levels(mytemp[,i]))<=10){
            varTocreate<-c(varTocreate,paste(colnames(mytemp)[i],levels(mytemp[,i]),sep = "__"))
            levelsToParsed<-c(levelsToParsed,levels(mytemp[,i]))
            corresp<-rbind(corresp,data.frame(modalite=levels(mytemp[,i]),myVar=colnames(mytemp)[i]))
          }else{
            warning(cat("La variable ",colnames(mytemp[,i]),"sera ignoree dans la caracterisation"))
            warning("Car plus de 10 modalites")
          }
        }
      }
      varTocreate<-setdiff(varTocreate,"Id_")
      # cat("")
      # print(levelsToParsed)
      caract2<-data.frame(Cluster=unique(tableForCaract[,c(clusterName)]))
      for(i in 1:length(varTocreate)){
        df<-data.frame()
        for(j in 1:nrow(caract2)){
          r<-subset(tableForCaract,cluster==caract2$Cluster[j])
          mySelectedVar<-strsplit(varTocreate[i],split = "__")[[1]][1]
          
          indexCol<-which(colnames(r)==mySelectedVar)
          #print(levelsToParsed[i])
          #print(head(r[r[,indexCol]==levelsToParsed[i],]))
          r<-nrow(r[r[,mySelectedVar]==levelsToParsed[i],])
          df<-rbind(df,r)
        }
        caract2[,c(varTocreate[i])]<-df
      }
      #caract2<-caract2[order(caract2$Cluster,decreasing = F),]
      caract2$totalLigne<-caract1$nb
      for(i in 1:nrow(caract2)){
        for (j in 2:ncol(caract2)-1){
          caract2[i,j]<-round(caract2[i,j]/caract2[i,ncol(caract2)],2)
        }
      }
      caract2$Cluster<-caract1[,1]
    }else{
      caract2<-NULL
    }
    
    caract1[,2:ncol(caract1)]<-round(caract1[,2:ncol(caract1)],2)
    return(list(Quanti=caract1,Quali=caract2))
  }
}


FitBestClustering<-function(dfMCA,dfPCA,dfFAMD,dfCART=dfFAMD,
                            methode=c('Kmeans','CAH','Kohonen'),
                            source=c('MCA','FAMD','PCA'),
                            nbFacteurs,nbClasses,CP=NA,xySOM_preClust=c(7,7),
                            taux_train,prp=F,Corr=F,cluster.name='',percentInName=F,tableForCaract=NULL){

  library(DT)
  library(formattable)
  prediction=c('CART')
  name1<-paste("S",seq(1:nbClasses),sep = "")
  
  # Analyse factorielle
  if(source=='MCA'){
    print('MCA')
    facto<-MCA(dfMCA,method="Burt",ncp=as.numeric(nbFacteurs),graph = F)
    coord<-data.frame(facto$ind$coord)
  } else if(source=='FAMD'){
    print('FAMD')
    facto<-FAMD(dfFAMD,ncp=as.numeric(nbFacteurs),graph = F)
    coord<-data.frame(facto$ind$coord)
  } else {
    print('PCA')
    facto<-PCA(dfPCA,ncp=as.numeric(nbFacteurs),graph = F)
    coord<-data.frame(facto$ind$coord)
  }
  
  
  if(methode=='Kmeans'){
    print('Kmeans')
    set.seed(123)
    classif<-kmeans(coord,centers = as.numeric(nbClasses))
    cluster<-classif$cluster
    new<-data.frame(dfCART,cluster=cluster)
    Method_classif<-classif
  } else if(methode=='CAH'){
    print('CAH')
    set.seed(123)
    carte <- som(scale(coord),grid=somgrid(xySOM_preClust[1],xySOM_preClust[2],"hexagonal")) # Carte topographique de Kohonen
    pre_cluster<-carte$unit.classif
    # Matrice des distances entre les noeuds (codebooks)
    dc <- dist(carte$codes[[1]])
    set.seed(123)
    cah <- hclust(dc,method="ward.D2")
    groupes <- cutree(cah,k=as.numeric(nbClasses)) # groupes de noeuds
    cluster <- groupes[carte$unit.classif] # cluster des indiv
    new<-data.frame(dfCART,cluster=cluster)
    Method_classif<-cah
  } else { # kohonen
    print('Kohonen')
    set.seed(123)
    if(nbClasses%%2==0){ # On parcourt seulement les nb de classes pairs soit 5 classif
      if(nbClasses==4){x<-2;y<-2}
      if(nbClasses==6){x<-3;y<-2}
      if(nbClasses==8){x<-4;y<-2}
      if(nbClasses==10){x<-5;y<-2}
      if(nbClasses==12){x<-4;y<-3}
      if(nbClasses==14){x<-2;y<-7}
      if(nbClasses==16){x<-4;y<-4}
      
      carte <- som(scale(coord),grid=somgrid(x,y,"hexagonal"))
      cluster<-carte$unit.classif
      new<-data.frame(dfCART,cluster=cluster)
      Method_classif<-carte
    }
  }
  # Prediction
  set.seed(123)
  ech<-echantillons_appr_test(new,taux_train,"T")
  train<-ech$apprentissage
  test<-ech$test
  
  # Renommage des noms des clusters en supprimant les percent
  if(length(cluster.name)>1){
    if(percentInName==T){
      for(i in 1:length(unique(cluster))){
        l<-nchar(cluster.name[i])
        cluster<-gsub(i,str_sub(cluster.name[i],1,l-3),cluster)
      }
    }
  }
  
  # Decision Tree
  
  print('CART')
  set.seed(123)
  CART<-rpart(cluster~.,data=train,method="class",parms = list(split="gini"),cp=0)
  CART<-prune(CART,cp=CP)
  print('Fin CART')
  modele<-CART
  predc1<-predict(CART,type="class",train)
  train$Predicted<-predict(CART,type="class",train)
  predc2<-predict(CART,type="class",test)
  test$Predicted<-predict(CART,type="class",test)
  e1<-Err.rate(train$cluster,predc1)
  e2<-Err.rate(test$cluster,predc2)
  confusion_test<-table(test$cluster,predc2)
  Importance<-data.frame(Importance_var=round(CART$variable.importance/sum(CART$variable.importance)*100))
  
  print('Fin prediction')
  
  # CARACTERISATION DES CLUSTERS
  if(is.null(tableForCaract)){
    tableForCaract<-cbind(dfCART,cluster)
  }
  if(!c("cluster")%in%colnames(tableForCaract)){
    tableForCaract<-cbind(tableForCaract,cluster)
  }
  
  c<-caracterisation(tableForCaract,"cluster")
  caract1<-c$Quanti
  caract2<-c$Quali
  print("Fin caract")
  
  # Correlation
  if(Corr==T){
    print("1")
    dfMCA$cluster<-as.factor(cluster)
    print("2")
    dfFAMD$cluster<-as.factor(cluster)
    print("3")
    dfPCA$cluster<-as.factor(cluster)
    print("4")
    if(source=='ACM'){
      print("5")
      corr<-suppressWarnings(matrix_correlation(dfMCA))
    } else if(source=='AFDM'){
      print("6")
      corr<-suppressWarnings(matrix_correlation(dfFAMD))
    }else {
      print("7")
      corr<-suppressWarnings(matrix_correlation(dfPCA))
    }
  } else{
    corr<-NULL
  }
  
  # Graphique Decision Tree
  
  if(prediction=='CART' & prp==T){
    x11()
    # Renommer les cluster dans l'affichage graphique
    if(length(cluster.name)>1){
      attr(CART, "ylevels")<-cluster.name
    }
    
    prp(CART,type=0,extra=4,nn=T,branch=0.7,varlen=0,faclen=35,split.box.col="gray70",
        box.palette =list('antiquewhite4','hotpink1','darkolivegreen2','mediumpurple1','gray82',
                          'orange','antiquewhite1','grey30','yellow2','deepskyblue','coral3',
                          'firebrick1','darkslategray3','greenyellow','khaki4','tan3'),cex=0.72,
        main=paste("Arbre de decision: \t","Taux d'erreur",round(e2*100,0),"%"))
    
    cat('---------------------------------------------------------------\n\n')
    print(Importance)
    cat('---------------------------------------------------------------\n\n')
    print(CART)
    cat('---------------------------------------------------------------\n\n')
  }
  
  library(DT)
  library(formattable)
  
  print(datatable(formattable(caract1)))
  # print(datatable(formattable(caract2)))
  
  library(beepr)
  beep("fanfare")
  
  return(list(Factorielle=facto,Methode_classif=Method_classif,
              Prediction=list(Modele=modele,DataBases=list(train,test),Importance=Importance),
              Evaluation=list(confusion_test=confusion_test,error_train=e1,error_test=e2),
              Caracterisation=list(data=tableForCaract,Caract=list(Quanti=caract1),Corr=corr)))
}

exploreClustering<-function(result,dfMCA,dfPCA,dfFAMD,dfCART=dfFAMD,taux_train,xySOM_preClust=c(7,7)){
  # Identifie les segmentations pour lesquelles il n'est pas possible de predire certaines classes
  
  library(DT)
  library(formattable)
  
  analyse<-result[,c(1:9)]
  analyse$Erreur_train2<-rep(NA,nrow(analyse))
  analyse$Erreur_test2<-rep(NA,nrow(analyse))
  analyse$Predict_all<-rep(NA,nrow(analyse))
  
  time1<-proc.time()
  pb <- winProgressBar("test progress bar", "Some information in %",
                       0, 100, 50)
  Sys.sleep(0.5)
  j<-0
  for(i in 1:nrow(analyse)){
    j<-j+1
    Sys.sleep(0.0)
    
    param<-analyse[i,]
    #print(param)
    eval<-FitBestClustering(dfMCA,dfPCA,dfFAMD,dfCART,
                            methode=param[,2],source = param[,3],nbFacteurs = param[,4],
                            nbClasses =param[,5],CP=param[,6],xySOM_preClust =xySOM_preClust,
                            taux_train = taux_train)
    
    analyse$Erreur_train2[i]<-eval$Evaluation$error_train
    analyse$Erreur_test2[i]<-eval$Evaluation$error_test
    # Les classes non predites par le classifieur
    analyse$Predict_all[i]<-ifelse(length(setdiff(c(seq(1:analyse$k[i])),
                                                  levels(eval$Prediction$DataBases[[1]]$Predicted)))>0,
                                   as.character(setdiff(c(seq(1:analyse$k[i])),
                                                        levels(eval$Prediction$DataBases[[1]]$Predicted))),'Full')
    cat('\n',analyse$k[i],eval$Evaluation$error_test,analyse$Predict_all[i],'\n\n')
    if(i==nrow(analyse)){print(analyse)}
    
    info <- sprintf("%d%% done", round((j/(nrow(analyse)))*100))
    setWinProgressBar(pb, round((j/(nrow(analyse)))*100), sprintf("Progress (%s)", info), info)
  }
  Sys.sleep(5)
  close(pb)
  library(beepr)
  beep("fanfare")
  proc.time()-time1
  
  print(datatable(formattable(analyse)))
  
  return(analyse)
}
taux_val_manq<-function(x){
  t<-round((length(is.na(x)[is.na(x)==TRUE]))/length(x)*100,7)
  return(t)
}

taux_val_manq_matrix<-function(x){
  Missing<-as.data.frame(apply(x,2,function(y) sum(is.na(y))))
  Taux<-as.data.frame(apply(x,2,function(y) round(sum(is.na(y))/nrow(x)*100,4)))
  retour<-data.frame(Missing,Taux)
  colnames(retour)<-c("Missing","Missing.Percent")
  return(retour)
}