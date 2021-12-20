source("antiguaED.R")
#library("optparse")

#Funcion principal del algoritmo Evolucion diferencial
DE <- function(datos, tamanoP, fMutacion, CR){
  #datos es el dataset, numCluster cantidad de cluste que se van a crear, numGeneraciones cantidad de generaciones que realizara el algoritmo
  #tamanoP es el tamano de la poblacion, fMutacion es el factor de mutacion y CR es factor de cruce 
  tiempo <- proc.time()
  #iban arriba
  numCluster <- 4
  numGeneraciones <- 20
  #fMutacion <- round(runif(1, 0.3, 0.9), 1)
  #CR <- round(runif(1, 0.8, 0.9), 1)
  
  p <- NULL
  g <- 1
  bestSolucion <- NULL
  pos <- NULL
  
  #obtener cantidad de muestras y genes
  dimensionesDatos <- dim(datos)
  features <- dimensionesDatos[1]
  cantDatos <- dimensionesDatos[2]
  
  #Crear lista de la poblacion
  pNew <- list()
  pNew <- vector("list", length = tamanoP)
  
  #lista donde se guardan todas las generaciones 
  #pT <- list()
  #pT <- vector("list", length =  (numGeneraciones + 2))
  
  #matrizBox <- numeric(tamanoP * (numGeneraciones + 1)) 
  #dim(matrizBox) <- c(tamanoP, (numGeneraciones + 1))  
  
  stopDE <- FALSE 
  noImprovement <- 0
  limite <- 50
  #crear poblacion inicial
  p <- setup(datos, p, numCluster, features, tamanoP, cantDatos)
  
  #pT[[1]] <- p
  
  #Grafico de la poblacion incial 
  #matrizBox[,1] <- individuosRend(datos, p, numCluster, tamanoP)

  while (g <= numGeneraciones & stopDE == FALSE) {
    
    for(i in 1:tamanoP){
      #Se obtienen los individuos aleatorios 
      individuosRandom <- sample(1:tamanoP, 6, replace = FALSE)
      if(i != individuosRandom[1] && i != individuosRandom[2] && i != individuosRandom[3]){
        #se realiza el crossover y se obtiene el mejoor individuo
        iPrueba <- crossover(p[[i]], p[[individuosRandom[1]]], p[[individuosRandom[2]]], p[[individuosRandom[3]]], datos, fMutacion, CR, numCluster, features, cantDatos)
        iNew <- reemplazo(datos, p[[i]], iPrueba, features, numCluster)
        pNew[[i]] <- iNew
      }else{
        #se realiza el crossover y se obtiene el mejoor individuo
        iPrueba <- crossover(p[[i]], p[[individuosRandom[4]]], p[[individuosRandom[5]]], p[[individuosRandom[6]]], datos, fMutacion, CR, numCluster, features, cantDatos)
        iNew <- reemplazo(datos, p[[i]], iPrueba, features, numCluster)
        pNew[[i]] <- iNew
      }
    }
    
    #Si el fitness de p y pNew es igual durante  n generaciones stopDE cambia a true 
    if(g > (0.7 * numGeneraciones)){
      individuoActual <- bestIndividuo(datos, p, features, numCluster, tamanoP)
      individuoNew <- bestIndividuo(datos, pNew, features, numCluster, tamanoP)
      if(individuoActual[1] == individuoNew[1]){
        noImprovement <- noImprovement + 1
        if(noImprovement > limite){
          stopDE <- TRUE
        }
      }
    }
    
    #Se guarda la generacion
    #pT[[g + 1]] <- pNew
    
    #Grafica de los individuos de la  generacion nueva
    #matrizBox[,g + 1] <- individuosRend(datos, pNew, numCluster, tamanoP) 
    
    p <- pNew
    pNew <- list()
    pNew <- vector("list", length = tamanoP)
    g <- g + 1
    
    
  }
  
 # boxplot(x = matrizBox, y = matrizBox, xlab = "Numero de generacion ", ylab = "Valor de los individuos")
  #legend("topright", c("Factor mutación = 0.5", "Cruce = 0.9", "Tamano poblacion = 70" ))
  
  print(paste0((proc.time() - tiempo)))
  
  mejorSolucion <- bestIndividuo(datos, p, features, numCluster, tamanoP)
  #return(p)
  return(mejorSolucion[1])
}

#option_list = list(
 # make_option(c("--seed"), type = "integer"),
  #make_option(c("--fMutacion"), type = "integer"),
  #make_option(c("--CR"), type = "double"),
  #make_option(c("--tamanoP"), type = "double"),
  #make_option(c("--tries"), type = "integer"),
  #make_option(c("--time"), type = "integer"),
  #make_option(c("--quiet"), type = "logical"),
  #make_option(c("-i", "--input"), type = "character")
#);

#opt_parser = OptionParser(option_list=option_list);
#opt = parse_args(opt_parser);
#set.seed(opt$seed)

#datos <- read.csv(opt$input)

#resultado <- DE(datos ,opt$tamanoP, opt$fMutacion, opt$CR)

#cat(resultado)