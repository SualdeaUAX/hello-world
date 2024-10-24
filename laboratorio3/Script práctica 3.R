library(ape)
library(phangorn)
library(phytools)
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta", 
                        format = "FASTA", type = "AA")
fraxatin

matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist

# A menor número más cercanos genéticamente están MÉTODO UPGMA
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)

# MÉTODO NJ, unión de vecinos
arbolNJ <- nj(matrizdist)
plot(arbolNJ)
  #modificar árbol característcias, colores, letras...
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolNJ, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)

# Personalizar árbol con pythools, otra librería par aconseguir el árbol.
plotTree(arbolNJ)
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)

#Ordenar de más tiempo a menos tiempo evolucionando.
plotTree(ladderize(arbolNJ))
  #Guardar el árbol
write.tree(arbolNJ, file = "arbolNJ.nex")
read.tree(file = "arbolNJ.nex")

#Ahora mismo los árboles que tenemos no están enraizados.
#Para enraizar ponemos nombre en el outgroup de la secuencia que nos interese.
arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)
  #Esto lo que nos permite es ordenar de una forma el árbol para que encontremos las similitudes en la secuencia

  #También con el árbol UPGMA, que nos da el árbol pero sin teimpos evolutivos
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)

#VISULAIZAR LOS ÁRBOLES A LA VEZ, en base a los árboles que tienen raíz
layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)

#ARBOLES DE PARSIMONIA
  #Buscar disminuir el número de pasos para hacer los árboles por lo tanto el árbol que se haya hecho con menos pasos es el óptimo y con ese nos quedamos
#No tenemos en cuenta las carcateres constantes o no informativos entre secuencias.
#Manera de estimar los árboles de máx parsimonia, número de pasos.
parsimony(arbolUPGMAraiz, fraxatin)
    ## [1] 313, Número de pasos por el cual se ha realizado este árbol

#Con el siguiente comando obtenemos el árbol con menos número de pasos que nos interesa.
  #Esté con raíz o no el número de pasos tiene que ser el mismo.
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)

  #Para visualizarlos 
plot(mejorNJ)
plot(mejorUPGMA)

#Cuando hay muchos taxones no es eficiente utilizar este, teniendo que realizar búsquedas heuríticas donde no se evaúan todos los posibles árboles


#Otro método para encontrar el mejor árbol de parsimonia, con el algoritmo.
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
fraxatin_parsimonia
  #= 4, posible árboles filogenéticos con mismos pasos e igualdad de longitud

#Para compararlos en un plot usamos el siguiente comando.
#Para poder compararlos hay que enraizarlos.
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)


#Vamos a realizar un árbol consenso el cual vamos a usar para obtener el arbol de parsimonia ideal.
#Si tenemos un consenso estricto nos incluirá los grupos monofiléticos presentes en todos los arboles supuestos.
#Podemos elegir además el grado de restricción.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)
  #Este árbol tiene una restricción del 100%

estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)
  #Este árbol tiene una restricción del 30%, por lo tanto esto nos proporciona menos politomías 
#Ya que no es tan estricto como el anterior.


#BOOTSTRAP
#Creación de réplicas de los alineaminetos eliminando parejas de la secuancia aleatoriamente
#Con cada réplica crearemos un árbol y esto nos dará la información de cuantas veces aparece un nodo específico en cada árbol
#Mientras más veces salga el nodo en cada réplpica diferente ese nodo tiene mayor probabilidad de ser real.

arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
plot(arbolesbootstrap, cex = .6)
#Cambiando el número de bs te da un nuevo número de génesis de las réplicas.
#Ahora generamos un consenso al 60%.
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)


#MODELOS PROBABILÍSTICOS
#ÁRBOLES DE MÁX VEROSIMILITUD

#Se genera un árbol enraizado de inicio de cualquier tipo, puede ser al azar.
#Se calcula la verosimilitud de casa sitio (de aminoácidos) usando un árbol.
#Se calcula la verosimilitud total del árbol por nodo. La verosimilitud es la suma de la reconstrucción de este sitio dado modelo de sustitución o de evolución.
# n=11, porque es el número de secuencias, árbol de 11 ramas.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)

#Para visualizarlo mejor lo enraizamos con secuencia ornitorrinco.
#Lo escalerizamos hacia la derecha y agragamos escala, longitud rama significativa, indica cantidad de cambios de aa.

arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()

#Con este árbol empezamos la búsqueda de un árbol mejor por máxima verosimilitud. Calculando verosimilitud dadas las secuencias.
#Con pml podemos computar la maxima verosimilitud.
ajustado <- pml(arbolazarR, fraxatin)
ajustado

#model: Mk 
#loglikelihood: -4318.749 
#unconstrained loglikelihood: -1479.871 
#Rate matrix: 

#Unconstrained loglikelihood reporta la verosimilitud de un modelo de substitución general que puede que no se ajuste bien a los datos.

#Ahora debemos encontrar un árbol que optimice la verosimilitud usando modelo de sustitución.
  #optim.pml computa árbol filogenético dado alineamiento multiple y modelo de evolución de aa.
  ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

#Para ver el árbol usamos $tree. 
ajustadoconDay$tree

#También lo enraizamos.
ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()

#Se puede usar otros modelos de sustitución.Como blosum o JTT
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")
#Comparamos los modelos según criterio de Akaike AIC:
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)

#df      AIC
#ajustadoconDay 19 5188.886
#ajustadoconBlo 19 5197.264
#ajustadoconJTT 19 5061.388

#Primera columna son los grados de libertad. El mejor modelo sería el JTT ya que tiene el AIC más bajo.

mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")

#model: JTT 
#loglikelihood: -2511.694 
#unconstrained loglikelihood: -1479.871 
#Rate matrix: JTT 

mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()

#ÁRBOL FINAL










