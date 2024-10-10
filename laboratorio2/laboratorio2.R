#PRÁCTICA LABORATORIO2 BIOINFORMÁTICA
library(tidyverse)
data(starwars)

#FILTRAR Y SELECCIONAR DATOS
#Ejemplos usando select
# Seleccionar todas las columnas menos el nombre
starwars %>% select(-name)

#Seleccionar sólo las columnas que tienen subraya (_)
starwars %>% select(contains("_"))

#Seleccionar sólo las columnas que empiezan con "s"
starwars %>% select(starts_with("s"))

#Crear un data frame con los nombres y planeta de origen (homeworld)
homeworld <- starwars %>% select(name, homeworld)

#Filtrar datos 
#Filtrar por especies: sólo humanos
human <- starwars %>% filter(species == "Human")

#Filtrar por especies: sólo humanos del planeta Tatooine
starwars %>% filter(species == "Human", homeworld == "Tatooine")

#Crear un nuevo datframe con todas las especies menos los Droides
starwars_nodroids <- starwars %>% filter(species != "Droid")

#¿Cuántos registros cumplen las condiciones finales?
#77 personajes que no son especie Droide


#SELECCIONAR Y AGRUPAR DATOS
#Usamos group_by y tally
starwars %>% group_by(species) %>% tally()

#Añadiendo otra variable
starwars %>% group_by(species, gender) %>% tally()

#Si lo quieres guardar en el environment recuerda asignarle un nombre
table_gender <- starwars %>% group_by(species, gender) %>% tally()


#CALCULAR ALGUNOS ESTADÍSTICOS
#Calculamos la altura y la masa media por especie
starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))

#¿Cómo calcularías la desviación estándar (sd) de esos parámetros?
tabla_estadistica <- starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm=TRUE), mean_mass = mean(mass, na.rm = TRUE), sd_height = sd(height, na.rm=TRUE), sd_mass = sd(mass, na.rm = TRUE))


#CREAR GRÁFICOS Y MODIFICAR ALGUNOS ELEMENTOS
#Hacer un gráfico de la altura vs. la masa de los personajes
ggplot(starwars, aes(height, mass)) + geom_point()

#Puedes modificar el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")

#Modificando el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)

#Modificando el color y el fondo 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()


#QUITAR MASA GRANDE
#Descubrimos quien es el que mas pesa
masas <- starwars %>% select(name, mass)

#Ahora creamos un dataframe sin este personaje
starwars_sin_gordo <- starwars %>% filter(name != "Jabba Desilijic Tiure")

#Creamos el plot del nuevo dataframe
ggplot(starwars_sin_gordo, aes(height,mass)) + geom_point(colour = "red") + theme_light()


#EJERCICIO
#Cargamos en R studio el dataset toy.csv
toy <- read_csv("C:/Users/34644/Downloads/toy.csv")

#Inspecciona el dataset, haz un resumen de la media (mean) de las variables (Peso, Altura,IMC, IAS, CCintura). Agrupando por sexo.
media_variables_sexo <- toy %>% group_by(Sex) %>% summarise(mean_peso = mean(Weight_Kg, na.rm=TRUE), mean_altura = mean(Height_cm, na.rm=TRUE), mean_IMC = mean(IMC_clas, na.rm=TRUE), mean_IAS = mean(IAS, na.rm=TRUE), mean_Ccintura = mean(Ccintura, na.rm=TRUE))

#Haz una tabla sólo con los pacientes femeninos
pacientes_femenino <- toy %>% filter(Sex != "Men")

#¿Cuántos registros cumplen las condiciones?
#Hay 58 pacientes femeninos

#¿De estos cuantos tienen Sobrepeso (Overweight)?
pacientes_femenino %>% filter(IMC_clas == "Overweight") %>% tally()

#Haz un gráfico usando ggplot relacionando el IMC (Indice de masa corporal) con el peso (Weight_Kg) de todos los pacientes.
ggplot(toy, aes(IMC, Weight_Kg)) + geom_point()

#Repítelo filtrando sólo los pacientes categorizados como "Overweight" y "Obesity".
pacientes_sobrepeso_obesidad <- toy %>% filter(IMC_clas != "Normal")
ggplot(pacientes_sobrepeso_obesidad, aes(IMC, Weight_Kg)) + geom_point()