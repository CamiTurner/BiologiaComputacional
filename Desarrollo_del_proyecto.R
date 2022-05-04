
# Características
# VIRUS: Sars-CoV-2
# NUCLEOTIDE COMPLETENESS: complete
# GEOGRAPHIC REGION: *especificar país*
# COLLECTION DATE: *registro más temprano*

# PANGO LINEAGE: B.1.1.519  ?????????????NO, saca las mismas :(????????????

rm(list=ls()) # Limpiar environment
cat("\014") # Limpiar consola

#install.packages("seqinr")
#install.packages("ggplot2")
library(seqinr)
library(ggplot2)

codonAamino = c(GAC="D",GAU="D", GAA="E",GAG="E", CGA="R",CGC="R",CGG="R",
                CGU="R",AGA="R",AGG="R", AAA="K",AAG="K", AAC="N",AAU="N", 
                CAC="H",CAU="H", CAA="Q",CAG="Q", UCA="S",UCC="S",UCG="S",
                UCU="S",AGC="S",AGU="S", ACA="T",ACC="T",ACG="T",ACU="T", 
                GCA="A",GCC="A",GCG="A",GCU="A", GGA="G",GGC="G",GGG="G",
                GGU="G", GUA="V",GUC="V",GUG="V",GUU="V", CCA="P",CCC="P",
                CCG="P",CCU="P", CUA="L",CUC="L",CUG="L",CUU="L",UUA="L",
                UUG="L", UUC="F",UUU="F", UAC="Y",UAU="Y", AUA="I",AUC="I",
                AUU="I", AUG="M", UGG="W", UGC="C",UGU="C", UAA="STOP",
                UAG="STOP",UGA="STOP")

# Encuentre la mutación dentro de un codón
diferencia = function(codonO, codonS){
  for(i in seq(1, 3, 1)){
    if(substr(codonO,i,i) != substr(codonS,i,i)){
      return (i);
    }
  }
}

dataFrame = data.frame(
  mutation = character(), 
  nucleotide = numeric(),
  codon = character(),
  protein = character(),
  index = numeric(),
  region = character()
)

original = read.fasta("Original_sequence_045512.txt") 
original = toupper( original[[3]] )
original[original == "T"] = "U" # Convertir de ADN a ARN

secuencias = list.files(path="./Secuencias", pattern=".txt", all.files=TRUE, full.names=TRUE)
cantMutaciones = 1

for(registro in seq_along(secuencias)){
  #print(secuencias[registro])
  sec = read.fasta( secuencias[registro] ) 
  sec = toupper( sec[[3]] )
  sec[sec == "T"] = "U"
  
  geo = substring(secuencias[registro], 14, 15) # AF:África, AS:Asia, EU:Europa, OC:Oceanía, AM:América
  cantCodones = 0 
  
  for(i in seq(1, min( length(original), length(sec) ), 3)){
    codonOr = paste(original[i], original[i+1], original[i+2], sep="")
    codonSec = paste(sec[i], sec[i+1], sec[i+2], sep="")
    
    if(codonOr != codonSec){
      posMutacion = diferencia(codonOr, codonSec)
      
      detectado = list( paste(substr(codonOr,posMutacion,posMutacion),"to",substr(codonSec,posMutacion,posMutacion)), (3*cantCodones+posMutacion), paste(codonOr,"to",codonSec), paste(codonAamino[codonOr],"to",codonAamino[codonSec]), cantCodones+1, geo )
      dataFrame[cantMutaciones, ] = detectado
      cantMutaciones = cantMutaciones+1
    }
    
    cantCodones = cantCodones+1
  }
}

#dataFrame

frecuencia = ggplot(dataFrame)
frecuencia = frecuencia + aes(x=mutation, fill=mutation, label=..count..)
frecuencia = frecuencia + ggtitle("Mutaciones de sustitución")
frecuencia = frecuencia + labs(x="Mutación", y="Frecuencia", fill="Mutaciones")
frecuencia = frecuencia + geom_bar(stat = "count")
frecuencia = frecuencia + geom_text(stat = "count", vjust = 1)
frecuencia = frecuencia + facet_grid(~"AS")
frecuencia
