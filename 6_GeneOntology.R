################
# Código para estudio de Gene Ontology  en todos los proyectos
################

library(limma)
library(GO.db)
library(topGO)

setwd("/media/scratch1/09MBIF/1_datos/3_anotacion/conteo/")

comparativas <- c("enfermedad","NASH_F2vsF0-1", "NASH_F3vsF2", "NASH_F4vsF3")
proyectos <- c("PRJNA542148", "PRJNA558102", "PRJNA794736")
for (i in comparativas){
  print(i)
  PRJNA542148 <- read.csv(file=paste0("PRJNA542148_DGE_",i,".csv"))
  PRJNA558102 <- read.csv(file=paste0("PRJNA558102_DGE_",i,".csv"))
  PRJNA794736 <- read.csv(file=paste0("PRJNA794736_DGE_",i,".csv"))
  
  PRJNA542148.UP <- PRJNA542148$ENTREZID[PRJNA542148$diffexpressed=="UP"]
  print(length(PRJNA542148.UP))
  PRJNA542148.DOWN <- PRJNA542148$ENTREZID[PRJNA542148$diffexpressed=="DOWN"]
  print(length(PRJNA542148.DOWN))
  PRJNA558102.UP <- PRJNA558102$ENTREZID[PRJNA558102$diffexpressed=="UP"]
  print(length(PRJNA558102.UP))
  PRJNA558102.DOWN <- PRJNA558102$ENTREZID[PRJNA558102$diffexpressed=="DOWN"]
  print(length(PRJNA558102.DOWN))
  PRJNA794736.UP <- PRJNA794736$ENTREZID[PRJNA794736$diffexpressed=="UP"]
  print(length(PRJNA794736.UP))
  PRJNA794736.DOWN <- PRJNA794736$ENTREZID[PRJNA794736$diffexpressed=="DOWN"]
  print(length(PRJNA794736.DOWN))
  

  #GeneOntology
  print("GO terms")
  for (x in c("BP", "MF", "CC")){
    print(x)
    go <- goana(PRJNA558102$ENTREZID, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    print(paste0("PRJNA558102 para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA558102_", i, "_",x,"_GO.txt"), sep = "\t")
    go <- goana(PRJNA794736$ENTREZID,species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    dim(top.go)
    print(paste0("PRJNA794736 para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA794736_", i, "_",x,"_GO.txt"), sep = "\t")
    go <- goana(PRJNA542148$ENTREZID, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    dim(top.go)
    print(paste0("PRJNA542148 para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA542148_", i, "_",x,"_GO.txt"), sep = "\t")
    go <- goana(PRJNA558102.UP, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    dim(top.go)
    print(paste0("PRJNA558102 UP para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA558102_", i, "_",x,"_UP_GO.txt"), sep = "\t")
    go <- goana(PRJNA794736.UP,species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    dim(top.go)
    print(paste0("PRJNA794736 UP para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA794736_", i, "_",x,"_UP_GO.txt"), sep = "\t")
    go <- goana(PRJNA542148.UP, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    print(paste0("PRJNA542148 UP para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA542148_", i, "_",x,"_UP_GO.txt"), sep = "\t")
    go <- goana(PRJNA558102.DOWN, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    print(paste0("PRJNA558102 DOWN para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA558102_", i, "_",x,"_DOWN_GO.txt"), sep = "\t")
    go <- goana(PRJNA794736.DOWN,species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    print(paste0("PRJNA794736 DOWN para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0("PRJNA794736_", i, "_",x,"_DOWN_GO.txt"), sep = "\t")
    go <- goana(PRJNA542148.DOWN, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 3000)
    write.table(top.go, file = paste0("PRJNA542148_", i, "_",x,"_DOWN_GO.txt"), sep = "\t")
    print(paste0("PRJNA542148 DOWN para ", i, " tiene ",dim(top.go)[1]))
  }
}
