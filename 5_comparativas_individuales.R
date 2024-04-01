##############
# Código para comparar de forma individual cada proyecto en los distintos 
# resultados de expresión diferencial
##############

library(limma)
library(UpSetR)
library(venn)
library(grid)
library(GO.db)
library(topGO)

setwd("/media/scratch1/09MBIF/1_datos/3_anotacion/conteo/")

#comparativas <- c("NASH_F2vsF0-1", "NASH_F3vsF2", "NASH_F4vsF3"
                  #,"enfermedad")

proyectos <- c("PRJNA542148", "PRJNA558102", "PRJNA794736")

for (i in proyectos){
  print(i)
  NASH_F2vsF0_1 <- read.csv(file=paste0(i,"_DGE_NASH_F2vsF0-1.csv"))
  NASH_F3vsF2   <- read.csv(file=paste0(i,"_DGE_NASH_F3vsF2.csv"))
  NASH_F4vsF3   <- read.csv(file=paste0(i,"_DGE_NASH_F4vsF3.csv"))
  
  NASH_F2vsF0_1.UP   <- NASH_F2vsF0_1$ENTREZID[NASH_F2vsF0_1$diffexpressed=="UP"]
  print(length(NASH_F2vsF0_1.UP))
  NASH_F2vsF0_1.DOWN <- NASH_F2vsF0_1$ENTREZID[NASH_F2vsF0_1$diffexpressed=="DOWN"]
  print(length(NASH_F2vsF0_1.DOWN))
  NASH_F3vsF2.UP     <- NASH_F3vsF2$ENTREZID[NASH_F3vsF2$diffexpressed=="UP"]
  print(length(NASH_F3vsF2.UP))
  NASH_F3vsF2.DOWN   <- NASH_F3vsF2$ENTREZID[NASH_F3vsF2$diffexpressed=="DOWN"]
  print(length(NASH_F3vsF2.DOWN))
  NASH_F4vsF3.UP     <- NASH_F4vsF3$ENTREZID[NASH_F4vsF3$diffexpressed=="UP"]
  print(length(NASH_F4vsF3.UP))
  NASH_F4vsF3.DOWN   <- NASH_F4vsF3$ENTREZID[NASH_F4vsF3$diffexpressed=="DOWN"]
  print(length(NASH_F4vsF3.DOWN))
  
  #myList <- list(`NASH_F3vsF2` = NASH_F3vsF2$ENTREZID,
  #                  `NASH_F4vsF3` = NASH_F4vsF3$ENTREZID,
  #                  `NASH_F2vsF0_1` = NASH_F2vsF0_1$ENTREZID)
  
  #myList.UP <- list(`NASH_F3vsF2.UP` = NASH_F3vsF2.UP,
  #                     `NASH_F4vsF3.UP` = NASH_F4vsF3.UP,
  #                     `NASH_F2vsF0_1.UP` = NASH_F2vsF0_1.UP)
  
  #myList.DOWN <- list(`NASH_F3vsF2.DOWN` = NASH_F3vsF2.DOWN,
  #                      `NASH_F4vsF3.DOWN` = NASH_F4vsF3.DOWN,
  #                      `NASH_F2vsF0_1.DOWN` = NASH_F2vsF0_1.DOWN)
  
  myList.FULL <- list(`NASH_F0_1.UP - NASH_F2.DOWN` = NASH_F2vsF0_1.UP,
                      `NASH_F0_1.DOWN - NASH_F2.UP` = NASH_F2vsF0_1.DOWN,
                      `NASH_F2.UP - NASH_F3.DOWN` = NASH_F3vsF2.UP,
                      `NASH_F2.DOWN - NASH_F3.UP` = NASH_F3vsF2.DOWN,
                      `NASH_F3.UP - NASH_F4.DOWN` = NASH_F4vsF3.UP,
                      `NASH_F3.DOWN - NASH_F4.UP` = NASH_F4vsF3.DOWN)


  # UpSetPlot
  print("UpsetPlots")

  #jpeg(paste0("Upsetplot ", i, " FULL.jpeg"), width = 969, height = 448, quality = 100)
  upset(fromList(myList.FULL),nintersects =  NA, 
        sets = c("NASH_F0_1.UP - NASH_F2.DOWN","NASH_F0_1.DOWN - NASH_F2.UP","NASH_F2.UP - NASH_F3.DOWN","NASH_F2.DOWN - NASH_F3.UP",
                 "NASH_F3.UP - NASH_F4.DOWN","NASH_F3.DOWN - NASH_F4.UP"),show.numbers = "yes",
        matrix.color = "blue", 
        shade.color = c("orange", "red", "blue"),
        sets.bar.color = c("orange", "orange", "red","red", "blue", "blue"), 
        main.bar.color = "blue", text.scale = 2, keep.order = T)
  grid.text(paste0("Upset Plot comparación \n", i), x = 0.85, y = 0.87,
            gp = gpar(fontsize = 18))
  dev.off()
  
  #jpeg(paste0("Upsetplot ", i, " FULL.jpeg"), width = 969, height = 448, quality = 100)
  #upset(fromList(myList.UP),nintersects =  NA, 
  #      sets = c("NASH_F3vsF2.UP", "NASH_F4vsF3.UP", "NASH_F2vsF0_1.UP"),show.numbers = "yes",
  #      matrix.color = "blue", sets.bar.color = "blue", main.bar.color = "blue", text.scale = 2, keep.order = T)
  #grid.text(paste0("Upset Plot comparación \n", i), x = 0.85, y = 0.87,
  #          gp = gpar(fontsize = 18))
  #dev.off()
  
  #jpeg(paste0("Upsetplot ", i, " FULL.jpeg"), width = 969, height = 448, quality = 100)
  #upset(fromList(myList.DOWN),nintersects =  NA, 
  #     sets = c("NASH_F3vsF2.DOWN", "NASH_F4vsF3.DOWN", "NASH_F2vsF0_1.DOWN"),show.numbers = "yes",
  #     matrix.color = "blue", sets.bar.color = "blue", main.bar.color = "blue", text.scale = 2, keep.order = T)
  #grid.text(paste0("Upset Plot comparación \n", i), x = 0.85, y = 0.87,
  #          gp = gpar(fontsize = 18))
  #dev.off()
  
  
  commons<-c(NASH_F2vsF0_1.UP, NASH_F2vsF0_1.DOWN, NASH_F3vsF2.DOWN, 
             NASH_F3vsF2.UP, NASH_F4vsF3.DOWN, NASH_F4vsF3.UP)[ table(c(NASH_F2vsF0_1.UP, 
             NASH_F2vsF0_1.DOWN, NASH_F3vsF2.DOWN, NASH_F3vsF2.UP, NASH_F4vsF3.DOWN, NASH_F4vsF3.UP))==3]
  
  print(length(commons))
  print(NASH_F2vsF0_1$SYMBOL[NASH_F2vsF0_1$ENTREZID %in% commons])
  print(NASH_F3vsF2$SYMBOL[NASH_F3vsF2$ENTREZID %in% commons])
  print(NASH_F4vsF3$SYMBOL[NASH_F4vsF3$ENTREZID %in% commons])
  
  #GeneOntology
  print("GO terms")
  for (x in c("BP", "MF", "CC")){
    print(x)
    go     <- goana(commons, species= "Hs")
    top.go <- topGO(go, ontology = x, number = 500)
    top.go <- top.go[top.go$P.DE<= 0.05,]
    print(paste0("para ", i, " tiene ",dim(top.go)[1]))
    write.table(top.go, file = paste0(i, "_",x,"_GO.txt"), sep = "\t")
  }
}
