#Loop para el Donut plot 
sample <- tbphylum$SampleID
donut_list <- list ()
for(i in sample){
  phylsample <- tbphylum %>% filter(SampleID==i)
  phylsample <- phylsample %>% select(Phylum,Abundance) %>% 
    mutate(Per=round(Abundance*100,1)) %>%
    mutate(Label=Per, Borrar="%") %>%
    unite(Label,Borrar,col="Label",sep="") %>% 
    mutate(Name=Phylum) %>% unite(Name,Label,col="Label",sep=" ") 
  phylsample <- arrange(phylsample,desc(Per)) %>% mutate(Phylum=factor(Phylum, levels=Phylum))  
  donut <- ggplot(phylsample, aes(x = 2, y = Per, fill = Phylum)) +
    geom_col(position="fill") + scale_fill_viridis(discrete=TRUE, labels=phylsample$Label, name="") +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+ theme_void() +
    xlim(-1.5, 2.5)
  donut <- donut + theme(legend.text = element_text(color = "black", size = 10, face = "italic"))
  donut_list[[i]] = donut
   } 

# Save plots to png. Makes a separate file for each plot.
for (i in sample) {
  file_name = paste("donut", i, ".png", sep="")
  png(file_name, width = 10, height = 5, units="in",res=300)
  print(donut_list[[i]])
  dev.off()
}

#Loop para las barras ####
sample <- tbphylum$SampleID
bar_list <- list ()
for(i in sample){
    intsample <- int10 %>% filter(SampleID==i) 
    B <- data.frame(Genus="Otros géneros", Per=summarise(intsample,Per=100-sum(Per)))
    intsample <- intsample %>% rows_upsert(B)
    #Cambiar nombres añadiendo spp. 
    intsample$Genus <- gsub("Prevotella", "Prevotella spp.", intsample$Genus)
    intsample$Genus <- gsub("Bacteroides", "Bacteroides spp.", intsample$Genus)
    intsample$Genus <- gsub("Faecalibacterium", "Faecalibacterium spp.", intsample$Genus)
    intsample$Genus <- gsub("Alistipes", "Alistipes spp.", intsample$Genus)
    intsample$Genus <- gsub("Roseburia", "Roseburia spp.", intsample$Genus)
    intsample$Genus <- gsub("Ruminococcus", "Ruminococcus spp.", intsample$Genus)
    intsample$Genus <- gsub("Bifidobacterium", "Bifidobacterium spp.", intsample$Genus)
    intsample$Genus <- gsub("Akkermansia", "Akkermansia spp.", intsample$Genus)
    intsample$Genus <- gsub("Lactobacillus", "Lactobacillus spp.", intsample$Genus)
        #Crear "label" con el símbolo de %
    intsample <- intsample %>% mutate(Label=Per, Borrar="%") %>% 
      unite(Label,Borrar,col="Label",sep="") 
        #Reordenar datos para gráfica de más abundante a menos
    intsample <- arrange(intsample,Per) %>% mutate(Genus=factor(Genus, levels=Genus))
        #Gráfica seleccionada####
    Plotint <- ggplot(intsample, aes(x=Genus,y=Per)) + geom_bar(stat='identity',aes(fill=Genus)) + 
      coord_flip() + scale_fill_viridis(option="magma", discrete=TRUE)
    Plotint <- Plotint + theme_void() + 
      theme(axis.text.y = element_text(color = "black", size = 10, angle = 0, hjust = 0.5, vjust = 0, face = "italic")) + 
      theme(legend.position="none")
    Plotint <- Plotint + geom_text(aes(label=Label),size=3, color="black", vjust=0.5, hjust=-0.15)
    bar <- Plotint + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
    bar_list[[i]] = bar
} 

# Save plots to png. Makes a separate file for each plot.
for (i in sample) {
  file_name = paste("bar", i, ".png", sep="")
  png(file_name, width = 8 , height = 3 ,units="in", res=300)
  print(bar_list[[i]])
  dev.off()
}


