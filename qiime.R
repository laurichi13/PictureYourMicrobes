library(tidyverse)
library(qiime2R)
library(phyloseq)
library(microbiome)
library(viridis)

#colores####
library(randomcoloR)
paleta20 <- distinctColorPalette(k = 20, altCol = FALSE, runTsne = FALSE)
paleta20 <- randomColor(20,luminosity="bright")
paleta15 <- distinctColorPalette(k = 15, altCol = FALSE, runTsne = FALSE)
paleta10 <- distinctColorPalette(k = 10, altCol = FALSE, runTsne = FALSE)
paleta11 <- distinctColorPalette(k = 11, altCol = FALSE, runTsne = FALSE)
paleta5 <- distinctColorPalette(k = 5, altCol = FALSE, runTsne = FALSE)
paleta6 <- distinctColorPalette(k = 6, altCol = FALSE, runTsne = FALSE)
paleta3 <- randomColor(3, hue="random")

#Ver tabla del DADA2####
tabla_DADA2 <- read_tsv("denoising-stats/dada2_stats.tsv")
#Escoger el num_seqs minimas a quitar
#mean(tabla_DADA2$non_chimeric)-2*(sd(tabla_DADA2$non_chimeric))
#Unir tabla DADA mas metadata
stats_dada2 <- left_join(tabla_DADA2,metadata)
stats_dada2 <- stats_dada2[grep("-", stats_dada2$SampleID, invert=T),]
stats_dada2$Grupo <- factor(stats_dada2$Grupo, levels = c("Control", "Control-HFD", "PT15", "PT30", "RSV30"), labels = c("Control", "Control-HFD", "PT15", "PT30", "RSV30"))
gdada <- ggplot(stats_dada2, aes(SampleID, non_chimeric, fill=Grupo)) + geom_col()
gdada <- gdada + scale_fill_brewer(palette="Set1") + theme_bw() + xlab("Sample ID") + ylab("Nº of Sequences")
gdada <- gdada + ggtitle ("Nº of filtered & non-chimeric sequences") + theme(axis.text.x = element_text(color = "black", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5)) +
  theme(axis.title = element_text(size=rel(1.5)))
gdada + facet_wrap(~Grupo, scales = "free_x")    

#Again para pdf 
seqs1 <- metadata_seqs %>% filter(Tiempo=="Inicio") %>% filter(Grupo=="RSV30")
mean(seqs1$num_seqs)
sd(seqs1$num_seqs)
quantile(seqs1$num_seqs, prob = c(0.25, 0.5, 0.75), na.rm = TRUE)

gdada <- ggplot(metadata_seqs, aes(metadata_seqs$Grupo, metadata_seqs$num_seqs, fill=metadata_seqs$Grupo))
gdada <- gdada + geom_col() + scale_color_manual() + scale_fill_brewer(palette="Set1")
gdada <- gdada + theme_minimal() + xlab("Study Group") + ylab("Nº of Sequences")
gdada <- gdada +  ggtitle ("Nº of filtered & non-chimeric sequences") + 
  theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5)) +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(axis.title = element_text(size=rel(1.5)))
gdada + facet_wrap(~Grupo, scales = "free_x") + theme(legend.position = "none")

#Tabla ASVs####
ASV<-read_qza("Data/Qiime/table.qza")
#Taxonomia
taxonomyq <-read_qza("Data/Qiime/taxonomy.qza")
#taxonomyc <-read_qza("Data/Qiime/collapsed.qza")
#Transformar tabla
taxtableq <- taxonomyq$data %>% as.tibble() %>% separate(Taxon, sep=";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
taxtableq$Kingdom <- gsub("d__", "",taxtableq$Kingdom)
taxtableq$Phylum <- gsub("p__", "",taxtableq$Phylum)
taxtableq$Class <- gsub("c__", "",taxtableq$Class)
taxtableq$Order <- gsub("o__", "",taxtableq$Order)
taxtableq$Family <- gsub("f__", "",taxtableq$Family)
taxtableq$Genus <- gsub("g__", "",taxtableq$Genus)
taxtableq$Species <- gsub("s__", "",taxtableq$Species)

#Arbol
tree<-read_qza("Data/Qiime/rooted-tree.qza")
#Metadata
metadata<-read_csv("Data/metadata_informe.csv")

#Limpiar taxa: Quitar Arqueas y unassigned en reino, quitar unculuterd y unidentified en familia y genero
#taxtableq <- taxtableq[grep("Archaea|Unassigned", taxtableq$Kingdom, invert=T),]
#taxtableq <- taxtableq[grep("unidentified|uncultured", taxtableq$Family, invert=T),]
#taxtableq <- taxtableq[grep("unidentified|uncultured", taxtableq$Genus, invert=T),]

#Objeto phyloseq
OTUs <- otu_table(ASV$data, taxa_are_rows = T) 
tree <- phy_tree(tree$data) 
TAXq <- tax_table(as.data.frame(taxtableq) %>% select_("-Confidence") %>% column_to_rownames("Feature.ID") %>% as.matrix()) #moving the taxonomy to the way phyloseq wants it
sample_metadata <- sample_data(metadata %>% as.data.frame()) 
sample_names(sample_metadata) <- paste(metadata$SampleID)
physeq<- merge_phyloseq(OTUs, tree, TAXq,sample_metadata)

#Diversidades
#PCoA beta-diversidad Bray Curtis 
ordinated_taxa = ordinate(physeq, method="PCoA", distance="bray")
#Por grupo 
#Comp 1 y 2
plot_ordination(physeq, ordinated_taxa, color="BMI", title = "Bray-Curtis Principal Coordinates Analysis") + theme_bw()
#Comp 1 y 3
plot_ordination(physeq, ordinated_taxa, color="BMI", title = "Bray-Curtis Principal Coordinates Analysis", axes=c(1,3)) + theme_bw()
#Comp 2 y 3 
plot_ordination(physeq, ordinated_taxa, color="BMI", title = "Bray-Curtis Principal Coordinates Analysis", axes=c(2,3)) + theme_bw()

#Unifrac de todo
ordinated_taxa_unifrac = ordinate(physeq, method="MDS", distance="unifrac",weighted=TRUE)
plot_ordination(physeq, ordinated_taxa_unifrac, color="BMI", shape="Género", title = "Unifrac MDS Analysis") + theme_bw()
#Alphas
alphas = c("Shannon", "Simpson", "InvSimpson")
plot_richness(physeq, measures = alphas)
pairs(estimate_richness(physeq, measures = alphas))

#Ojo diversidades comparacion por grupo NO 
ptero <- prune_taxa(taxa_sums(physeq)>0, physeq)
plot_richness(ptero, x="Género", color = "BMI")
plot_richness(ptero, measures=c("Chao1", "Shannon", "Simpson"))
rich <- plot_richness(ptero, x="Género", color = "BMI", measures=c("Chao1", "Shannon", "Simpson")) + theme_bw() 
rich + geom_boxplot()
#Comparaciones alfa diversidad
alpha_ptero <- estimate_richness(ptero, measures=c("Chao1", "Shannon", "Simpson"))
alpha_ptero <- alpha_ptero %>% as.data.frame() %>% rownames_to_column("SampleID") %>% left_join(metadata)
#Diferencias en varianzas 
library(car)
leveneTest(y = alpha_ptero_final$Chao1, group = alpha_ptero_final$Grupo, center = "median")
#Todo junto
Chao1 <- kruskal.test(Chao1 ~ Grupo, data = alpha_ptero)
Shannon <- kruskal.test(Shannon ~ Grupo, data = alpha_ptero)
Simpson <- kruskal.test(Simpson ~ Grupo, data = alpha_ptero)

#Datos por muestra####

#Alfa diversidad 
rich <- alpha(physeq, index = "all")
pseq.rarified <- rarefy_even_depth(physeq)
rich_rare <- alpha(pseq.rarified, index = "all")

#Compositional 
pseq.compositional <- transform(physeq, "compositional")

#Unifrac
ordinated_taxa_unifrac = ordinate(physeq, method="MDS", distance="unifrac",weighted=TRUE)
plot_ordination(physeq, ordinated_taxa_unifrac, color="BMI", title = "Unifrac MDS Analysis") + theme_bw()

#Porcentaje Phylum 
Phyl <- tax_glom(pseq.compositional, taxrank = "Phylum") 
tb_phyl <- psmelt(Phyl) 
tbphylum <- tb_phyl %>% group_by(Phylum) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=5)
tbphylum <- tbphylum %>% inner_join(tb_phyl)
Plotf <- ggplot(tbphylum, aes(SampleID, Abundance ,fill=Phylum)) + geom_col(position="fill") + scale_fill_manual(values=paleta5)
Plotf <- Plotf + theme_minimal() + xlab("SampleID") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=11, face="italic"))
Plotf <- Plotf +  ggtitle ("") + theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plotf

#Cambiar nombres por nombre común
tbphylum$Phylum <- gsub("Actinobacteriota", "Actinobacteria",tbphylum$Phylum)
tbphylum$Phylum <- gsub("Bacteroidota", "Bacteroidetes",tbphylum$Phylum)
tbphylum$Phylum <- gsub("Verrucomicrobiota", "Verrucomicrobia",tbphylum$Phylum)
#Phylum por muestra####
phylsample <- tbphylum %>% filter(SampleID=="PYM000")
#Plot apilado####
Plotf <- ggplot(phylsample, aes(SampleID, Abundance ,fill=Phylum)) + geom_col(position="fill") + scale_fill_manual(values=paleta5)
Plotf <- Plotf + theme_minimal() + xlab("") + ylab("Abundancia Relativa") +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=11, face="italic"))


#Preprar datos y Label para elplot####
phylsample <- phylsample %>% select(Phylum,Abundance) %>% 
               mutate(Per=round(Abundance*100,1)) %>%
               mutate(Label=Per, Borrar="%") %>%
              unite(Label,Borrar,col="Label",sep="") %>% 
              mutate(Name=Phylum) %>% unite(Name,Label,col="Label",sep=" ") 

#Donut chart####              
#Reordenar datos 
phylsample <- arrange(phylsample,desc(Per)) %>% mutate(Phylum=factor(Phylum, levels=Phylum))  
donut <- ggplot(phylsample, aes(x = 2, y = Per, fill = Phylum)) +
  geom_col(position="fill") + scale_fill_viridis(discrete=TRUE, labels=phylsample$Label, name="") +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0)+ theme_void() +
  xlim(-1.5, 2.5) 
  donut <- donut + theme(legend.text = element_text(color = "black", size = 10, face = "italic"))

#Guardar plot####
png("Resultados/PYM000_phy2.png", width = 10, height = 5, units="in", res =300)
donut
dev.off()  

#Lollipopo chart####
 #Reordenar datos 
  phylsample <- arrange(phylsample,Per) %>% mutate(Phylum=factor(Phylum, levels=Phylum))  
  
#Grafico
loli <- ggplot(phylsample, aes(x=Phylum, y=Per)) +
        geom_segment(aes(x=Phylum, xend=Phylum, y=0, yend=Per), color="gray") +
        geom_point(aes(colour = factor(Phylum), size=4, alpha=0.8)) +
        scale_color_viridis(discrete=TRUE)   +
        theme_light() +
        coord_flip(expand = TRUE) +
        theme(panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank()) + 
        geom_text(aes(label=Label),hjust=-0.3) + 
        theme(axis.text.y = element_text(color = "black", size = 11, angle = 0, hjust = 0.5, vjust = 0, face = "italic")) + 
        theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
 
loli <- loli + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
loli <- loli + theme(legend.position = "none")  
loli


#Ratio Firmicutes/Bacteroidetes#### 
firm <- tb_phyl %>% select(SampleID, Abundance, Phylum) %>%
        filter(Phylum ==" Firmicutes")
colnames(firm) <- c("SampleID","Firmicutes","Phylum")
bact <- tb_phyl %>% select(SampleID, Abundance, Phylum) %>%
  filter(Phylum ==" Bacteroidota")
colnames(bact) <- c("SampleID","Bacteroidetes","Phylum")

ratio <- firm %>% full_join(bact,by="SampleID") %>% select(SampleID,Firmicutes,Bacteroidetes) %>% 
         mutate (ratio=Firmicutes/Bacteroidetes)

#Guardar resultados
#write_csv(ratio, "Resultados/ratio.csv")

#Porcentaje de géneros seleccionados (Top 10 - o específicos)####
#Convert to relative abundance
#physeq.rel  = transform_sample_counts(physeq.tree, function(x) x / sum(x))
#physeq.rel.fil = filter_taxa(physeq.rel, function(x) mean(x) > 1e-5, TRUE)
Gen <- tax_glom(pseq.compositional, taxrank = "Genus") 
tb <- psmelt(Gen) 
tb15genus <- tb %>% group_by(Genus) %>% summarise(M=median(Abundance)) %>% ungroup() %>% unique() %>% top_n(M,n=15)
tb15genus <- tb15genus %>% inner_join(tb)
Plot15g <- ggplot(tb15genus, aes(SampleID, Abundance ,fill=Genus)) + geom_col(position="fill") + scale_fill_manual(values=paleta15)
Plot15g <- Plot15g + theme_minimal() + xlab("SampleID") + ylab("Relative Abundance") +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 45, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=11, face="italic"))
Plot15g <- Plot15g +  ggtitle ("15 most abundant Genera") + theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plot15g

#Loop por muestra 
tb15sample <- tb15genus %>% filter(SampleID=="PYM000")
Plot15g <- ggplot(tb15sample, aes(SampleID, Abundance ,fill=Genus)) + geom_col(position="fill") + scale_fill_manual(values=paleta15)
Plot15g <- Plot15g + theme_minimal() + xlab("") + ylab("Abundancia Relativa") +
  theme(axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain")) +
  theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "plain")) + 
  theme(legend.title=element_text(color= "black", size=12), legend.text=element_text(size=11, face="italic"))
Plot15g



#Géneros de interes (Akkermasia, Lactobacillus y Bifido)####
int <- tb %>% select(SampleID,Genus,Abundance) %>% 
      filter(Genus == " Akkermansia" | Genus == " Bifidobacterium" | Genus == " Lactobacillus") 
#%>% spread(Genus,Abundance)
intsample <- int %>% filter(SampleID=="PYM000")
Plotint <- ggplot(intsample, aes(x=Genus,y=Abundance)) + geom_bar(stat='identity',aes(fill=Genus)) + 
 coord_flip() + scale_fill_manual(values=paleta3)
Plotint <- Plotint + theme_minimal() + xlab("") + ylab("Abundancia relativa") +
 theme(axis.text.x = element_text(color = "black", size = 14, angle = 0, hjust = .5, vjust = .5, face = "plain")) +
theme(axis.text.y = element_text(color = "black", size = 14, angle = 0, hjust = 1, vjust = 0, face = "italic")) + 
theme(legend.position="none")
Plotint <- Plotint +  ggtitle ("Géneros interesantes") + theme(plot.title = element_text(family="Arial", size=rel(2), vjust=2,face="plain", color="black", lineheight=1.5))  
Plotint

#Generos seleccionados####
int10 <- tb %>% select(SampleID,Genus,Abundance) %>% 
  filter(Genus == " Akkermansia" | Genus == " Bifidobacterium" | Genus == " Lactobacillus" |
         Genus == " Alistipes" | Genus == " Bacteroides" |
         Genus == " Prevotella" | Genus == " Faecalibacterium" | Genus == " Roseburia" | Genus == " Ruminococcus") %>% 
    mutate(Per=round(Abundance*100,1))

intsample <- int10 %>% filter(SampleID=="PYM000") 
#Añadir nueva fila
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
Plotint <- Plotint + scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))

#Guardar####
png("Resultados/PYM000.png", width = 8, height = 3, units="in", res=300)
Plotint
dev.off()

