library(VennDiagram)
library(ggvenn)
library(ggsci)
library(RColorBrewer)
library(ggrepel)
library(ggforce) # for 'geom_arc_bar'

########### draw venna grame

source("scratch/18.NASH/Rscripts/NASH_preprocessing_data.R")
metabolites_CM <- CM.meta.dict$Compound
metabolites_liver <-liver.meta.dict$Compound 
metabolites_plasma <- plasma.meta.dict$Compound

liver.meta.count$`412` %>% sum()


liver.meta.dict$Batch %>% as.data.frame() %>% mutate_if(is.character, as.factor) %>% table()

list_batch <- c("HILIC_Neg","HILIC_Pos","RP_Pos","RP_Neg")

liver.meta.dict %>% filter(Batch == "HILIC_Neg" ) %>% dplyr::select(super.class..HMDB.) -> liver.meta.dict.1
liver.meta.dict %>% filter(Batch == "HILIC_Pos" ) %>% dplyr::select(super.class..HMDB.) -> liver.meta.dict.1
liver.meta.dict %>% filter(Batch == "RP_Neg" ) %>% dplyr::select(super.class..HMDB.) -> liver.meta.dict.1
liver.meta.dict %>% filter(Batch == "RP_Pos" ) %>% dplyr::select(super.class..HMDB.) -> liver.meta.dict.1

list_batch <- c("HILIC_Neg", "HILIC_Pos")
list_batch <- c("RP_Neg", "RP_Pos")


get_total <- function(meta.dict, batch) {
  liver.meta.dict %>%
    filter(Batch == batch) %>% dplyr::select(super.class..HMDB.) %>%
    mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- ")) %>%
    unnest(super.class..HMDB.) %>%
    group_by(super.class..HMDB.) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    head(n = 8) -> Total.1
  
  liver.meta.dict %>%
    filter(Batch == batch) %>% dplyr::select(super.class..HMDB.) %>% 
    mutate(super.class..HMDB. = strsplit(super.class..HMDB., " --- ")) %>%
    unnest(super.class..HMDB.) %>%
    group_by(super.class..HMDB.) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    mutate(super.class..HMDB. = ifelse(row_number() > 8, "SOther", as.character(super.class..HMDB.))) %>%
    filter(super.class..HMDB. == "SOther") %>%  group_by(super.class..HMDB.) %>% 
    summarise(n = sum(n)) -> Total.2
  
  Total <- rbind(Total.1, Total.2)
  
  return(Total)
}

Total

plots <- lapply(list_batch, function(batch) {
  Total_table <- get_total(liver.meta.dict, "HILIC_Pos")
  sumMetabolites <- Total_table %>% dplyr::select(n) %>% sum()
  Total_table
  ggplot(Total_table,
         #mutate(super.class..HMDB. = factor(super.class..HMDB., levels = unique(super.class..HMDB.))),
         aes(x = "", y = n, fill = super.class..HMDB.)) +
    geom_bar(width = 0.2, stat = "identity", color = "white") +
    geom_text(aes(label = n, y = n),
              position = position_stack(vjust = 0.5),
              size = 3.4) +
    coord_polar("y", start = 0) +
    theme_void() +
    scale_fill_npg() +
    theme(legend.position = "right") +
    ggtitle(paste(batch, "\n", sumMetabolites, "metabolites"))
})

# Display the four plots
multiplot <- function(...) {
  plot_list <- list(...)
  gridExtra::grid.arrange(grobs = plot_list, ncol = 2)
}

pdf("HILIC_metabolite_pie.pdf", width = 10, height = 3)
multiplot(plots[[1]], plots[[2]])
dev.off()
pdf("RP_metabolite_pie.pdf", width = 10, height = 3)
multiplot(plots[[1]], plots[[2]])
dev.off()



########################################################
24+1+38+36+26+6+75+18+62
14+36+108+39+17+130+15+58+59

# Calculate the percentage of each super class
percentage_table <- liver.meta.dict_separated %>%
  group_by(super.class..HMDB.) %>%
  summarise(Percentage = n() / nrow(liver.meta.dict_separated) * 100)

# Print the result
percentage_table %>% as.data.frame() %>% arrange(desc(Percentage)) %>% head(n=10)



# Order the data frame by Percentage in descending order and select the top 10
top_10 <- percentage_table %>%
  arrange(desc(Percentage)) %>%
  head(n = 10) %>% arrange(super.class..HMDB.)


color <- brewer.pal(length(10), "Set1") 
# Create a pie chart with percentage labels
pie(top_10$Percentage, labels = paste(top_10$super.class..HMDB., scales::percent(top_10$Percentage / 100)), main = "Top 10 Super Classes")
pie(top_10$Percentage, col=color,radius = 1, srt=15,
    labels = paste0(round(top_10$Percentage,2),"%"), main = "Top 10 Super Classes")


# Order the data frame by Percentage in descending order and select the top 10
top_10 <- top_10 %>%
  arrange(desc(Percentage)) %>%
  head(n = 10) %>% arrange(super.class..HMDB.)




################################################################################
################################################################################
################################################################################


# Create a color vector for the pie chart
color <- rainbow(length(top_10$super.class..HMDB.))

# Create a pie chart with labels using ggplot2 and ggrepel
top_10 %>%  mutate(csum = rev(cumsum(rev(Percentage))), 
                   pos = Percentage/2 + lead(csum, 1),
                   pos = if_else(is.na(pos), Percentage/2, pos)) %>%
  ggplot( aes(x = "", y = Percentage, fill = super.class..HMDB.)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") +
  geom_label_repel(aes(y=pos,
                       label = ifelse(Percentage >= 4, paste0(round(Percentage, 2), "%"), "")), 
                   label.size=NA,
                   label.r=0.7, 
                   size = 5, box.padding = -0.5) +  # Adjust label size and box.padding
  labs(title = "Top 10 Super Classes")



liver.meta.dict$Compound %>% sample(12)
liver.meta.dict$Compound %>% unique() %>% length()
liver.meta.dict$Compound %>% length()
# Create a list of sets
venn_list <- list(
  CM = metabolites_CM,
  Liver = metabolites_liver,
  Plasma = metabolites_plasma
)


# Create Venn diagram
venn.plot <- venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("CM", "Liver", "Plasma"),
  filename = NULL,
  output = TRUE
)
metabolites_CM %>% length()


# Display the diagram
grid.draw(venn.plot)
help(ggvenn)
pdf("vennaDiagram.pdf")
ggvenn(venn_list, c("Liver", "CM"), fill_alpha=0.5, set_name_size=8, fill_color=c("yellow", "pink"), 
       stroke_color = "white", digits=0, 
       text_size = 6)
dev.off()

liver.meta.dict %>% mutate_if(is.character, as.factor) %>% dplyr::select(`super.class..HMDB.`) %>% summary()
liver.meta.dict %>% mutate_if(is.character, as.factor) %>% dplyr::select(contains("class")) %>% dplyr::select(c(1,4)) %>% summary()
