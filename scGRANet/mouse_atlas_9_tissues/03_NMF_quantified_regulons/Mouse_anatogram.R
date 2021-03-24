devtools::install_github("jespermaag/gganatogram")
library(gganatogram)
library(dplyr)
library(viridis)
library(gridExtra)





mmMale <- gganatogram(data=mmMale_key, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour") + theme_void()
mmMale
mmMale <- gganatogram(data=mmMale_key, fillOutline='#440154FF', organism='mouse', sex='male', fill="value") + theme_void() +  scale_fill_viridis()

p1 <- mmMale_key %>% 
  dplyr::filter(organ %in% c("bone_marrow", "brain", "kidney", "liver", "lung", "small_intestine", "spleen", "testis", "thymus")) %>% 
  gganatogram(outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +
  theme_void()  +
  coord_fixed()
p1

p2 <- mmMale_key %>% 
  dplyr::filter(organ %in% c("bone_marrow", "brain", "kidney", "liver", "lung", "small_intestine", "spleen", "testis", "thymus")) %>% 
  gganatogram(outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +
  facet_wrap(~type, ncol=4) +
  theme_void()  +
  coord_fixed()
p2

l <- lapply(c("bone_marrow", "brain", "kidney", "liver", "lung", "small_intestine", "spleen", "testis", "thymus"), function(tissue){
  mmMale_key %>% 
    dplyr::filter(organ %in% tissue) %>% 
    gganatogram(outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +
    ggtitle(tissue) +
    theme_void()  +
    coord_fixed()
})

p3 <- patchwork::wrap_plots(l)
p3

p4 <- mmMale_key %>% 
  mutate(colour = if_else(organ %in% c("bone_marrow", "brain", "kidney", "liver", "lung", "small_intestine", "spleen", "testis", "thymus"),
                         "red", "steelblue")) %>% 
  gganatogram(outline = T, fillOutline='#a6bddb', organism='mouse', sex='male', fill="colour")  +
  theme_void()  +
  coord_fixed()
p4


pdf(file = "/home/bq_aquintero/projects/mouse_atlas_TF_activity/adult_9tissues/mouse_anatogram.pdf", width=10, height=15)
for (p in list(p1, p2, p3, p4)) {
  plot(p)
}
dev.off()
