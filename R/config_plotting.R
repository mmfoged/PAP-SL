
# Color palettes ----

## Gradients ----
#' Heatmap white to red 101-step gradient
color1 <- colorRampPalette(c("white", "#9F1931"))(101)

#' Heatmap colors blue-white-red 101-step gradient
color2 <- colorRampPalette(c("#003F7D", "white", "#9F1931"))(101)

#' color2 alternative heatmap colors
color2Alt <- colorRampPalette(c("#8C510A", "white", "#01665E"))(101)

## Headgroups ----
#' Headgroup colors
col_headgroups<-c(choline = "#BF6196",
                  ethanolamine = "#0F8F61",
                  inositol = "#4DA2DA",
                  serine = "#575756",
                  glycerol = "#DE8E0A")

## Organelles ----
#' Organelle colors
col_organelle<- c('Whole cells' = "#1D1D1B",
                  'LMF' = "#575756",
                  'Mitochondria' = "#288455",
                  'ER' = "#24425c",
                  'Cis-Golgi' = "#3d7199",
                  'Trans-Golgi' = "#6497c0",
                  'Plasma membrane' = "#a0c0d9",
                  'Lysosomes' = "#ae304f")

## Cells ----
#' Grey cell colors HKh-2 and HCT 116
col_cell<-c(HKH2 = "light grey", 
            HCT116 = "dark grey")

## Cells + organelles ----
#' cell_organelle colors
col_cell_organelle<- c("HKH2_Whole cells" = "#9D9D9C",
                       "HKH2_LMF" = "#C6C6C6",
                       "HKH2_Mitochondria"= "#F1C992",
                       "HKH2_ER" = "#A7C6AF",
                       "HKH2_Golgi" = "#9BAAD4",
                       "HKH2_Cis-Golgi" = "#9BAAD4",
                       "HKH2_Trans-Golgi"="#B6D1EF",
                       "HKH2_Plasma membrane"="#EAAB9E",
                       "HKH2_Lysosomes"="#DEB7CF",
                       "HCT116_Whole cells"="#1D1D1B",
                       "HCT116_LMF"="#575756",
                       "HCT116_Mitochondria" ="#DE8E0A",
                       "HCT116_ER"="#0F8F61",
                       "HCT116_Golgi"="#0E5DA3",
                       "HCT116_Cis-Golgi"="#0E5DA3",
                       "HCT116_Trans-Golgi"="#4DA2DA",
                       "HCT116_Plasma membrane"="#C94A17",
                       "HCT116_Lysosomes"="#BF6196")

## Categories  ----
#' Lipid category colors
col_cat<-c("GL" = "#B99AAD",
           "GPL"= "#5A5877",
           "GPLO-" = "#588C73",
           "LGPL" = "#F2AE72",
           "SL" = "#D96459",
           "ST" = "#8C4747",
           "LGPLO-" = "#E58965",
           "FA" = "grey")

## Classes  ----
#' Lipid class colors
col_class<- c("Chol" = "#8C4747",
              "CE" = "#A36B6B",
              "LCB" = "#D96459",
              "DihLCB" = "#DB6E64",
              "SM" = "#DE786F",
              "Cer" = "#E0837A",
              "CerP" = "#E38D85",
              "HexCer" = "#E59790",
              "diHexCer" = "#E8A29B",
              "triHexCer" = "#EAACA6",
              "GM3" = "#EDB6B1",
              "GM2" = "#EFC1BC",
              "GM1" = "#F2CBC7",
              "SHexCer" = "#F4D5D2",
              "DAG" = "#B99AAD",
              "TAG" = "#C0A5B6",
              "PA" = "#5A5877",
              "PC" = "#6A6884",
              "PI" = "#7B7992",
              "PE" = "#8B8A9F",
              "PS" = "#9C9AAD",
              "PG" = "#ACABBB",
              "CL" = "#BDBCC8",
              "LPA" = "#F2AE72",
              "LPC" = "#F3B781",
              "LPI" = "#F4C091",
              "LPE" = "#F6C9A1",
              "LPS" = "#F7D2B0",
              "LPG" = "#F9DBC0",
              "PAO-" = "#588C73",
              "PCO-" = "#689781",
              "PIO-" = "#79A38F",
              "PEO-" = "#8AAE9D",
              "PSO-" = "#9ABAAB",
              "PGO-" = "#ABC5B9",
              "LPAO-" = "#E58965",
              "LPCO-" = "#E79676",
              "LPIO-" = "#EAA387",
              "LPEO-" = "#EDB098",
              "LPSO-" = "#F0BDA9",
              "LPGO-" = "#F3CABA",
              "FA" = "#BEBEBE",
              "FAlc" = "#CBCBCB",
              "MAGO-" = "#C8B0BF",
              "DAGO-" = "#D0BBC8",
              "TAGO-" = "#D8C6D1",
              "ADAG" = "#DFD2DA",
              "PAF" = "#BCD1C7")

## Clusters  ----
#' Cluster colors for 8 clusters
col_clusters8 <- c('1'="#a56595",
                   '2'="#d06643",
                   '3'="#d783b1",
                   '4'="#b1552f",
                   '5'="#7cb691",
                   '6'="#c78c33",
                   '7'="#56936c",
                   '8'="#87223c")


# Plotting theme ----

#' ggplot theme function.
theme_print_function <- function(fontsize = 6) {
  theme_out <- ggplot2::theme_minimal() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = fontsize, face = "bold") ,
                   legend.text = ggplot2::element_text(size = fontsize, color = "black"),
                   legend.title = ggplot2::element_text(size = fontsize, color = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = fontsize, color = "black"),
                   axis.title.x = ggplot2::element_text(size = fontsize, color = "black"),
                   axis.text.y = ggplot2::element_text(size = fontsize, color = "black"),
                   axis.title.y = ggplot2::element_text(size = fontsize, color = "black"),
                   plot.title = ggplot2::element_text(size = fontsize, color = "black"),
                   plot.margin = ggplot2::margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "cm"),
                   axis.line.x.bottom =ggplot2::element_line(linewidth = .2),
                   axis.line.y.left = ggplot2::element_line(linewidth = .2),
                   plot.background = ggplot2::element_blank(),legend.key.size = ggplot2::unit(fontsize,units = "points"),
                   panel.grid.major.y = ggplot2::element_line(linewidth = 0.15, color = "grey", linetype = 2),
                   panel.background = ggplot2::element_blank()
    )
  
  return(theme_out)
  
}




