Plot_pca <- function(pca_result, rango_axis=1:2, vector_grupos=NULL, 
                     italic_labels=TRUE, scale_factor=1, use_ellipse=TRUE,
                     ellipse_alpha=0.2, ellipse_lwd=0.75, name_groups="Grupos", 
                     show_group_shape=TRUE, point_size = 3, point_alpha=0.4) {
  
  # Librerías
  require(ggrepel)
  require(tidyverse)
  
  # Calculando los porcentajes de varianza explicada
  explained_var <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  xlab <- paste0("PC", rango_axis[1], " (", round(explained_var[rango_axis[1]], 2), "%)")
  ylab <- paste0("PC", rango_axis[2], " (", round(explained_var[rango_axis[2]], 2), "%)")
  
  # Bases de datos para graficar
  datos_samples <- pca_result$x %>% as.data.frame() %>% select(all_of(rango_axis))
  datos_variables <- pca_result$rotation %>% as.data.frame() %>% select(all_of(rango_axis))
  
  # Multiplicar las cargas por scale_factor
  datos_variables <- datos_variables * scale_factor
  
  # Condición para las etiquetas en cursiva
  if(italic_labels) {
    label_face <- "italic"
  } else {
    label_face <- "plain"
  }
  
  # Gráfico base
  plot <- ggplot() + theme_bw() + labs(x=xlab, y=ylab)
  
  # Condicional en vector_grupos
  if(is.null(vector_grupos)) {
    plot <- plot + 
      geom_point(data=datos_samples, aes(x=PC1, y=PC2), size = point_size, alpha = point_alpha) 
  } else {
    datos_samples$grupo <- vector_grupos
    
    if(show_group_shape) {
      if(use_ellipse) {
        plot <- plot + 
          stat_ellipse(data=datos_samples, aes(x=PC1, y=PC2, color=grupo), geom="path", lwd=ellipse_lwd)
      } else {
        # Cálculo de convex hull
        hull <- datos_samples %>%
          group_by(grupo) %>%
          slice(chull(PC1, PC2))
        
        plot <- plot + 
          geom_polygon(data=hull, aes(x=PC1, y=PC2, fill=grupo, color=grupo), alpha=ellipse_alpha)
      }
    }
    
    plot <- plot +
      geom_point(data=datos_samples, aes(x=PC1, y=PC2, color=grupo), size = point_size, alpha = point_alpha) 
  }
  
  # Dibujar flechas para las variables y las etiquetas
  plot <- plot +
    geom_segment(data=datos_variables, aes(x=0, y=0, xend=PC1, yend=PC2), 
                 arrow=arrow(length=unit(0.2, "cm")), color="gray30") +
    ggrepel::geom_text_repel(data=datos_variables, aes(x=PC1, y=PC2, label=rownames(datos_variables)), fontface = label_face) +
    guides(color=guide_legend(override.aes=list(label=rep("", length(unique(vector_grupos))))))
  
  plot <- plot + labs(color=name_groups, fill=name_groups)
  
  return(plot)
}
