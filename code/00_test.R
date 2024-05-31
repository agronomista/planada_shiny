library(shiny)
library(ggplot2)
library(sf)
library(dplyr)
library(readr)
library(gstat)
library(raster)
library(sp)
library(akima)
library(plotly)
library(shinyjs)

# Vector de nombres de malezas
weed_names <- c(
  "Bidens pilosa L.",
  "Caperonia palustris (L.) A.St.-Hil.",
  "Commelina diffusa Burm.f.",
  "Cynodon nlemfuensis Vanderyst",
  "Cyperus esculentus L.",
  "Digitaria sp.",
  "Echinochloa colonum (L.) Link",
  "Eleusine indica (L.) Gaertn.",
  "Emilia fosbergii Nicolson",
  "Euphorbia heterophylla L.",
  "Euphorbia hirta L.",
  "Ixophorus unisetus (J.Presl) Schltdl.",
  "Kyllinga sp.",
  "Melampodium divaricatum (Rich.) DC.",
  "Melampodium perfoliatum (Cav.) Kunth",
  "Portulaca oleraceae L.",
  "Rottboellia cochinchinensis (Lour.) Clayton"
)

ui <- fluidPage(
  useShinyjs(),  # Inicializar shinyjs
  titlePanel("Interpolaciones de Malezas"),
  sidebarLayout(
    sidebarPanel(
      selectInput("interpolation", "Seleccione el tipo de interpolación:",
                  choices = c("IDW" = "idw", "Spline" = "spline", "Moda" = "mode")),
      actionButton("generate", "Generar Interpolación"),
      br(),
      textOutput("description"),
      br(),
      fluidRow(
        column(6, downloadButton("downloadRaster", "Descargar Raster")),
        column(6, actionButton("showStats", "Estadísticas"))
      ),
      br(),
      plotlyOutput("histogram", height = "400px")
    ),
    mainPanel(
      plotOutput("plot", height = "800px")  # Aumentar el tamaño del plot
    )
  )
)

server <- function(input, output, session) {
  
  # Variable reactiva para almacenar el raster generado
  rasterData <- reactiveVal(NULL)
  
  observeEvent(input$generate, {
    
    # Cargar los datos
    data <- read_delim("C:/Users/usuario/Documents/PAMELA_puntos_28_5_2024/levantamiento_malezas/data/Merged_Data_with_Coordinates.csv", delim = ";", show_col_types = FALSE)
    
    # Asegúrate de que las columnas de latitud y longitud tengan los nombres correctos
    points_sf <- st_as_sf(data, coords = c("x_proj", "y_proj"), crs = 4326)
    points_sf <- st_transform(points_sf, crs = 32632) # CRS UTM zona 32N, ajústalo según tu región
    points_sp <- as(points_sf, "Spatial")
    
    # Crear un grid para la interpolación con mayor resolución
    grid <- as.data.frame(spsample(points_sp, "regular", n = 500000))
    names(grid) <- c("x", "y")
    coordinates(grid) <- c("x", "y")
    gridded(grid) <- TRUE
    fullgrid(grid) <- TRUE
    proj4string(grid) <- proj4string(points_sp)
    
    if (input$interpolation == "idw") {
      # Realizar la interpolación IDW
      idw_result <- idw(ID ~ 1, points_sp, newdata = grid, idp = 2.0)
      r <- raster(idw_result)
      crs(r) <- CRS(st_crs(points_sf)$proj4string)
      interp_df <- as.data.frame(r, xy = TRUE)
      names(interp_df)[3] <- "ID"
      rasterData(r)
      
    } else if (input$interpolation == "spline") {
      # Convertir a data frame para la interpolación
      df <- as.data.frame(st_coordinates(points_sf))
      df$ID <- points_sf$ID
      
      # Realizar la interpolación spline
      interpolated <- interp(x = df$X, y = df$Y, z = df$ID, duplicate = "mean", linear = TRUE, extrap = FALSE)
      r <- rasterFromXYZ(data.frame(x = rep(interpolated$x, length(interpolated$y)),
                                    y = rep(interpolated$y, each = length(interpolated$x)),
                                    z = as.vector(interpolated$z)))
      crs(r) <- CRS(st_crs(points_sf)$proj4string)
      interp_df <- as.data.frame(r, xy = TRUE)
      names(interp_df)[3] <- "ID"
      rasterData(r)
      
    } else if (input$interpolation == "mode") {
      # Crear un grid de mayor resolución
      x_range <- range(points_sp@coords[,1])
      y_range <- range(points_sp@coords[,2])
      res <- 100  # Ajustar resolución del grid
      grid <- expand.grid(x = seq(from = x_range[1], to = x_range[2], length.out = res),
                          y = seq(from = y_range[1], to = y_range[2], length.out = res))
      
      # Función para encontrar la moda
      mode_function <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      
      # Aplicar la interpolación por moda
      grid$ID <- apply(grid, 1, function(row) {
        dists <- sqrt((points_sp@coords[,1] - row['x'])^2 + (points_sp@coords[,2] - row['y'])^2)
        nearest_indices <- order(dists)[1:5]  # Usar los 5 puntos más cercanos
        mode_function(points_sp@data$ID[nearest_indices])
      })
      
      # Asignar nombres en lugar de IDs
      grid$Weed <- weed_names[grid$ID + 1]  # +1 porque los IDs comienzan desde 0
      r <- rasterFromXYZ(grid)
      crs(r) <- CRS(st_crs(points_sf)$proj4string)
      interp_df <- as.data.frame(r, xy = TRUE)
      names(interp_df)[3] <- "Weed"
      rasterData(r)
    }
    
    output$plot <- renderPlot({
      ggplot() +
        geom_raster(data = interp_df, aes(x = x, y = y, fill = Weed)) +
        geom_sf(data = points_sf, color = "red") +
        scale_fill_manual(values = viridis::viridis(length(weed_names)), labels = weed_names, name = "Maleza") +
        theme_minimal() +
        labs(title = ifelse(input$interpolation == "idw", "Interpolación IDW de IDs de Malezas",
                            ifelse(input$interpolation == "spline", "Interpolación Spline de IDs de Malezas",
                                   "Interpolación por Moda de IDs de Malezas")),
             x = "Longitud", 
             y = "Latitud")
    })
    
    output$description <- renderText({
      if (input$interpolation == "idw") {
        "IDW (Inverse Distance Weighting): Este método de interpolación asigna valores a los puntos desconocidos basado en una combinación ponderada de valores de puntos conocidos cercanos. Los puntos más cercanos tienen un peso mayor que los más alejados. La interpretación de los resultados es que las áreas con más datos cercanos tendrán predicciones más precisas."
      } else if (input$interpolation == "spline") {
        "Spline: Este método de interpolación utiliza funciones matemáticas suaves para generar una superficie que pasa exactamente por los puntos de datos conocidos. Es útil para obtener una superficie suavizada y continua. La interpretación de los resultados es que las variaciones en los datos se suavizan, lo que puede ser útil para identificar tendencias generales."
      } else if (input$interpolation == "mode") {
        "Moda: Este método de interpolación asigna a cada punto del grid la categoría (ID de maleza) que más frecuentemente aparece entre los puntos de datos más cercanos. Es útil para clasificaciones categóricas donde se desea identificar la categoría más común en una región específica."
      } else {
        ""
      }
    })
  })
  
  observeEvent(input$showStats, {
    # Cargar los datos
    data <- read_delim("C:/Users/usuario/Documents/PAMELA_puntos_28_5_2024/levantamiento_malezas/data/Merged_Data_with_Coordinates.csv", delim = ";", show_col_types = FALSE)
    
    # Cambiar los valores de ID a nombres de malezas
    data$Weed <- factor(data$ID, levels = 0:16, labels = weed_names)
    
    output$histogram <- renderPlotly({
      plot_ly(data, x = ~Weed, type = 'histogram',
              marker = list(color = 'white', line = list(color = 'black', width = 1)),
              xbins = list(size = 1)) %>%
        layout(title = "Histograma de IDs de Malezas",
               xaxis = list(title = "Maleza", tickangle = -90),
               yaxis = list(title = "Frecuencia"),
               bargap = 0.1)
    })
  })
  
  output$downloadRaster <- downloadHandler(
    filename = function() {
      paste0("interpolacion_", input$interpolation, ".tif")
    },
    content = function(file) {
      writeRaster(rasterData(), file)
    }
  )
}

shinyApp(ui = ui, server = server)
