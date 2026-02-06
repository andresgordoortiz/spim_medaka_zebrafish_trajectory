library(dplyr)
library(readr)
library(plotly)
library(shiny)
library(pracma)

options(shiny.maxRequestSize = 500 * 1024^2)

# -----------------------
# CORE FUNCTIONS
# -----------------------

load_trackmate_data <- function(spots_file, tracks_file) {
  spots_raw <- read_csv(spots_file, show_col_types = FALSE)[-c(1:3),]
  tracks_raw <- read_csv(tracks_file, show_col_types = FALSE)[-c(1:3),]
  
  spots <- spots_raw %>%
    mutate(across(c(ID, TRACK_ID, POSITION_X, POSITION_Y, POSITION_Z, FRAME, QUALITY), ~as.numeric(.))) %>%
    filter(!is.na(TRACK_ID))
  
  tracks <- tracks_raw %>% mutate(across(-LABEL, ~as.numeric(.)))
  
  list(spots = spots, tracks = tracks)
}

detect_animal_pole_robust <- function(spots) {
  early_spots <- spots %>% filter(FRAME == min(FRAME, na.rm = TRUE))
  if(nrow(early_spots) == 0) early_spots <- spots 
  c(mean(early_spots$POSITION_X, na.rm = TRUE),
    mean(early_spots$POSITION_Y, na.rm = TRUE),
    mean(early_spots$POSITION_Z, na.rm = TRUE))
}

get_rotation_matrix <- function(u, v) {
  u <- u / sqrt(sum(u^2)); v <- v / sqrt(sum(v^2))
  if(sum(abs(u - v)) < 1e-6) return(diag(3))
  if(sum(abs(u + v)) < 1e-6) return(diag(c(1, -1, -1)))
  axis <- pracma::cross(u, v)
  axis <- axis / sqrt(sum(axis^2))
  angle <- acos(max(-1, min(1, sum(u * v))))
  K <- matrix(c(0, -axis[3], axis[2],
                axis[3], 0, -axis[1],
                -axis[2], axis[1], 0), 3, 3, byrow = TRUE)
  diag(3) + sin(angle) * K + (1 - cos(angle)) * (K %*% K)
}

apply_transform <- function(spots, ap, dorsal) {
  mid_x <- mean(spots$POSITION_X, na.rm = TRUE)
  mid_y <- mean(spots$POSITION_Y, na.rm = TRUE)
  mid_z <- mean(spots$POSITION_Z, na.rm = TRUE)
  
  mat <- as.matrix(spots[, c("POSITION_X", "POSITION_Y", "POSITION_Z")])
  mat_centered <- sweep(mat, 2, c(mid_x, mid_y, mid_z), FUN = "-")
  
  # Map animal pole to +Y
  ap_c <- ap - c(mid_x, mid_y, mid_z)
  dorsal_c <- dorsal - c(mid_x, mid_y, mid_z)
  
  R1 <- get_rotation_matrix(ap_c, c(0, 1, 0))
  mat_r1 <- t(R1 %*% t(mat_centered))
  
  # Rotate dorsal to +X while preserving +Y
  dorsal_r1 <- as.numeric(R1 %*% dorsal_c)
  dorsal_xz <- c(dorsal_r1[1], 0, dorsal_r1[3])
  if(sum(abs(dorsal_xz)) < 1e-8) {
    R2 <- diag(3)
  } else {
    R2 <- get_rotation_matrix(dorsal_xz, c(1, 0, 0))
  }
  mat_final <- t(R2 %*% t(mat_r1))
  
  spots$POSITION_X <- mat_final[, 1]
  spots$POSITION_Y <- mat_final[, 2]
  spots$POSITION_Z <- mat_final[, 3]
  
  list(data = spots, R1 = R1, R2 = R2, center = c(mid_x, mid_y, mid_z))
}

transform_tracks <- function(tracks_df, center, R1, R2) {
  if(is.null(tracks_df) || nrow(tracks_df) == 0) return(tracks_df)
  cols <- colnames(tracks_df)
  candidate <- grep("_(?:POSITION_)?[XYZ]$", cols, perl = TRUE, value = TRUE)
  if(length(candidate) == 0) return(tracks_df)
  stems <- unique(sub("_(?:POSITION_)?[XYZ]$", "", candidate, perl = TRUE))
  transformed_any <- FALSE
  
  for(stem in stems) {
    xname_opts <- c(paste0(stem, "_POSITION_X"), paste0(stem, "_X"))
    yname_opts <- c(paste0(stem, "_POSITION_Y"), paste0(stem, "_Y"))
    zname_opts <- c(paste0(stem, "_POSITION_Z"), paste0(stem, "_Z"))
    
    xname <- xname_opts[which(xname_opts %in% cols)[1]]
    yname <- yname_opts[which(yname_opts %in% cols)[1]]
    zname <- zname_opts[which(zname_opts %in% cols)[1]]
    if(is.na(xname) || is.na(yname) || is.na(zname)) next
    
    mat_raw <- tracks_df[, c(xname, yname, zname)]
    mat <- tryCatch({
      apply(mat_raw, 2, function(col) as.numeric(as.character(col)))
    }, warning = function(w) {
      apply(mat_raw, 2, function(col) suppressWarnings(as.numeric(col)))
    }, error = function(e) NULL)
    
    if(is.null(mat) || nrow(mat) == 0) next
    if(is.null(dim(mat))) mat <- matrix(mat, ncol = 3)
    if(ncol(mat) != 3) next
    
    mat_centered <- sweep(mat, 2, center, FUN = "-")
    mat_r1 <- t(R1 %*% t(mat_centered))
    mat_final <- t(R2 %*% t(mat_r1))
    
    tracks_df[[xname]] <- mat_final[, 1]
    tracks_df[[yname]] <- mat_final[, 2]
    tracks_df[[zname]] <- mat_final[, 3]
    
    transformed_any <- TRUE
  }
  if(!transformed_any) return(tracks_df)
  tracks_df
}

# -----------------------
# UI
# -----------------------

ui <- fluidPage(
  titlePanel("Interactive Embryo Orientation Tool"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("1. Upload Data"),
      fileInput("spots_upload", "Upload Spots .csv", accept = ".csv"),
      fileInput("tracks_upload", "Upload Tracks .csv", accept = ".csv"),
      hr(),
      
      h4("2. Pick Landmarks"),
      p("Click 'Set...' then click a dot on the 3D plot."),
      uiOutput("ap_button"),
      textOutput("val_ap"),
      br(),
      uiOutput("dorsal_button"),
      textOutput("val_dorsal"),
      hr(),
      
      h4("Display options"),
      sliderInput("plot_n_spots", "Number of spots to display", min = 100, max = 50000, value = 5000, step = 100),
      sliderInput("marker_size", "Marker size", min = 1, max = 12, value = 2, step = 1),
      hr(),
      
      h4("Standardized coloring (choose axis)"),
      checkboxGroupInput("color_axes", label = NULL,
                         choices = c("Z (depth)", "Dorsal–Ventral (X)", "Animal–Vegetal (Y)")),
      helpText("If multiple boxes are checked, the first selected choice is used for coloring."),
      hr(),
      
      h4("3. Run"),
      actionButton("btn_process", "Standardize Orientation", class = "btn-success", width = "100%"),
      br(), br(),
      downloadButton("download_data", "Download Result (ZIP)", style = "width:100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Interactive View", plotlyOutput("plot_interactive", height = "700px")),
        tabPanel("Standardized Result", plotlyOutput("plot_result", height = "700px"))
      )
    )
  )
)

# -----------------------
# SERVER
# -----------------------

server <- function(input, output, session) {
  vals <- reactiveValues(
    spots = NULL, tracks = NULL,
    animal_pole = NULL, dorsal_landmark = NULL,
    pick_mode = "none",
    standardized_data = NULL, standardized_tracks = NULL,
    camera = NULL
  )
  
  observe({
    req(input$spots_upload, input$tracks_upload)
    tryCatch({
      data <- load_trackmate_data(input$spots_upload$datapath, input$tracks_upload$datapath)
      vals$spots <- data$spots
      vals$tracks <- data$tracks
      vals$animal_pole <- detect_animal_pole_robust(vals$spots)
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) showNotification(paste("Error loading data:", e$message), type = "error"))
  })
  
  output$ap_button <- renderUI({
    if (isolate(vals$pick_mode) == "ap") actionButton("btn_pick_ap", "Set Animal Pole (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_ap", "Set Animal Pole", width = "100%", class = "btn-primary")
  })
  output$dorsal_button <- renderUI({
    if (isolate(vals$pick_mode) == "dorsal") actionButton("btn_pick_dorsal", "Set Dorsal (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_dorsal", "Set Dorsal", width = "100%", class = "btn-info")
  })
  
  observeEvent(input$btn_pick_ap, {
    vals$pick_mode <- ifelse(vals$pick_mode == "ap", "none", "ap")
    # update UI labels (these do not affect the plot because plot reads pick_mode via isolate)
    output$ap_button <- renderUI({
      if (vals$pick_mode == "ap") actionButton("btn_pick_ap", "Set Animal Pole (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_ap", "Set Animal Pole", width = "100%", class = "btn-primary")
    })
    output$dorsal_button <- renderUI({
      if (vals$pick_mode == "dorsal") actionButton("btn_pick_dorsal", "Set Dorsal (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_dorsal", "Set Dorsal", width = "100%", class = "btn-info")
    })
  })
  observeEvent(input$btn_pick_dorsal, {
    vals$pick_mode <- ifelse(vals$pick_mode == "dorsal", "none", "dorsal")
    output$ap_button <- renderUI({
      if (vals$pick_mode == "ap") actionButton("btn_pick_ap", "Set Animal Pole (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_ap", "Set Animal Pole", width = "100%", class = "btn-primary")
    })
    output$dorsal_button <- renderUI({
      if (vals$pick_mode == "dorsal") actionButton("btn_pick_dorsal", "Set Dorsal (SELECTING)", width = "100%", class = "btn-warning") else actionButton("btn_pick_dorsal", "Set Dorsal", width = "100%", class = "btn-info")
    })
  })
  
  # handle clicks for picking landmarks
  observeEvent(event_data("plotly_click", source = "raw_plot"), {
    click_data <- event_data("plotly_click", source = "raw_plot")
    if(is.null(click_data)) return()
    coords <- c(click_data$x, click_data$y, click_data$z)
    if(vals$pick_mode == "ap") {
      vals$animal_pole <- coords
      vals$pick_mode <- "none"
      output$ap_button <- renderUI({ actionButton("btn_pick_ap", "Set Animal Pole", width = "100%", class = "btn-primary") })
      output$dorsal_button <- renderUI({ actionButton("btn_pick_dorsal", "Set Dorsal", width = "100%", class = "btn-info") })
      showNotification("Animal Pole Updated")
    } else if(vals$pick_mode == "dorsal") {
      vals$dorsal_landmark <- coords
      vals$pick_mode <- "none"
      output$ap_button <- renderUI({ actionButton("btn_pick_ap", "Set Animal Pole", width = "100%", class = "btn-primary") })
      output$dorsal_button <- renderUI({ actionButton("btn_pick_dorsal", "Set Dorsal", width = "100%", class = "btn-info") })
      showNotification("Dorsal Point Updated")
    }
  })
  
  # capture camera on relayout
  observeEvent(event_data("plotly_relayout", source = "raw_plot"), {
    rel <- event_data("plotly_relayout", source = "raw_plot")
    if(is.null(rel)) return()
    names_rel <- names(rel)
    ex <- if("scene.camera.eye.x" %in% names_rel) rel[["scene.camera.eye.x"]] else if("camera.eye.x" %in% names_rel) rel[["camera.eye.x"]] else NULL
    ey <- if("scene.camera.eye.y" %in% names_rel) rel[["scene.camera.eye.y"]] else if("camera.eye.y" %in% names_rel) rel[["camera.eye.y"]] else NULL
    ez <- if("scene.camera.eye.z" %in% names_rel) rel[["scene.camera.eye.z"]] else if("camera.eye.z" %in% names_rel) rel[["camera.eye.z"]] else NULL
    if(!is.null(ex) && !is.null(ey) && !is.null(ez)) vals$camera <- list(eye = list(x = ex, y = ey, z = ez))
  })
  
  # interactive plot — IMPORTANT: use isolate(vals$pick_mode) so toggling pick-mode DOES NOT re-render the plot
  output$plot_interactive <- renderPlotly({
    req(vals$spots)
    nshow <- min(nrow(vals$spots), max(1, as.integer(input$plot_n_spots)))
    plot_data <- vals$spots %>% sample_n(nshow)
    msize <- max(1, as.numeric(input$marker_size))
    
    p <- plot_ly(source = "raw_plot") %>%
      add_trace(data = plot_data,
                x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                type = "scatter3d", mode = "markers",
                marker = list(size = msize, color = plot_data$POSITION_Z, colorscale = "Viridis", showscale = TRUE),
                text = ~ID, hoverinfo = "text", showlegend = FALSE)
    
    if(!is.null(vals$animal_pole)) {
      p <- p %>% add_trace(x = vals$animal_pole[1], y = vals$animal_pole[2], z = vals$animal_pole[3],
                           type = "scatter3d", mode = "markers",
                           marker = list(size = max(4, msize + 6), color = "red", symbol = "diamond"),
                           name = "Animal Pole")
    }
    if(!is.null(vals$dorsal_landmark)) {
      p <- p %>% add_trace(x = vals$dorsal_landmark[1], y = vals$dorsal_landmark[2], z = vals$dorsal_landmark[3],
                           type = "scatter3d", mode = "markers",
                           marker = list(size = max(4, msize + 6), color = "cyan", symbol = "square"),
                           name = "Dorsal")
    }
    
    # use isolate() so this read does NOT create a reactive dependency
    pick_mode_val <- isolate(vals$pick_mode)
    title_str <- "Original Data"
    if(pick_mode_val == "ap") title_str <- "CLICK A DOT TO SET ANIMAL POLE"
    if(pick_mode_val == "dorsal") title_str <- "CLICK A DOT TO SET DORSAL SIDE"
    
    scene_list <- list(aspectmode = "data")
    if(!is.null(vals$camera)) scene_list$camera <- vals$camera
    
    p %>% layout(title = title_str, scene = scene_list) %>% event_register("plotly_click")
  })
  
  output$val_ap <- renderText({ if(is.null(vals$animal_pole)) "Not set" else paste("AP set at Y=", round(vals$animal_pole[2], 1)) })
  output$val_dorsal <- renderText({ if(is.null(vals$dorsal_landmark)) "Not set" else paste("Dorsal set at X=", round(vals$dorsal_landmark[1], 1)) })
  
  observeEvent(input$btn_process, {
    req(vals$spots, vals$animal_pole, vals$dorsal_landmark)
    withProgress(message = "Rotating...", {
      res <- apply_transform(vals$spots, vals$animal_pole, vals$dorsal_landmark)
      vals$standardized_data <- res$data
      vals$standardized_tracks <- transform_tracks(vals$tracks, res$center, res$R1, res$R2)
      if(!is.null(vals$standardized_tracks)) {
        candidate_after <- grep("_(?:POSITION_)?[XYZ]$", colnames(vals$standardized_tracks), perl = TRUE, value = TRUE)
        if(length(candidate_after) == 0) showNotification("Standardized spots saved. No coordinate triplets found in tracks to transform.", type = "warning") else showNotification("Standardized spots + tracks saved.", type = "message")
      } else showNotification("Standardized spots saved. Tracks not available.", type = "warning")
    })
  })
  
  output$plot_result <- renderPlotly({
    req(vals$standardized_data)
    nshow <- min(nrow(vals$standardized_data), max(1, as.integer(input$plot_n_spots)))
    dat <- vals$standardized_data %>% sample_n(nshow)
    msize <- max(1, as.numeric(input$marker_size))
    
    # Determine color axis safely
    chosen <- NULL
    if(!is.null(input$color_axes) && length(input$color_axes) > 0) chosen <- input$color_axes[[1]]
    if(!is.null(chosen)) {
      color_var <- switch(as.character(chosen),
                          "Z (depth)" = "POSITION_Z",
                          "Dorsal–Ventral (X)" = "POSITION_X",
                          "Animal–Vegetal (Y)" = "POSITION_Y",
                          NULL)
    } else {
      color_var <- NULL
    }
    subtitle_txt <- if(is.null(color_var)) "No color gradient" else paste("Color by:", chosen)
    
    scene_list <- list(aspectmode = 'data',
                       xaxis = list(title = "Dorsal (+X)"),
                       yaxis = list(title = "AP (+Y)"),
                       zaxis = list(title = "Depth (+Z)"))
    if(!is.null(vals$camera)) scene_list$camera <- vals$camera
    
    # Use annotation in title to show subtitle safely
    title_html <- if(is.null(color_var)) {
      list(text = "Standardized Data")
    } else {
      list(text = paste0("Standardized Data<br><sup>", subtitle_txt, "</sup>"))
    }
    
    if(is.null(color_var)) {
      p <- plot_ly(dat, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                   type = "scatter3d", mode = "markers",
                   marker = list(size = msize, color = "grey40"),
                   text = ~ID, hoverinfo = "text")
    } else {
      # safe subset of the color vector (do not call dat[[NULL]])
      color_vec <- dat[[color_var]]
      p <- plot_ly(dat, x = ~POSITION_X, y = ~POSITION_Y, z = ~POSITION_Z,
                   type = "scatter3d", mode = "markers",
                   marker = list(size = msize, color = color_vec, colorscale = "Plasma", showscale = TRUE),
                   text = ~ID, hoverinfo = "text")
    }
    
    p %>% layout(title = title_html, scene = scene_list)
  })
  
  output$download_data <- downloadHandler(
    filename = function() { paste0("oriented_data_", Sys.Date(), ".zip") },
    content = function(file) {
      req(vals$standardized_data)
      tmpdir <- tempdir()
      spots_file <- file.path(tmpdir, "oriented_spots.csv")
      write_csv(vals$standardized_data, spots_file)
      files_to_zip <- c("oriented_spots.csv")
      if(!is.null(vals$standardized_tracks)) {
        tracks_file <- file.path(tmpdir, "oriented_tracks.csv")
        write_csv(vals$standardized_tracks, tracks_file)
        files_to_zip <- c(files_to_zip, "oriented_tracks.csv")
      }
      oldwd <- setwd(tmpdir); on.exit(setwd(oldwd), add = TRUE)
      utils::zip(zipfile = file, files = files_to_zip)
    }
  )
}

shinyApp(ui, server)
