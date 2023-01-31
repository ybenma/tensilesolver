library(shiny)
library(tidyverse)

solve_tensile <- function(raw_tt, width, side, length, error, YM_time) {
    require(tidyverse)
    require(plotly)
    if (class(raw_tt) != 'data.frame') {
        stop('raw_tt must be a data frame')
    }
    null_area = width*side*10^(-6)
    
    colnames(raw_tt) = c('Force', 'Distance', 'Time')
    
    ttdf <- raw_tt %>% mutate(Area = (length/(length+Distance))*null_area) %>% 
        mutate(True_stress = Force/Area/1000) %>% 
        mutate(True_strain = log((Distance+length)/length))
    
    # Fracture calculation
    Fracture_info <- ttdf[which.max(ttdf$True_stress),] %>% as_vector()
    Fracture_output <- Fracture_info[c(3,5:6)] %>% as_tibble_row()
    colnames(Fracture_output) <- c("Fracture time (s)", "Fracture stress (Pa)", "Fracture strain")
    F_time <- Fracture_info[3]
    
    # Young's modulus calculation
    YMstress <- ttdf %>% filter(Time <= YM_time)  %>% select(True_stress)
    YMstrain <- ttdf %>% filter(Time <= YM_time)  %>% select(True_strain)
    YM <- lm(YMstress$True_stress~YMstrain$True_strain)$coefficients[2]
    YM_intercept <- lm(YMstress$True_stress~YMstrain$True_strain)$coefficients[1]
    YM_info <- as_tibble_row(c(a = YM, b = YM_intercept))
    colnames(YM_info) <- c("Young's modulus (Pa)", "Intercept")
    
    # Transition point calculation
    ttdf <- ttdf %>%  
        mutate(Offset_stress = True_stress * error/100) %>% 
        mutate(stress_prop = True_strain*YM+YM_intercept) %>% 
        mutate(diff_stress = abs(True_stress - stress_prop)-Offset_stress) %>% 
        as_tibble()
    
    
    transition_output <- ttdf %>% filter(Time >= YM_time & Time <= F_time) %>% 
        mutate(diff_stress = abs(diff_stress)) %>% 
        arrange(diff_stress) %>% 
        select(c(True_stress, True_strain)) %>% slice(1)
    
    colnames(transition_output) <- c('Transition Stress (Pa)', 'Transition Strain')
    
    
    toughness <- ttdf %>% mutate(dif_stress = True_stress - lag(True_stress)) %>% 
        mutate(dif2_stress = dif_stress - lag(dif_stress)) %>% 
        arrange(desc(dif2_stress)) %>% 
        slice(1) %>% select(Time, True_stress, True_strain)
    
    toughness_output <- ttdf %>% filter(Time <= toughness$Time) %>% 
        summarise(toughness = round(sum(True_stress)/1000, 2)) %>% 
        bind_cols(toughness)
    colnames(toughness_output)[1] <- c('Toughness (kPa)')
    
    results <- list(Fracture_output, YM_info, transition_output, 
                    toughness_output, ttdf)
    names(results) <- c('Fracture_info', 'Youngs_modulus', 
                        'Transition','Toughness', 'Processed_data')
    return(results)
    
}


ui <- pageWithSidebar(
    headerPanel('Tensile Solver'),
    sidebarPanel(
        fileInput("file1", "Choose a CSV File", accept = ".csv"),
        
        checkboxInput("header", "Header", TRUE), 
        numericInput(inputId = "width",
                     label = "Dogbone width (mm)",
                     value = 2.7),
        
        numericInput(inputId = "side",
                     label = "Dogbone side (mm)",
                     value = 4.15),
        
        numericInput(inputId = "length",
                     label = "Dogbone length (mm)",
                     value = 8.5),
        
        sliderInput(inputId = "YM_time",
                    label = "Young's modulus cutoff (s)",
                    min = 0.5,
                    max = 5,
                    value = 1.5,
                    step = 0.5),
        
        sliderInput(inputId = "transition_error",
                    label = "Transition Error (%)",
                    min = 1,
                    max = 5,
                    value = 2)
    ),
    mainPanel(
        plotOutput('tensile_plot'),
        tableOutput("Fracture_info"),
        tableOutput("YM_info")
    )
)

server <- function(input, output, session) {
    
    
    tensile_file <- reactive({
        file <- input$file1
        ext <- tools::file_ext(file$datapath)
        
        req(file)
        validate(need(ext == "csv", "Please upload a csv file"))
        
        read.csv(file$datapath, skip = 3, header = T,
                 col.names = c('Force', 'Distance', 'Time'))
    })
    
    all_results <- reactive({
        solve_tensile(tensile_file(), input$width, 
                      input$side, input$length, input$transition_error, 
                      input$YM_time)
    })
    
    
    output$Fracture_info <- renderTable(bind_cols(all_results()$Fracture_info[,-1], 
                                                  all_results()$Toughness[,1]))
    
    output$YM_info <- renderTable(bind_cols(all_results()$Youngs_modulus[,1],
                                            all_results()$Transition))
    
    #output$Toughness <- renderTable(all_results()$Toughness[,1])
    
    output$tensile_plot <- renderPlot({
        all_results()$Processed_data %>% 
            filter(Time <= all_results()$Fracture_info$`Fracture time (s)` + 20) %>% 
            ggplot(aes(True_strain, True_stress))+
            geom_line(lwd = 1.2, aes(color = 'Tensile Test'))+
            geom_point(aes(all_results()$Transition$`Transition Strain`, 
                           all_results()$Transition$`Transition Stress (Pa)`, 
                           color = 'Transition Point'), 
                       color = 'darkred',
                       size = 2)+
            geom_abline(slope = all_results()$Youngs_modulus$`Young's modulus (Pa)`, 
                        intercept = all_results()$Youngs_modulus$Intercept, 
                        lty = 2, aes(color = "Young's Modulus Trend"),
                        color = 'blue')+
            scale_color_manual(name=' ',
                               breaks=c('Tensile test', 'Transition point', "Young's modulus"),
                               values=c('Tensile test'='black', 
                                        'Transition point'='darkred', 
                                        "Young's modulus"='blue'))+
            theme_bw(base_size = 16)+
            xlab('True Strain')+
            ylab('True Stress')
    })
    
}


# Run the application 
shinyApp(ui = ui, server = server)
