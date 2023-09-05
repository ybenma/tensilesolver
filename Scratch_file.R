getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

solve_tensile <- function(raw_tt, width, side, length, error, YM_time) {
   require(tidyverse)
   if (class(raw_tt) != 'data.frame') {
      stop('raw_tt must be a data frame')
   }
   null_area = width*side*10^(-6) #m2
   
   colnames(raw_tt) = c('Force', 'Distance', 'Time')
   
   mode <- round(raw_tt$Force, 3) %>% getmode()
   #cutoff_time <- raw_tt %>% mutate(norm_force = abs(Force-mode)) %>%
   #   filter(Time >= 0.5) %>% 
   #   arrange(norm_force) %>% 
   #   slice(1)
   cutoff_time <- raw_tt %>% filter(Force > mode) %>% 
      mutate(norm_force = abs(Force-mode)) %>%
      arrange(norm_force) %>% slice(1)
   
   process_df <- raw_tt %>% filter(Time < cutoff_time$Time+5)
   
   ttdf <- process_df %>% mutate(Area = (length/(length+Distance))*null_area) %>% 
      mutate(True_stress = Force/Area/1000) %>% #kPa
      mutate(True_strain = log((Distance+length)/length))
   
   # Fracture calculation
   Fracture_info <- ttdf[which.max(ttdf$True_stress),] %>% as_vector()
   Fracture_output <- Fracture_info[c(3,5:6)] %>% as_tibble_row()
   colnames(Fracture_output) <- c("Fracture time (s)", "Fracture stress (kPa)", "Fracture strain")
   F_time <- Fracture_info[3]
   
   # Young's modulus calculation
   YMstress <- ttdf %>% filter(Time <= YM_time)  %>% select(True_stress)
   YMstrain <- ttdf %>% filter(Time <= YM_time)  %>% select(True_strain)
   YM <- lm(YMstress$True_stress~YMstrain$True_strain)$coefficients[2]
   YM_intercept <- lm(YMstress$True_stress~YMstrain$True_strain)$coefficients[1]
   YM_info <- as_tibble_row(c(a = YM, b = YM_intercept))
   colnames(YM_info) <- c("Young's modulus (kPa)", "Intercept")
   
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
   
   colnames(transition_output) <- c('Transition Stress (kPa)', 'Transition Strain')
   
   #toughness calculation
   toughness <- ttdf %>% filter(Time == cutoff_time$Time)

   #toughness <- ttdf %>% 
   #   mutate(dif_stress = True_stress - lag(True_stress, 10)) %>% 
   #   mutate(dif2_stress = dif_stress - lag(dif_stress, 10)) %>% 
   #   arrange(desc(dif2_stress)) %>% 
   #   slice(1) %>% select(Time, True_stress, True_strain)
   
   toughness_output <- ttdf %>% filter(Time <= toughness$Time) %>% 
      summarise(toughness = round(sum(True_stress)/1000, 4)) %>% 
      bind_cols(toughness)
   
   
   colnames(toughness_output)[1] <- c('Toughness (MPa)')
   
   results <- list(Fracture_output, YM_info, transition_output, 
                   toughness_output, ttdf)
   names(results) <- c('Fracture_info', 'Youngs_modulus', 
                       'Transition','Toughness', 'Processed_data')
   return(results)
   
}



file <- read.csv('test_1_1_SPI_per03.csv', skip = 3, header = T,
                 col.names = c('Force', 'Distance', 'Time'))

mode <- round(file$Force, 3) %>% getmode()

file %>% mutate(Force = round(Force, 3)) %>% 
   filter(Force > mode) %>% 
   #mutate(norm_force = abs(Force-mode)) %>%
   ggplot(aes(Time, Force))+
   geom_path()

file %>% mutate(Force = round(Force, 3)) %>% 
   filter(Force > mode) %>% 
   #mutate(norm_force = abs(Force-mode)) %>% 
   arrange(Force) %>% slice(4)

results <- solve_tensile(file, 3.19, 5.11, 8.5, 1.5, 2)


results$Processed_data %>% ggplot(aes(Time, Force))+
   geom_path()

file %>% ggplot(aes(Time, Force))+
   geom_path()

             