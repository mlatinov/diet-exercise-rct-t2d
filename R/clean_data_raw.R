
#### Function to clean the Raw data ####
clean_data_raw <- function(data_raw){

  data_clean <- data_raw %>%
    ## Rename all the variables that will be used in the analysis
    rename(
      # Identefiers 
      "id" = Sr.,
      "treatment" = PIP..1.,

      # Demografics
      "sex" = Sex..F.0.,
      "age" = Age,
      "marital" = Marital.status..Single.0...Marries.1,

      # anthropometrics 
      "height_pre"  = HeighPre,
      "height_post" = HeighPost,
      "weight_pre"  = WtPre,
      "weight_post" = WtPost,
      "bmi_pre"     = BMIPre,
      "bmi_post"    = BMIPost,
      "waist_pre"   = WaistPre,
      "waist_post"  = WaistPost,
      "hip_pre"     = HipPre,
      "hip_post"    = HipPost,

      # Diet intake 
      "carb_pre"    = Pre_Carb,
      "carb_post"   = Post_Carb,
      "energy_pre"  = Pre_energy,
      "energy_post" = Post_energy,

      # Diet Quality 
      "diet_score_pre"  = Pre.Diet,
      "diet_score_post" = Post.Diet,

      # Adherence
      "diet_adherence"     = Adhere..1,
      "exercise_adherence" = EXR_ADH,
      
      # Exercise frequency
      "exercise_type_1_freq_pre"  = Pre.Ex_1,
      "exercise_type_1_freq_post" = Post.Ex_1,
      "exercise_type_2_freq_pre"  = Pre.Ex_2,
      "exercise_type_2_freq_post" = Post.Ex_2,
      "exercise_post_total"       = Post.Exr,

      # Exercise intensity
      "exercise_intensity_pre"  = Pre1.low.2..moderate.3..high,
      "exercise_intensity_post" = Post1.low.2..moderate.3..high,

      # Exersize Duration 
      "exercise_duration_pre"  = Pre1....15.minutes...........2..30.minutes.......................3...60.minutes,
      "exercise_duration_post" = Post1....15.minutes...........2..30.minutes.......................3...60.minutes,
      
      # Self Care 
      "self_care_score_pre"  = Pre_r.SDSCA,
      "self_care_score_post" = Post_r.SDSCA,

      # Knowledge
      "q_score_pre"  = Qpre,
      "q_score_post" = Qpost,
      
      # Barriers 
      "lim_time_food_prep"  = Limited.time.to.prepare.healthy.food,
      "diff_count_calories" = Difficulty.of.counting.carbohydrates,
      "lack_motivation"     = Lack.of.motivation,
      "freq_social_interactions" = Frequent.social.invitations,
      "fast_food_delivery"       = Fast.food.delivery,
      "home_load"                = Children.and.home.workload
    ) %>%
    
    # Drop Derived Columns
    select(-any_of(c(
      "D_Heigh", "D_WT", "D_BMI", "D_Waist", "D_Hip", "D_W.H",
      "Waist.Hip_Pre", "Waist.Hip_Post","Pre.Ex_Total"
  ))) %>%

    # Convert Categorical Features to factors 
    mutate(
      # Demografics
      sex     = factor(sex, levels = c(0, 1), labels = c("female", "male")),
      marital = factor(marital, levels = c(0, 1), labels = c("single", "married")),
      treatment       = factor(treatment, levels = c(0, 1), labels = c("control", "intervention")),
      # Diet 
      diet_adherence  = factor(diet_adherence, levels = c(0, 1), labels = c("no", "yes")),
      # Exercise 
      exercise_intensity_pre  = factor(exercise_intensity_pre,
        levels = 0:3, 
        labels = c(NA, "low", "moderate", "high"),
        ordered = TRUE
      ),
      exercise_intensity_post = factor(exercise_intensity_post,
        levels = 0:3,
        labels = c(NA,"low", "moderate", "high"),
        ordered = TRUE
      ),
      exercise_duration_pre   = factor(exercise_duration_pre,
        levels = 0:3, 
        labels = c(NA,"15min", "30min", "60min"),
        ordered = TRUE
      ),
      exercise_duration_post = factor(exercise_duration_post,
        levels = 0:3,
        labels = c(NA,"15min", "30min", "60min"),
        ordered = TRUE
      ),
      exercise_adherence = factor(exercise_adherence, levels = c(0, 1), labels = c("no", "yes")),

      # Barriers
      lim_time_food_prep = factor(lim_time_food_prep, levels = c(0, 1), labels = c("no", "yes")),
      lack_motivation    = factor(lack_motivation,    levels = c(0, 1), labels = c("no", "yes")),
      fast_food_delivery = factor(fast_food_delivery, levels = c(0, 1), labels = c("no", "yes")),
      home_load          = factor(home_load,          levels = c(0, 1), labels = c("no", "yes")), 
      freq_social_interactions = factor(freq_social_interactions, levels = c(0, 1), labels = c("no", "yes"))
    )
}

