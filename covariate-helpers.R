###
## helper functions related to wrangling the covariates, such as conversion to
## dummy matrix format
###

## wrapper to get cov matrix with consistent reference levels
## if as_tibble true, will keep geospatial identifiers and infestation status
## TODO: maybe center and scale continuous covariates?
make_dummy_mat <- function(df_org,
                           inter=TRUE, # include interecept col?
                           as_tibble=FALSE) {
    require(tidyverse)
    require(caret)
    
    dat_disc <- df_org %>%
        select(bed_hygiene:infestation) %>%
        select(-c(num_humans:num_pigs)) # both pigs and cats had enough lvls to remain cont
    
    dat_cont <- df_org %>%
        select(num_humans:num_pigs)
    
    dat_disc %>%
        select_if(is.factor) %>%
        map(levels)
    
    ## relevel/collapse each factor to remove any low level counts
    dat_disc$bed_hygiene <- dat_disc$bed_hygiene %>%
        relevel('buena')
    
    dat_disc$bird_nests_inside <- dat_disc$bird_nests_inside %>%
        relevel('no')
    
    dat_disc$chicken_coop_location <- dat_disc$chicken_coop_location %>%
        relevel('no')
    
    dat_disc$bedroom_clutter <- dat_disc$bedroom_clutter %>%
        relevel('no acumula')
    
    dat_disc$dark_main_room <- dat_disc$dark_main_room %>%
        relevel('no')
    
    dat_disc$education_level <- dat_disc$education_level %>%
        relevel('ninguno')
    
    dat_disc$floor_type <- dat_disc$floor_type %>%
        relevel('ladrillo')
    
    dat_disc$construction_pile_type <- dat_disc$construction_pile_type %>%
        relevel('ninguno')
    
    dat_disc$firewood_location <- dat_disc$firewood_location %>%
        relevel('ninguno')
    
    dat_disc$grain_storage_in_house <- dat_disc$grain_storage_in_house %>%
        relevel('no')

    if (any(levels(dat_disc$house_age) == 'menos de 1 año')) {
        dat_disc$house_age <- dat_disc$house_age %>%
            relevel('menos de 1 año')
    }
    
    dat_disc$house_hygiene <- dat_disc$house_hygiene %>%
        relevel('buena')
    
    ## TODO: for now, try add'l lumping?
    dat_disc$kitchen_location <- dat_disc$kitchen_location %>%
        fct_collapse(shared=c('no tiene', 'cocina compartida'))
    
    dat_disc$land_for_agriculture <- dat_disc$land_for_agriculture %>%
        relevel('no')
    
    dat_disc$material_house_wall <- dat_disc$material_house_wall %>%
        relevel('brick_block_oth')
    
    dat_disc$material_roof <- dat_disc$material_roof %>%
        relevel('alum_cement')
    
    dat_disc$ownership <- dat_disc$ownership %>%
        relevel('owned')
    
    dat_disc$sign_animals <- dat_disc$sign_animals %>%
        relevel('FALSE')
    
    dat_disc$sign_rats <- dat_disc$sign_rats %>%
        relevel('FALSE')
    
    dat_disc$windows_in_bedroom <- dat_disc$windows_in_bedroom %>%
        relevel('FALSE')
    
    dat_disc$condition_house_wall <- dat_disc$condition_house_wall %>%
        relevel('buen estado')
    
    dat_disc$condition_bedroom_wall <- dat_disc$condition_bedroom_wall %>%
        relevel('buen estado')
    
    if (inter) {
        cov_matrix <- matrix(rep(1, nrow(dat_cont)), ncol=1) %>% # make the intercept
            cbind(predict(dummyVars(
                infestation ~ .,
                data=dat_disc,
                fullRank=TRUE
            ), newdata=dat_disc)) %>%
            cbind(as.matrix(dat_cont))
        colnames(cov_matrix)[1] <- 'inter'
    }
    else {
        cov_matrix <- predict(dummyVars(
            infestation ~ .,
            data=dat_disc,
            fullRank=TRUE
        ), newdata=dat_disc) %>%
            cbind(as.matrix(dat_cont))
    }
    if (as_tibble) {
        cov_matrix <- cov_matrix %>%
            as_tibble(.name_repair='minimal')
        cov_matrix <- bind_cols(df_org[,1:4], cov_matrix, infestation=df_org$infestation)
    }
    return(cov_matrix)
}
